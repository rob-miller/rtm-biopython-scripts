#!/usr/local/bin/python3
# Copyright (c) 2019 Robert T. Miller.  All rights reserved.

# -*- coding: latin-1 -*-
"""Interconvert PDB internal and external coordinates.

Interconvert PDB Structure data between external (X, Y, Z cartesian)
coordinates and internal (bond length, angle and dihedral angle) measurements.
"""

#
# replicate buildprot with biopython
#

import argparse
import os
import sys
import re
import cProfile
import pstats
import timeit

# print(sys.path)

# import pyximport

# pyximport.install(language_level=3)

# import itertools
# print(sys.path)

import gzip
import warnings

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from io import StringIO

from Bio.PDB.ic_rebuild import (
    report_IC,
    structure_rebuild_test,
    write_PDB,
    IC_duplicate,
    compare_residues,
)

from Bio.PDB.PICIO import write_PIC, read_PIC
from Bio.PDB.internal_coords import IC_Residue, AtomKey, IC_Chain
from Bio.PDB.SCADIO import write_SCAD

PDB_repository_base = None

pdbDirs = ["/Users/rob/data/pdb/", "/media/data/rcsb/pdb", "/Volumes/data/pdb/"]
for d in pdbDirs:
    if os.path.isdir(d):
        PDB_repository_base = d
        break

if PDB_repository_base is not None:
    print("pdb base = ", PDB_repository_base)

print(sys.version)

scale_val = 2
scadCode_val = True
ti_iter = 3


def doCpti(args, cmd, pdb_structure):
    """Profile handler for cmd."""
    if args.ti:
        tr = timeit.repeat(stmt=cmd, repeat=int(ti_iter), number=1, globals=globals(),)
        print(min(tr))

    elif args.cp:
        cmd2 = "for i in range(0, " + str(ti_iter) + "): " + cmd
        cProfile.runctx(
            cmd2, globals(), locals(), "Profile.prof",
        )
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("cumulative").print_stats()


def adj60(ric, angle):
    """Add 60 degrees to angle of residue internal coord."""
    ang = ric.get_angle(angle)
    if ang is not None:
        ric.set_angle(angle, ang + 60.0)


def test66(pdb_structure, alist):
    """change dihedral angles 6x60."""
    for i in range(6):
        for r in pdb_structure.get_residues():
            ric = r.internal_coord
            if ric:
                for a in alist:
                    adj60(ric, a)
        pdb_structure.internal_to_atom_coordinates()


arg_parser = argparse.ArgumentParser(
    description="Interconvert .pic (protein internal coordinates) and "
    ".pdb (protein data bank) files."
)
arg_parser.add_argument(
    "file",
    nargs="*",
    help="a .pdb, .cif or .pic path/filename to read, or a PDB idCode with "
    "optional chain ID to read from {0} as .ent.gz".format(
        (
            PDB_repository_base
            or "[PDB resource not defined - please configure before use]"
        )
    ),
)

arg_parser.add_argument(
    "-f",
    dest="filelist",
    help="a Dunbrack cullPDB pdb ID list to read from {0} as .ent.gz".format(
        (
            PDB_repository_base
            or "[PDB resource not defined - please configure before use]"
        )
    ),
)
arg_parser.add_argument(
    "-skip", dest="skip_count", help="count of pdb ID list entries to skip"
)
arg_parser.add_argument(
    "-limit",
    dest="limit_count",
    help="stop after processing this many pdb ID list entries",
)

arg_parser.add_argument(
    "-chain", dest="sel_chain", help="chain to pick from PDB/mmCIF file"
)

arg_parser.add_argument(
    "-wp", help="write pdb file with .pdb extension", action="store_true"
)
arg_parser.add_argument(
    "-wip",
    help="convert to internal coords and back, write .pdb file",
    action="store_true",
)
arg_parser.add_argument(
    "-wpi", help="convert to atom coords and back, write .pic file", action="store_true"
)
arg_parser.add_argument(
    "-wc", help="write MMCIF file with .cif extension", action="store_true"
)
arg_parser.add_argument(
    "-wi",
    help="write internal coordinates file with .pic extension",
    action="store_true",
)
arg_parser.add_argument(
    "-ws", help="write OpenSCAD file with .scad extension", action="store_true"
)
arg_parser.add_argument(
    "-scale",
    dest="scale",
    help="OpenSCAD output: units (usually mm) per "
    "angstrom, default " + str(scale_val),
)
arg_parser.add_argument(
    "-nsc", help="OpenSCAD output: no OpenSCAD code, just arrays", action="store_true"
)
arg_parser.add_argument(
    "-flex", help="OpenSCAD output: rotatable backbone CA bonds", action="store_true"
)
arg_parser.add_argument(
    "-hb", help="OpenSCAD output: magnetic backbone H bonds", action="store_true"
)
arg_parser.add_argument(
    "-maxp",
    dest="maxp",
    help="max N-C peptide bond length for chain breaks, "
    "default " + str(IC_Chain.MaxPeptideBond),
)
arg_parser.add_argument(
    "-backbone", help="OpenSCAD output: skip sidechains", action="store_true"
)
arg_parser.add_argument(
    "-t", help="test conversion pdb/pic to pic/pdb", action="store_true"
)
arg_parser.add_argument(
    "-tv", help="verbose test conversion pdb<>pic", action="store_true"
)
arg_parser.add_argument(
    "-cp", help="cprofile option for -t, -tv, -t6x, -trp", action="store_true"
)
arg_parser.add_argument(
    "-ti", help="timeit option for -t, -tv, -t6x, -trp", action="store_true"
)
arg_parser.add_argument(
    "-t6a",
    help="read pdb, change all dihedrals 6x60 deg & assemble, back to pdb and test",
    action="store_true",
)
arg_parser.add_argument(
    "-t6p",
    help="read pdb, change every psi 6x60 deg & assemble, back to pdb and test",
    action="store_true",
)
arg_parser.add_argument(
    "-t6s",
    help="read pdb, change every sidechain dihedral 6x60 deg & assemble, back to pdb and test",
    action="store_true",
)
arg_parser.add_argument(
    "-t6c",
    help="read pdb, change every chi1 6x60 deg & assemble, back to pdb and test",
    action="store_true",
)
arg_parser.add_argument(
    "-trp", help="read pdb, compute internal coords (profile)", action="store_true"
)
arg_parser.add_argument(
    "-ti_iter",
    dest="ti_iter",
    help="cProfile, timeit iteration count option for -t, -tv, -t6x; default "
    + str(ti_iter),
)
arg_parser.add_argument("-nh", help="ignore hydrogens on PDB read", action="store_true")
arg_parser.add_argument(
    "-d2h", help="swap D (deuterium) for H on PDB read", action="store_true"
)
arg_parser.add_argument(
    "-amide", help="only amide proton, skip other Hs on PDB read", action="store_true"
)
arg_parser.add_argument("-gcb", help="generate GLY C-beta atoms", action="store_true")
arg_parser.add_argument(
    "-rama", help="print psi, phi, omega values", action="store_true"
)
arg_parser.add_argument(
    "-rama2", help="test setting dihedral angles", action="store_true"
)

arg_parser.add_argument(
    "-i", help="test just convert to internal coordinates", action="store_true"
)

arg_parser.add_argument("-p", help="print header info", action="store_true")

args = arg_parser.parse_args()

# print(args)

if args.nh:
    IC_Residue.accept_atoms = IC_Residue.accept_mainchain
if args.amide:
    IC_Residue.accept_atoms = IC_Residue.accept_mainchain + ("H",)
if args.d2h:
    IC_Residue.accept_atoms += IC_Residue.accept_deuteriums
    AtomKey.d2h = True
if args.gcb:
    IC_Residue.gly_Cbeta = True
if args.maxp:
    IC_Chain.MaxPeptideBond = float(args.maxp)
if args.scale:
    scale_val = float(args.scale)
if args.skip_count:
    args.skip_count = int(args.skip_count)
if args.limit_count:
    args.limit_count = int(args.limit_count)

toProcess = args.file
pdbidre = re.compile(r"(^\d(\w\w)\w)(\w)?$")

if args.filelist:
    flist = open(args.filelist, "r")
    for aline in flist:
        fields = aline.split()
        pdbidMatch = pdbidre.match(fields[0])
        if pdbidMatch:
            # print(m.group(1) + ' ' + m.group(2))
            # toProcess.append(PDB_repository_base + m.group(2)
            # + '/pdb' + m.group(1) + '.ent.gz' )
            toProcess.append(pdbidMatch.group(0))


if len(toProcess):
    print(len(toProcess), "entries to process")
else:
    print("no files to process. use '-h' for help")
    sys.exit(0)

PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)
CIF_parser = MMCIFParser(QUIET=True)

fileNo = 1

for target in toProcess:
    if args.skip_count and fileNo <= args.skip_count:
        fileNo += 1
        continue
    if args.limit_count is not None:
        if args.limit_count <= 0:
            # sys.exit(0)
            break
        args.limit_count -= 1

    pdb_input = False
    pic_input = False
    cif_input = False
    pdb_structure = None
    # pdb_chain = None
    prot_id = ""
    outfile = os.path.basename(target)

    pdbidMatch = pdbidre.match(target)
    if pdbidMatch is not None:
        assert PDB_repository_base, "PDB repository base directory missing, "
        "please configure for this host"
        pdb_input = True
        filename = (
            PDB_repository_base
            + pdbidMatch.group(2).lower()
            + "/pdb"
            + pdbidMatch.group(1).lower()
            + ".ent.gz"
        )
        prot_id = pdbidMatch.group(1)
    else:
        filename = target
        pdbidre = re.compile(r"(\d\w\w\w)(\w)?")  # find anywhere in string
        pdbidMatch2 = pdbidre.search(target)
        if pdbidMatch2:
            prot_id = pdbidMatch2.group(0)
        else:
            prot_id = target

    if not pdb_input:
        try:
            pdb_structure = read_PIC(
                gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename
            )
        except FileNotFoundError:  # OSError python2
            # print(print(os.getcwd()), filename, "not found")
            print(os.getcwd(), filename, "not found")
            fileNo += 1
            continue

        if pdb_structure is not None:
            pic_input = True

    if pdb_structure is None:
        pdb_input = True
        try:
            pdb_structure = PDB_parser.get_structure(
                prot_id,
                gzip.open(filename, mode="rt")
                if filename.endswith(".gz")
                else filename,
            )
        except FileNotFoundError:  # OSError python2
            # print(print(os.getcwd()), filename, "not found")
            print(os.getcwd(), filename, "not found")
            fileNo += 1
            continue
        except Exception:
            pass  # continue  # try as cif below

        if pdb_structure is None or pdb_structure.child_dict == {}:
            # try:
            pdb_structure = CIF_parser.get_structure(
                prot_id,
                gzip.open(filename, mode="rt")
                if filename.endswith(".gz")
                else filename,
            )
            # except Exception:
            #    print("unable to open ", filename, " as PDB, MMCIF or PIC file format.")
            #    fileNo += 1
            #    continue

        if args.p:
            for k, v in pdb_structure.header.items():
                print(k, ":", v)

    # get specified chain if given

    if pdbidMatch is not None and pdbidMatch.group(3) is not None:
        # have chain specifier for PDBid
        # if pdb_structure[0][pdbidMatch.group(3)] is not None:
        if pdbidMatch.group(3) in pdb_structure[0]:
            pdb_chain = pdb_structure[0][pdbidMatch.group(3)]
            pdb_structure = pdb_chain
        else:
            print("chain " + pdbidMatch.group(3) + " not found in " + filename)
            continue
    elif args.sel_chain is not None:
        if args.sel_chain in pdb_structure[0]:
            pdb_chain = pdb_structure[0][args.sel_chain]
            pdb_structure = pdb_chain
        else:
            print("chain " + args.sel_chain + " not found in " + filename)
            continue

    if pdb_input:
        if not args.t:
            print(fileNo, "parsed pdb input ID", prot_id, "file:", filename)
        # print('header:', pdb_structure.header.get('head', 'NONE'))
        # print('idcode:', pdb_structure.header.get('idcode', 'NONE'))
        # print('deposition date:', pdb_structure.header.get(
        #    'deposition_date', 'NONE'))
    #    for res in pdb_chain.get_residues():   # pdb_structure.get_residues():
    #        print(res.get_full_id(), res.resname,
    #              'disordered' if res.disordered else '')
    else:
        print(fileNo, "parsed pic input ", filename)
        report_IC(pdb_structure, verbose=True)

    # print(pdb_structure.header['idcode'], pdb_chain.id, ':',
    #      pdb_structure.header['head'])

    if args.wp:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates()
        write_PDB(pdb_structure, outfile + ".pdb")
        print("wrote pdb output for", outfile)

    if args.wip:
        if pdb_input:
            pdb2 = IC_duplicate(pdb_structure)
            pdb2.internal_to_atom_coordinates()
            write_PDB(pdb2, outfile + ".pdb")
            print("wrote pdb output for converted", outfile)
        else:
            print("-wip requires PDB or MMCIF input")

    if args.wpi:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates()
            pdb_structure.atom_to_internal_coordinates()
            write_PIC(pdb_structure, outfile + ".pic")
        else:
            print("-wpi requires PIC input")

    if args.wc:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates()
        io = MMCIFIO()
        io.set_structure(pdb_structure)
        io.save(outfile + ".cif")
        print("wrote cif output for", outfile)

    if args.i:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates()

    if args.wi:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates()
        write_PIC(pdb_structure, outfile + ".pic")
        print("wrote pic output for", outfile)

    if args.t or args.tv:
        sp = StringIO()
        if pdb_input:
            with warnings.catch_warnings(record=True) as w:
                # warnings.simplefilter("error")
                # if True:  # try:
                try:
                    if args.cp or args.ti:
                        doCpti(
                            args,
                            "structure_rebuild_test(pdb_structure, args.tv)",
                            pdb_structure,
                        )
                    else:
                        r = structure_rebuild_test(pdb_structure, args.tv)
                    warns = len(w) > 0
                    if args.tv and warns:
                        for wrn in w:
                            print(wrn.message)
                    if not (args.ti or args.cp):
                        print(
                            prot_id, fileNo, r["report"], ("WARNINGS" if warns else "")
                        )
                except Exception as e:
                    print(prot_id, fileNo, "EXCEPTION:", type(e), e)

        elif pic_input:
            pdb_structure.internal_to_atom_coordinates()
            write_PDB(pdb_structure, sp)
            sp.seek(0)
            pdb2 = PDB_parser.get_structure(prot_id, sp)
            pdb2.atom_to_internal_coordinates()
            sp2 = StringIO()
            write_PIC(pdb2, sp2)
            sp2.seek(0)
            inf = open(filename, "r")
            lineCount = 0
            matchCount = 0
            diffCount = 0
            # for line1, line2 in itertools.zip_longest(inf, sp2):
            for line1, line2 in zip(inf, sp2):
                lineCount += 1
                if line1 == line2:
                    matchCount += 1
                else:
                    diffCount += 1
                    if args.tv:
                        print(line1, "!=", line2)
            print(lineCount, matchCount, diffCount)

    if args.trp:
        if not pdb_input:
            print("-trp only for pdb input")
        else:
            doCpti(args, "pdb_structure.atom_to_internal_coordinates()", pdb_structure)

    if args.t6p or args.t6a or args.t6s or args.t6c:
        if pdb_input:
            pdb_copy = pdb_structure.copy()
            pdb_structure.atom_to_internal_coordinates()
        else:
            print("-t6x only for pdb input")
            continue
        alist = []
        if args.t6p:
            alist = ["psi"]
        elif args.t6a:
            alist = ["psi", "phi", "omg", "chi1", "chi2", "chi3", "chi4", "chi5"]
        elif args.t6s:
            alist = ["chi1", "chi2", "chi3", "chi4", "chi5"]
        elif args.t6c:
            alist = ["chi1"]
        if args.ti or args.cp:
            doCpti(args, "test66(pdb_structure, alist)", pdb_structure)
        else:
            test66(pdb_structure, alist)
        r = compare_residues(pdb_structure, pdb_copy, True)
        print(prot_id, fileNo, r["report"])

    if args.rama:
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates()
        for r in pdb_structure.get_residues():
            ric = r.internal_coord
            if ric:
                print(
                    r,
                    ric.get_angle("psi"),
                    ric.get_angle("phi"),
                    ric.get_angle("omg"),
                    ric.get_angle("tau"),
                    ric.get_angle("chi2"),
                    ric.get_length("0C:1N"),
                )
                print(
                    r,
                    ric.get_angle("N:CA:C:1N"),
                    ric.get_angle("-1C:N:CA:C"),
                    ric.get_angle("-1CA:-1C:N:CA"),
                    ric.get_angle("N:CA:C"),
                    ric.get_angle("CA:CB:CG:CD"),
                    None
                    if not ric.rnext
                    else ric.get_length((ric.rak("C"), ric.rnext[0].rak("N"))),
                )
            # print(r.internal_coord.get_dihedral('N:CA:C:O'))

    if args.rama2:
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates()
        nvc1 = {}
        nvpsi = {}
        nvtau = {}
        c1count = 0
        psicount = 0
        taucount = 0

        for r in pdb_structure.get_residues():
            # print(r.internal_coord.get_dihedral('N:CA:C:O'))
            if r.internal_coord:
                ric = r.internal_coord
                chi1 = ric.get_angle("chi1")
                if chi1 is not None:
                    c1count += 1
                    nv = chi1 + 90
                    if nv > 180.0:
                        nv -= 360.0
                    ric.set_angle("chi1", nv)
                    nvc1[str(r)] = nv
                psi = ric.get_angle("psi")
                if psi is not None:
                    psicount += 1
                    nv = psi - 90
                    if nv < -180.0:
                        nv += 360.0
                    ric.set_angle("psi", nv)
                    nvpsi[str(r)] = nv
                tau = ric.get_angle("tau")
                if tau is not None:
                    taucount += 1
                    nv = tau - 5
                    ric.set_angle("tau", nv)
                    nvtau[str(r)] = nv

        pdb2 = IC_duplicate(pdb_structure)
        # write_PIC(pdb2, "foo.pic")
        pdb2.internal_to_atom_coordinates()
        sf = StringIO()
        write_PDB(pdb2, sf)
        sf.seek(0)
        new_pdb = PDB_parser.get_structure("1NEW", sf)
        new_pdb.atom_to_internal_coordinates()
        c1tcount = 0
        psitcount = 0
        tautcount = 0
        for r in new_pdb.get_residues():
            ric = r.internal_coord
            if ric:
                chi1 = ric.get_angle("chi1")
                if chi1 is not None:
                    c1tcount += 1
                    print(
                        f"chi1 {chi1} should be {nvc1[str(r)]} rslt= {chi1 == nvc1[str(r)]}"
                    )
                psi = ric.get_angle("psi")
                if psi is not None:
                    psitcount += 1
                    print(
                        f"psi {psi} should be {nvpsi[str(r)]} rslt= {psi == nvpsi[str(r)]}"
                    )
                tau = ric.get_angle("tau")
                if tau is not None:
                    tautcount += 1
                    print(
                        f"tau {tau} should be {nvtau[str(r)]} rslt= {tau == nvtau[str(r)]}"
                    )

        print(f"chi1count {c1count} {c1tcount} {c1count == c1tcount}")
        print(f"psicount {psicount} {psitcount} {psicount == psitcount}")
        print(f"taucount {taucount} {tautcount} {taucount == tautcount}")

    if args.ws:
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates()
        if args.flex:
            if args.gcb and pic_input:
                pdb_structure.internal_to_atom_coordinates()
                pdb_structure.atom_to_internal_coordinates()  # build C-beta's
            for r in pdb_structure.get_residues():
                if r.internal_coord:
                    r.internal_coord.set_flexible()
        if args.hb:
            for r in pdb_structure.get_residues():
                if r.internal_coord:
                    r.internal_coord.set_hbond()

        write_SCAD(
            pdb_structure,
            outfile + ".scad",
            scale_val,
            pdbid=prot_id,
            backboneOnly=args.backbone,
            maxPeptideBond=args.maxp,
            includeCode=(not args.nsc),
        )
        write_PIC(pdb_structure, "foo.pic")

    fileNo += 1

print("normal termination")
# print(AtomKey.icd)
