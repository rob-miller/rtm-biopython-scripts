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
import numpy as np
import copy

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

# handle old version of rebuild test
fn = structure_rebuild_test
srt_args = fn.__code__.co_varnames[: fn.__code__.co_argcount]
old_srt = False if len(srt_args) > 2 else True
srt_string = (
    "structure_rebuild_test(pdb_structure, args.tv)"
    if old_srt
    else "structure_rebuild_test(pdb_structure, args.tv, args.cc)"
)


# def doCpti(args, cmd, pdb_structure):
def doCpti(args, cmd):
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


def rmse(a0, a1):
    return np.sqrt(((a0 - a1) ** 2).mean())


def test66(pdb_structure, alist, pdb_copy):
    """change dihedral angles 6x60."""

    for chnTpl in zip(pdb_structure.get_chains(), pdb_copy.get_chains()):
        cic0, cic1 = chnTpl[0].internal_coord, chnTpl[1].internal_coord
        aNdxList = []
        for r in chnTpl[0].get_residues():
            ric = r.internal_coord
            if not old_srt:
                if ric is not None:
                    for a in alist:
                        edron = ric.pick_angle(a)
                        if edron is not None:
                            aNdxList.append(edron.ndx)
        r0 = rmse(
            cic0.atomArray[:, 0:3],
            cic1.atomArray[:, 0:3],
            # cic0.dihedraAngle[aNdxList],
            # cic1.dihedraAngle[aNdxList],
        )
        print(f"{r0:.3}", end="")
        for i in range(6):
            for r in chnTpl[0].get_residues():
                ric = r.internal_coord
                if ric:
                    for a in alist:
                        adj60(ric, a)
            cic0.internal_to_atom_coordinates()
            cic0.atom_to_internal_coordinates()
            if not old_srt:
                r1 = rmse(
                    cic0.atomArray[:, 0:3],
                    cic1.atomArray[:, 0:3],
                    # cic0.dihedraAngle[aNdxList],
                    # cic1.dihedraAngle[aNdxList],
                )
                print(f" -> {r1:.3} ", end="")
            pass
            # pdb_structure.internal_to_atom_coordinates()

        if not old_srt:
            print()


def testTau(pdb_structure, pdb_copy):
    """change tau angles 1 degree steps -3, +3, back to 0."""

    for chnTpl in zip(pdb_structure.get_chains(), pdb_copy.get_chains()):
        cic0, cic1 = chnTpl[0].internal_coord, chnTpl[1].internal_coord
        aNdxList = []
        for r in chnTpl[0].get_residues():
            ric = r.internal_coord
            if ric is not None:
                edron = ric.pick_angle("tau")
                if edron is not None:
                    aNdxList.append(edron.ndx)
        if not old_srt:
            r0 = rmse(
                cic0.atomArray[:, 0:3],
                cic1.atomArray[:, 0:3],
                # cic0.hedraAngle[aNdxList],
                # cic1.hedraAngle[aNdxList],
            )
            print(f"{r0:.3}", end="")
        for i in (-3, 1, 1, 1, 1, 1, 1, -1, -1, -1):
            first = True
            for r in chnTpl[0].get_residues():
                ric = r.internal_coord
                if ric and not first:
                    a = ric.get_angle("tau")
                    ric.set_angle("tau", a + i)
                first = False
            cic0.internal_to_atom_coordinates()
            cic0.atom_to_internal_coordinates()
            if not old_srt:
                r1 = rmse(
                    cic0.atomArray[:, 0:3],
                    cic1.atomArray[:, 0:3],
                    # cic0.hedraAngle[aNdxList],
                    # cic1.hedraAngle[aNdxList],
                )
                print(f" -> {r1:.3} ", end="")
            pass
            # pdb_structure.internal_to_atom_coordinates()
    if not old_srt:
        print()


def copyCheck(pdb_src, pdb_cpy):
    for schn, cchn in zip(pdb_src.get_chains(), pdb_cpy.get_chains()):
        schni = schn.internal_coord
        cchni = cchn.internal_coord
        print(
            cchni.hedraL12.base is schni.hedraL12,
            np.may_share_memory(cchni.hedraL12, schni.hedraL12),
        )
        print(
            cchni.atomArray.base is schni.atomArray,
            np.may_share_memory(cchni.atomArray, schni.atomArray),
        )


def color(val, tval, accept, fmt, tol):
    if val is None:
        return f"{'':{fmt}}"
    if tval is None:
        return "\u001b[34m" + f"{val:{fmt}}" + "\u001b[0m"
    if val == tval:
        return f"{val:{fmt}}"
    diff = abs(val - tval)
    # vstr = f"{val:{fmt}}"
    # tvstr = f"{tval:{fmt}}"
    if diff < (tol * 0.1):
        return f"{val:{fmt}}"
    if diff < tol:
        return "\u001b[33m" + f"{val:{fmt}}" + "\u001b[0m"
    if accept:
        return "\u001b[32m" + f"{val:{fmt}}" + "\u001b[0m"
    return "\u001b[31m" + f"{val:{fmt}}" + "\u001b[0m"


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
    "-t", help="test conversion pdb/pic to pic/pdb", action="store_true"
)
arg_parser.add_argument(
    "-cc",
    help="rebuild test copies coordspace matrices (faster, less rigorous)",
    action="store_true",
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
    "-tt",
    help="read pdb, change every tau -3 to +3 and back to 0 deg & assemble, back to pdb and test",
    action="store_true",
)
arg_parser.add_argument(
    "-t33",
    help="read pdb, residue 33 add 33 to each dihedral and check, add 3.3 to tau and check",
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
arg_parser.add_argument(
    "-ser", help="revert to serial algorithm for residue assembly", action="store_true"
)
arg_parser.add_argument("-nh", help="ignore hydrogens on PDB read", action="store_true")
arg_parser.add_argument(
    "-d2h", help="swap D (deuterium) for H on PDB read", action="store_true"
)
arg_parser.add_argument(
    "-backbone", help="only backbone heavy atoms on PDB read", action="store_true"
)
arg_parser.add_argument(
    "-amidep",
    help="add only amide proton for PDB read (use with -nh, -backbone)",
    action="store_true",
)
arg_parser.add_argument(
    "-cb",
    help="accept C-beta for PDB read (use with -backbone, -gcb)",
    action="store_true",
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

arg_parser.add_argument(
    "-tip",
    help="test just convert to internal coordinates and back",
    action="store_true",
)

arg_parser.add_argument("-p", help="print header info", action="store_true")

args = arg_parser.parse_args()

# print(args)

if args.nh:
    IC_Residue.accept_atoms = IC_Residue.accept_mainchain
if args.d2h:
    IC_Residue.accept_atoms += IC_Residue.accept_deuteriums
    AtomKey.d2h = True
if args.backbone:
    IC_Residue.accept_atoms = IC_Residue.accept_backbone
if args.amidep:
    IC_Residue.accept_atoms += ("H",)
if args.cb:
    IC_Residue.accept_atoms += ("CB",)
if args.gcb:
    IC_Residue.gly_Cbeta = True
if args.maxp:
    IC_Chain.MaxPeptideBond = float(args.maxp)
if args.ser:
    IC_Chain.ParallelAssembleResidues = False
if args.scale:
    scale_val = float(args.scale)
if args.skip_count:
    args.skip_count = int(args.skip_count)
if args.limit_count:
    args.limit_count = int(args.limit_count)
if args.ti_iter:
    ti_iter = int(args.ti_iter)

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

    if args.tip:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates()
            for chn in pdb_structure.get_chains():
                cic = chn.internal_coord
                cic.atomArrayValid[:] = False
                cic.atomArrayValid[0:3] = True
            pdb_structure.internal_to_atom_coordinates()
        elif pic_input:
            pdb_structure.internal_to_atom_coordinates()
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
                if True:  # try:
                    if args.cp or args.ti:
                        doCpti(args, srt_string)
                    else:
                        if old_srt:
                            r = structure_rebuild_test(pdb_structure, args.tv)
                        else:
                            r = structure_rebuild_test(pdb_structure, args.tv, args.cc)
                    warns = len(w) > 0
                    if args.tv and warns:
                        for wrn in w:
                            print(wrn.message)
                    if not (args.ti or args.cp):
                        print(
                            prot_id, fileNo, r["report"], ("WARNINGS" if warns else "")
                        )
                # except Exception as e:
                #    print(prot_id, fileNo, "EXCEPTION:", type(e), e)

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

    if args.t6p or args.t6a or args.t6s or args.t6c or args.tt:
        if pic_input:
            # pdb_copy = read_PIC(
            #    gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename
            # )
            # pdb_copy = copy.deepcopy(pdb_structure)
            pdb_structure.internal_to_atom_coordinates()
            # pdb_copy.internal_to_atom_coordinates()
            pdb_copy = copy.deepcopy(pdb_structure)
        elif pdb_input:
            pdb_copy = pdb_structure.copy()
            pdb_structure.atom_to_internal_coordinates()
            pdb_copy.atom_to_internal_coordinates()
        # else:
        #    print("-t6x only for pdb input")
        #    continue
        alist = []
        if args.t6p:
            alist = ["psi"]
        elif args.t6a:
            alist = ["psi", "omg", "phi", "chi1", "chi2", "chi3", "chi4", "chi5"]
        elif args.t6s:
            alist = ["chi1", "chi2", "chi3", "chi4", "chi5"]
        elif args.t6c:
            alist = ["chi1"]
        if args.tt:
            if args.ti or args.cp:
                doCpti(args, "testTau(pdb_structure, pdb_copy)")
            else:
                testTau(pdb_structure, pdb_copy)
        elif args.ti or args.cp:
            doCpti(args, "test66(pdb_structure, alist, pdb_copy)")
        else:
            test66(pdb_structure, alist, pdb_copy)
        r = compare_residues(pdb_structure, pdb_copy, True)
        print(prot_id, fileNo, r["report"])

    if args.t33:
        resTarg = 32  # 32  # account for 0-start child_list
        delta = 33  # degrees to change
        tdelta = delta  # / 10.0  # more realistic for bond angle
        if pic_input:
            pdb_structure.internal_to_atom_coordinates()
            # pdb_structure.atom_to_internal_coordinates()
            pdb_copy = copy.deepcopy(pdb_structure)
        elif pdb_input:
            pdb_copy = pdb_structure.copy()
            pdb_structure.atom_to_internal_coordinates()
            pdb_copy.atom_to_internal_coordinates()
        alist = ["omg", "phi", "psi", "chi1", "chi2", "chi3", "chi4", "chi5", "tau"]
        print()
        for chnTpl in zip(pdb_structure.get_chains(), pdb_copy.get_chains()):
            cic0, cic1 = chnTpl[0].internal_coord, chnTpl[1].internal_coord
            c = 3
            for a in alist:
                if c > 0:
                    c -= 1
                else:
                    pass
                    # break

                ric = chnTpl[0].child_list[resTarg].internal_coord
                ric1 = chnTpl[1].child_list[resTarg].internal_coord
                print(
                    ric.rbase[0],
                    ric.rbase[2],
                    ric.rbase[1] if ric.rbase[1] is not None else "",
                    a,
                )
                print("         ", end="")
                for ar in alist:
                    print(f"{ar: >8}", end=" ")
                print()
                print("start:   ", end="")
                for ar in alist:
                    ang = ric.get_angle(ar)
                    print(f"{ang if ang is not None else '':8.4}", end=" ")
                print()
                print("target:  ", end="")

                setter = True  # use edra.setter or direct array manipulation
                if setter:
                    ang0 = ric.get_angle(a)
                    if ang0 is not None:
                        targ = ang0 + (tdelta if a == "tau" else delta)
                        if targ > 180.0:
                            targ -= 360.0
                        ric.set_angle(a, targ)
                else:
                    try:
                        andx = ric.pick_angle(a).ndx
                        if a == "tau":
                            cic0.hedraAngle[andx] += tdelta
                            cic0.hAtoms_needs_update[andx] = True
                            cic0.atomArrayValid[cic0.h2aa[andx]] = False
                        else:
                            cic0.dihedraAngle[andx] += delta
                            if cic0.dihedraAngle[andx] > 180.0:
                                cic0.dihedraAngle[andx] -= 360.0
                            cic0.dihedraAngleRads[andx] = np.deg2rad(
                                cic0.dihedraAngle[andx]
                            )
                            cic0.dAtoms_needs_update[andx] = True
                            cic0.atomArrayValid[cic0.d2aa[andx]] = False
                    except AttributeError:
                        pass  # if residue does not have e.g. chi5

                # if a == "omg":
                #    write_PIC(pdb_structure, "omg1.pic")

                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")

                print()
                print("result:  ", end="")
                cic0.internal_to_atom_coordinates()
                cic0.atom_to_internal_coordinates()

                # if a == "omg":
                #    write_PIC(pdb_structure, "omg1.pic")

                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")
                print()
                print("reset:   ", end="")
                if setter:
                    ric.set_angle(a, ang0)
                else:
                    if a == "tau":
                        cic0.hedraAngle[andx] = cic1.hedraAngle[andx]
                        cic0.hAtoms_needs_update[andx] = True
                        cic0.atomArrayValid[cic0.h2aa[andx]] = False
                    else:
                        cic0.dihedraAngle[andx] = cic1.dihedraAngle[andx]
                        cic0.dihedraAngleRads[andx] = np.deg2rad(
                            cic0.dihedraAngle[andx]
                        )
                        cic0.dAtoms_needs_update[andx] = True
                        cic0.atomArrayValid[cic0.d2aa[andx]] = False
                        # cic0.dcsValid[andx] = False

                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")
                print()
                cic0.internal_to_atom_coordinates()
                cic0.atom_to_internal_coordinates()

                print("finish:  ", end="")
                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")
                print()

                print()

        r = compare_residues(pdb_structure, pdb_copy, True)
        print(prot_id, fileNo, r["report"])
        pass

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
