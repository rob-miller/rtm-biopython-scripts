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
import ast

import matplotlib.pyplot as plt

import gzip
import warnings

# import time

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

from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
from Bio.PDB.internal_coords import IC_Residue, AtomKey, IC_Chain, Dihedron
from Bio.PDB.SCADIO import write_SCAD
from Bio.PDB.Chain import Chain
from Bio import SeqIO

PDB_repository_base = None

pdbDirs = ["/Users/rob/data/pdb/", "/media/data/rcsb/pdb", "/Volumes/data/pdb/"]
for d in pdbDirs:
    if os.path.isdir(d):
        PDB_repository_base = d
        break

resnameNdx = AtomKey.fields.resname
atmNdx = AtomKey.fields.atm

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
    else "structure_rebuild_test(pdb_structure, args.tv, args.q)"
)


# def doCpti(args, cmd, pdb_structure):
def doCpti(args, cmd):
    """Profile handler for cmd."""
    if args.ti:
        tr = timeit.repeat(
            stmt=cmd,
            repeat=int(ti_iter),
            number=1,
            globals=globals(),
        )
        print(min(tr))

    elif args.cp:
        cmd2 = "for i in range(0, " + str(ti_iter) + "): " + cmd
        cProfile.runctx(
            cmd2,
            globals(),
            locals(),
            "Profile.prof",
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
        if not old_srt:
            if not hasattr(cic0, "atomArray"):
                continue

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
            cic0.internal_to_atom_coordinates(verbose=args.v)
            cic0.atom_to_internal_coordinates(verbose=args.v)
            if not old_srt:
                r1 = rmse(
                    cic0.atomArray[:, 0:3],
                    cic1.atomArray[:, 0:3],
                    # cic0.dihedraAngle[aNdxList],
                    # cic1.dihedraAngle[aNdxList],
                )
                print(f" -> {r1:.3} ", end="")
            pass
            # pdb_structure.internal_to_atom_coordinates(verbose=args.v)

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
            cic0.internal_to_atom_coordinates(verbose=args.v)
            cic0.atom_to_internal_coordinates(verbose=args.v)
            if not old_srt:
                r1 = rmse(
                    cic0.atomArray[:, 0:3],
                    cic1.atomArray[:, 0:3],
                    # cic0.hedraAngle[aNdxList],
                    # cic1.hedraAngle[aNdxList],
                )
                print(f" -> {r1:.3} ", end="")
            pass
            # pdb_structure.internal_to_atom_coordinates(verbose=args.v)
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
        return f"{'  -  ':{fmt}}"
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


def checkAllDihedra(a, ric, cic0, cic1, setter):
    andx = 0
    if not setter:
        try:
            if a is not None:
                andx = ric.pick_angle(a).ndx
            diff = np.abs(cic0.dihedraAngle - cic1.dihedraAngle)
            diffmask = diff < 0.001
            if a is not None:
                diff[andx] = 0
            diff[diffmask] = 0
            rslt = np.sum(diff)
            # print(f" : {rslt:{8.4}}", end=" ")
            print(color(rslt, 0.0, False, 8.4, 0.001), end=" ")
            if rslt > 0.001:
                print(" ", end=" ")
        except AttributeError:
            pass


def pic_load(filename):
    return read_PIC(
        gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename,
        verbose=args.v,
        defaults=args.default,
    )


def cif_load(filename, prot_id):
    return CIF_parser.get_structure(
        prot_id,
        gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename,
    )


def pdb_load(filename, prot_id):
    return PDB_parser.get_structure(
        prot_id,
        gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename,
    )


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
    help="text file each line beginning PDB id, optional chain id (7RSAA) to read from {0} as .ent.gz".format(
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


arg_parser.add_argument("-hstat", help="collect hedra stats", action="store_true")
arg_parser.add_argument("-dstat", help="collect dihedra stats", action="store_true")

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
    "default " + str(IC_Chain.MaxPeptideBond) + " angstroms",
)
arg_parser.add_argument(
    "-t", help="test conversion pdb/pic to pic/pdb", action="store_true"
)
arg_parser.add_argument(
    "-tok", help="test conversion pdb<>pic output only OK", action="store_true"
)
arg_parser.add_argument(
    "-q", help="rebuild test only checks atomArray match", action="store_true"
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
arg_parser.add_argument("-na", help="no altloc atoms on PDB read", action="store_true")
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
    "-stdAA", help="only std 20 AAs on PDB read - no UNK etc", action="store_true"
)
arg_parser.add_argument(
    "-cb",
    help="accept C-beta for PDB read (use with -backbone, -gcb)",
    action="store_true",
)
arg_parser.add_argument("-gcb", help="generate GLY C-beta atoms", action="store_true")
arg_parser.add_argument(
    "-rama", help="print psi, phi, omega values (sample code)", action="store_true"
)
arg_parser.add_argument(
    "-rama2", help="test setting dihedral angles (sample code)", action="store_true"
)

arg_parser.add_argument(
    "-i", help="test just convert to internal coordinates", action="store_true"
)

arg_parser.add_argument(
    "-tip",
    help="test just convert to internal coordinates and back",
    action="store_true",
)

arg_parser.add_argument(
    "-tdup",
    help="populate atom/internal coords as needed, deep copy and compare",
    action="store_true",
)

arg_parser.add_argument("-p", help="print header info", action="store_true")
arg_parser.add_argument("-pri", help="list primary angles", action="store_true")
arg_parser.add_argument(
    "-default", help="read .pic with defaults as needed", action="store_true"
)

arg_parser.add_argument(
    "-nwh", help="no hedra written to .pic files", action="store_true"
)
arg_parser.add_argument(
    "-n2nd", help="no secondary dihedra written to .pic files", action="store_true"
)
arg_parser.add_argument(
    "-pCut",
    dest="pCut",
    help="minimum std dev in ref db to write primary angle to pic files",
)
arg_parser.add_argument(
    "-hCut",
    dest="hCut",
    help="minimum angle std dev in ref db to write hedron to pic files",
)
arg_parser.add_argument(
    "-nbf", help="no b-factor records written to .pic files", action="store_true"
)
arg_parser.add_argument(
    "-pfc",
    dest="pfc",
    help=f"specify picFlags class ({IC_Residue.pic_flags})",
)

arg_parser.add_argument("-v", help="verbose", action="store_true")

arg_parser.add_argument(
    "-start",
    dest="start_offset",
    help="start position param chain.i2ac, -wip chain only",
)
arg_parser.add_argument(
    "-fin", dest="fin_offset", help="fin position param chain.i2ac, -wip chain only"
)

arg_parser.add_argument(
    "-e",
    dest="execfile",
    help="python file to eecute first, e.g. override hedra_defaults",
)
arg_parser.add_argument("-ramaOut", help="ouput dihedrals", action="store_true")
arg_parser.add_argument(
    "-ramaIn",
    dest="ramaIn_file",
    help="load dihedrals",
)
arg_parser.add_argument("-seq", help="print sequence from atoms", action="store_true")
arg_parser.add_argument(
    "-s2s",
    help="read protein sequence fasta file as PIC file with defaults",
    action="store_true",
)
arg_parser.add_argument(
    "-dinfo", help="print info about each dihedron", action="store_true"
)
arg_parser.add_argument(
    "-dpd", help="build copy of input from distance plot", action="store_true"
)
arg_parser.add_argument(
    "-dpc",
    dest="dpc",
    help="display distplot difference with comparison structure DPC",
)

args = arg_parser.parse_args()

# print(args)

if args.execfile:
    exec(open(args.execfile).read())
if args.nh:
    IC_Residue.accept_atoms = IC_Residue.accept_mainchain
if args.stdAA:
    IC_Residue.accept_resnames = ()  # no UNK etc
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
if args.na:
    IC_Residue.no_altloc = True

pCut = float(args.pCut) if args.pCut is not None else None
hCut = float(args.hCut) if args.hCut is not None else None

picFlags = IC_Residue.picFlagsDefault
if args.nwh:
    picFlags &= ~IC_Residue.pic_flags.hedra
if args.n2nd:
    picFlags &= ~IC_Residue.pic_flags.secondary
if args.nbf:
    picFlags &= ~IC_Residue.pic_flags.bFactors

if args.pfc:
    picFlags = int(args.pfc)


start_offset = int(args.start_offset) if args.start_offset else None
fin_offset = int(args.fin_offset) if args.fin_offset else None

toProcess = args.file
pdbidre = re.compile(r"(^\d(\w\w)\w)(\w)?$")

hedra_stats = {}
dihedra_stats = {}
primary_angles = {}  # set()

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

if args.v or (len(toProcess) == 0):
    if PDB_repository_base is not None:
        print("pdb base = ", PDB_repository_base)
    print(sys.version)

if len(toProcess):
    if args.v:
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
    outfile = os.path.basename(target)

    prot_id = ""
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
            pdb_structure = pic_load(filename)
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
            pdb_structure = pdb_load(filename, prot_id)
        except FileNotFoundError:  # OSError python2
            # print(print(os.getcwd()), filename, "not found")
            print(os.getcwd(), filename, "not found")
            fileNo += 1
            continue
        except Exception:
            pass  # continue  # try as cif below

        if pdb_structure is None or pdb_structure.child_dict == {}:
            try:
                pdb_structure = cif_load(filename, prot_id)
                cif_input = True
            except Exception:
                # print("unable to open ", filename, " as PDB, MMCIF or PIC file format.")
                pdb_input = False

        if pdb_structure is not None:
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

    if args.seq:
        for record in SeqIO.parse(filename, "pdb-atom"):
            print(f">{record.id} {record.annotations['name']}")
            out = [(record.seq[i : i + 80]) for i in range(0, len(record.seq), 80)]
            for lin in out:
                print(lin)

    if args.s2s:
        seqIter = SeqIO.parse(filename, "fasta")
        for record in seqIter:
            pdb_structure = read_PIC_seq(record)
        if pdb_structure is not None:
            pic_input = True

    # ok have a structure file of some sort to work with
    if args.v:
        if pdb_input:
            # if not (args.t or args.tok or args.hstat or args.dstat or args.pri):
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

    if pdb_input:
        if args.hstat:
            IC_Residue.accept_resnames = ()  # no UNK etc
            IC_Residue.accept_atoms = IC_Residue.accept_mainchain + ("H",)
            IC_Residue.no_altloc = True
            # r = structure_rebuild_test(pdb_structure, False, False)
            # if r["pass"]:
            #    print(target, fileNo, r["report"])
            pdb_structure.atom_to_internal_coordinates()
            if not isinstance(pdb_structure, Chain):  # get 1st chain if needed
                pdb_structure = pdb_structure[0].child_list[0]
            for h in pdb_structure.internal_coord.hedra.values():  # chain
                hc = h.xrh_class if hasattr(h, "xrh_class") else h.rdh_class
                # hc = h.crdh_class
                try:
                    hedra_stats[hc][0].append(h.len12)
                    hedra_stats[hc][1].append(h.angle)
                    hedra_stats[hc][2].append(h.len23)
                    hedra_stats[hc][3].add(h.dh_class)

                except KeyError:
                    hedra_stats[hc] = [
                        [h.len12],
                        [h.angle],
                        [h.len23],
                        {h.dh_class},
                    ]
        if args.pri:
            pdb_structure.atom_to_internal_coordinates()
            if not isinstance(pdb_structure, Chain):  # get 1st chain if needed
                pdb_structure = pdb_structure[0].child_list[0]
            cic = pdb_structure.internal_coord
            for d in cic.dihedra.values():  # chain
                if d.primary:
                    # primary_angles.add(d.dh_class)
                    """rtm
                    if hasattr(d, "psi_class") and d.psi_class not in primary_angles:
                        # primary_angles.add(d.psi_class)
                        primary_angles[d.psi_class] = tuple(
                            d.aks[x].akl[atmNdx] for x in range(4)
                        )
                    elif hasattr(d, "omg_class") and d.omg_class not in primary_angles:
                        # primary_angles.add(d.omg_class)
                        primary_angles[d.omg_class] = tuple(
                            d.aks[x].akl[atmNdx] for x in range(4)
                        )
                    elif hasattr(d, "phi_class") and d.phi_class not in primary_angles:
                        # primary_angles.add(d.phi_class)
                        primary_angles[d.phi_class] = tuple(
                            d.aks[x].akl[atmNdx] for x in range(4)
                        )
                    elif (
                        hasattr(d, "altCB_class") and d.phi_class not in primary_angles
                    ):
                        # primary_angles.add(d.phi_class)
                        primary_angles[d.altCB_class] = tuple(
                            d.aks[x].akl[atmNdx] for x in range(4)
                        )
                    elif d.rdh_class not in primary_angles:
                        primary_angles[d.rdh_class] = tuple(
                            d.aks[x].akl[atmNdx] for x in range(4)
                        )
                        # primary_angles.add(d.rdh_class)
                    """

        if args.dstat:
            IC_Residue.accept_resnames = ()  # no UNK etc
            IC_Residue.accept_atoms = IC_Residue.accept_mainchain + ("H",)
            IC_Residue.no_altloc = True
            if False:
                r = structure_rebuild_test(pdb_structure, verbose=False, quick=True)
                if not r["pass"]:
                    continue
            else:
                pdb_structure.atom_to_internal_coordinates()
            if not isinstance(pdb_structure, Chain):  # get 1st chain if needed
                pdb_structure = pdb_structure[0].child_list[0]
            cic = pdb_structure.internal_coord
            for d in cic.dihedra.values():  # chain
                if d.primary:
                    dclass = d.pclass
                    try:
                        dihedra_stats[dclass][0].append(d.angle)
                        dihedra_stats[dclass][1][int(np.round(d.angle + 180))] += 1
                    except KeyError:
                        dihedra_stats[dclass] = [[d.angle], np.zeros(361, int), {}]
                        dihedra_stats[dclass][1][int(np.round(d.angle + 180))] = 1
                        dihedra_stats[dclass].append(
                            tuple(d.aks[x].akl[atmNdx] for x in range(4))
                        )
                    """
                    if dclass in [
                        "YCGYCD1YCE1YCZ",
                        "WCE2WCZ2WCH2WCZ3",
                        "WCD2WCE2WCZ2WCH2",
                        "RCDRNERCZRNH1",
                        "FCGFCD1FCE1FCZ",
                    ]:
                        print(target, d.ric.residue, dclass, d.angle)
                    """
                    for d2k in cic.id3_dh_index[d.id3]:
                        d2base = cic.dihedra[d2k]
                        if d2base.dh_class == d.dh_class:
                            continue

                        d2list = [d2base]
                        try:
                            for d2r in cic.id3_dh_index[d2base.id32[::-1]]:
                                d2list.append(cic.dihedra[d2r])
                        except KeyError:
                            pass

                        for d2 in d2list:
                            try:
                                d2class = d2.altCB_class
                            except AttributeError:
                                d2class = d2.rdh_class

                            # know angles do not overlap so diff !~ 0
                            delta = d.difference(d2)
                            try:
                                dihedra_stats[dclass][2][d2class][0].append(delta)
                            except KeyError:
                                dihedra_stats[dclass][2][d2class] = [[delta]]

                # a = int(np.round(d.angle))
                # if a < 0:
                #    a += 360
                # dihedra_stats[d.crdh_class][1][a] += 1

    if args.ramaIn_file:
        dihedlist = open(args.ramaIn_file, "r")
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
        for zlist in zip(pdb_structure.get_residues(), dihedlist):
            ric = zlist[0].internal_coord
            if ric:
                rdlist = [
                    ric.dihedra[dk]
                    for dk in sorted(ric.dihedra)
                    if ric.dihedra[dk].primary
                ]
                # rdlist = [d for d in ric.dihedra.values() if d.primary]
                dlist = ast.literal_eval(zlist[1])
                for zlist2 in zip(rdlist, dlist):
                    ric.bond_set(zlist2[0].aks, zlist2[1])

    if args.ramaOut:
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
        for r in pdb_structure.get_residues():
            ric = r.internal_coord
            if ric:
                # rdlist = [d.angle for d in ric.dihedra.values() if d.primary]
                dlist = [
                    ric.dihedra[dk].angle
                    for dk in sorted(ric.dihedra)
                    if ric.dihedra[dk].primary
                ]
                # for d in ric.dihedra.values():
                #    if d.primary:
                #        dlist.append(d.angle)
                print(dlist)

    if args.wp:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
        write_PDB(pdb_structure, outfile + ".pdb")
        print("wrote pdb output for", outfile)

    if args.wip:
        if pdb_input:
            pdb2 = IC_duplicate(pdb_structure)
            if isinstance(pdb_structure, Chain):
                pdb2[0].child_list[0].internal_to_atom_coordinates(
                    verbose=args.v, start=start_offset, fin=fin_offset
                )
            else:
                pdb2.internal_to_atom_coordinates(verbose=args.v)
            write_PDB(pdb2, outfile + ".pdb")
            print("wrote pdb output for converted", outfile)
        else:
            print("-wip requires PDB or MMCIF input")

    if args.wpi:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            write_PIC(pdb_structure, outfile + ".pic", picFlags, hCut, pCut)
        else:
            print("-wpi requires PIC input")

    if args.wc:
        if pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
        io = MMCIFIO()
        io.set_structure(pdb_structure)
        io.save(outfile + ".cif")
        print("wrote cif output for", outfile)

    if args.i:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)

    if args.dinfo:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
        pdb_structure.internal_to_atom_coordinates(verbose=args.v)
        for r in pdb_structure.get_residues():
            ric = r.internal_coord
            cic = ric.cic
            if ric:
                for k, d in ric.dihedra.items():
                    # print(k, d.angle)
                    h1 = d.hedron1
                    h2 = d.hedron2
                    oa = h1.len12 if d.reverse else h1.len23
                    ob = h1.len23 if d.reverse else h1.len12
                    ac = h2.len12 if d.reverse else h2.len23
                    oc = np.sqrt(
                        oa * oa + ac * ac - 2 * oa * ac * np.cos(np.deg2rad(h2.angle))
                    )
                    ab = np.sqrt(
                        oa * oa + ob * ob - 2 * oa * ob * np.cos(np.deg2rad(h1.angle))
                    )
                    bndx = (
                        cic.atomArrayIndex[h1.aks[2]]
                        if d.reverse
                        else cic.atomArrayIndex[h1.aks[0]]
                    )
                    ondx = cic.atomArrayIndex[h1.aks[1]]
                    cndx = (
                        cic.atomArrayIndex[h2.aks[0]]
                        if d.reverse
                        else cic.atomArrayIndex[h2.aks[2]]
                    )
                    b = cic.atomArray[bndx]
                    o = cic.atomArray[ondx]
                    c = cic.atomArray[cndx]
                    bc = np.linalg.norm(c - b)

                    Ws = (ab + ac + bc) / 2
                    Xs = (ob + bc + oc) / 2
                    Ys = (oa + ac + oc) / 2
                    Zs = (oa + ob + ab) / 2
                    Wsqr = Ws * (Ws - ab) * (Ws - ac) * (Ws - bc)
                    Xsqr = Xs * (Xs - ob) * (Xs - bc) * (Xs - oc)
                    Ysqr = Ys * (Ys - oa) * (Ys - ac) * (Ys - oc)
                    Zsqr = Zs * (Zs - oa) * (Zs - ob) * (Zs - ab)
                    Hsqr = (
                        4 * oa * oa * bc * bc
                        - np.square((ob * ob + ac * ac) - (oc * oc + ab * ab))
                    ) / 16
                    Jsqr = (
                        4 * ob * ob * ac * ac
                        - np.square((oc * oc + ab * ab) - (oa * oa + bc * bc))
                    ) / 16
                    Ksqr = (
                        4 * oc * oc * ab * ab
                        - np.square((oa * oa + bc * bc) - (ob * ob + ac * ac))
                    ) / 16

                    Y = np.sqrt(Ysqr)
                    Z = np.sqrt(Zsqr)
                    X = np.sqrt(Xsqr)
                    W = np.sqrt(Wsqr)

                    cosOA = (Ysqr + Zsqr - Hsqr) / (2 * Y * Z)
                    cosOB = (Zsqr + Xsqr - Jsqr) / (2 * Z * X)
                    cosOC = (Xsqr + Ysqr - Ksqr) / (2 * X * Y)
                    cosBC = (Wsqr + Xsqr - Hsqr) / (2 * W * X)
                    cosCA = (Wsqr + Ysqr - Jsqr) / (2 * W * Y)
                    cosAB = (Wsqr + Zsqr - Ksqr) / (2 * W * Z)

                    OA = np.rad2deg(np.arccos(cosOA))
                    OB = np.rad2deg(np.arccos(cosOB))
                    OC = np.rad2deg(np.arccos(cosOC))
                    BC = np.rad2deg(np.arccos(cosBC))
                    CA = np.rad2deg(np.arccos(cosCA))
                    AB = np.rad2deg(np.arccos(cosAB))

                    diff = np.absolute(OA - d.angle)
                    err = True if diff > 1.0 else False
                    """
                    cosOA2 = (
                        Wsqr
                        - Xsqr
                        - Ysqr
                        - Zsqr
                        + (2 * Z * X * cosOB)
                        + (2 * X * Y * cosOC)
                    ) / (-2 * Y * Z)

                    OA2 = np.rad2deg(np.arccos(cosOA2))
                    """
                    print(
                        d.dh_class,
                        ("B" if err else "X"),
                        ("R" if d.reverse else "F"),
                        d.angle,
                        OA,
                        # OA2,
                        OB,
                        OC,
                        BC,
                        CA,
                        AB,
                    )

                    # print(oa, ob, ac, ab, oc, bc)
                    # print(Ys, Zs)
                    # print(Ysqr, Zsqr, Hsqr)
                    # print(Y, Z)
                    # print(cosOA, OA)
                    pass

    if args.dpd:
        # IC_Residue.no_altloc = True
        # IC_Residue.accept_atoms = IC_Residue.accept_mainchain
        if pic_input:
            pdb2 = pic_load(filename)
            pdb_structure.internal_to_atom_coordinates()
        elif pdb_input:
            if cif_input:
                pdb2 = cif_load(filename, prot_id)
            else:
                pdb2 = pdb_load(filename, prot_id)
            pdb_structure.atom_to_internal_coordinates()
        # just get a chain
        if isinstance(pdb_structure, Chain):
            chn = pdb_structure
        else:
            for chn in pdb_structure.get_chains():
                break
        cic = chn.internal_coord
        distplot = cic.distance_plot()
        dihedra_signs = cic.dihedral_signs()

        if isinstance(pdb2, Chain):
            chn2 = pdb2
        else:
            if pdbidMatch is not None and pdbidMatch.group(3) is not None:
                chn2 = pdb2[0][pdbidMatch.group(3)]
            elif args.sel_chain is not None:
                chn2 = pdb2[0][args.sel_chain]
            else:
                for chn2 in pdb2.get_chains():
                    break
        if pic_input:
            # read internal_coords, do not create new IC_Chain
            cic2 = chn2.internal_coord
        else:
            # read atoms, need internal_coords but not initialised
            cic2 = chn2.internal_coord = IC_Chain(chn2)

        cic2.init_edra()
        cic2.distplot_to_dh_arrays(distplot)
        cic2.distance_to_internal_coordinates(dihedra_signs)

        # works:
        # pdb3 = IC_duplicate(chn2)
        # for chn3 in pdb3.get_chains():
        #    break
        # cic2 = chn3.internal_coord

        # also works:
        # write_PIC(cic2.chain, "foo.pic")

        # cic2.build_atomArray()  # copies from Atom coords
        # cic2.clear_ic()
        cic2.atomArrayValid[:] = False
        # cic2.dcsValid[:] = False
        # cic2.hAtoms_needs_update[:] = True
        # cic2.dAtoms_needs_update[:] = True
        cic2.atomArray = np.zeros((cic2.AAsiz, 4), dtype=np.float64)  # redundant
        cic2.atomArray[:, 3] = 1.0

        cic2.copy_initNCaCs(cic)
        cic2.internal_to_atom_coordinates()

        dp2 = cic2.distance_plot()
        dpdiff = np.abs(distplot - dp2)
        maxd = np.amax(dpdiff)
        chnbrk = len(cic.initNCaCs) - 1
        # print(np.size(distplot, 0), np.size(dp2, 0), np.size(dpdiff, 0))
        print(
            target,
            fileNo,
            "size= ",
            np.size(dpdiff, 0),
            "max difference= ",
            maxd,
            f"breaks= {chnbrk}" if maxd > 0.0001 else "",
        )

        # r = compare_residues(chn, chn2, verbose=True)
        # print(r)
        # plt.imshow(distplot, cmap="hot", interpolation="nearest")
        # plt.show()
        # plt.imshow(dp2, cmap="hot", interpolation="nearest")
        # plt.show()

        # plt.imshow(dpdiff, cmap="hot", interpolation="nearest")
        # plt.show()

        """
        klist = list(cic2.atomArrayIndex.keys())
        for i in range(32):
            print(f"{klist[i]}[{dpdiff[3,i]:.2f}] ", end="")
        print()
        """
        # """

        dpdiff[dpdiff < 0.001] = 0.0

        if False:
            aai = cic2.atomArrayIndex
            for ric in cic2.ordered_aa_ic_list:
                # rmask = np.asarray([aai[ak] for ak in sorted(ric.ak_set)])
                # print(sorted(ric.ak_set))
                # https://www.semicolonworld.com/question/55316/subsetting-a-2d-numpy-array
                # rdp_orig = distplot[np.ix_(rmask, rmask)]
                # rdp2 = dp2[np.ix_(rmask, rmask)]
                # rdpdiff = np.abs(rdp_orig - rdp2)
                # print(rdpdiff[3, :])
                for d in ric.dihedra:
                    d1 = cic.dihedra[d]
                    d2 = cic2.dihedra[d]
                    diff = np.absolute(d1.angle - d2.angle)
                    err = True if diff > 0.000001 else False
                    # h2 = d2.h1key if d2.reverse else d2.h2key
                    try:
                        d2lst = cic2.id3_dh_index[d2.id32]
                    except KeyError:
                        d2lst = []
                    try:
                        d2rlst = cic2.id3_dh_index[d2.id32[::-1]]
                    except KeyError:
                        d2rlst = []
                    a4lst = [dx[3] for dx in d2lst]
                    a4rlst = [dx[3] for dx in d2rlst]
                    print(
                        f"{ric.rbase[0]:4d}{ric.rbase[2]} {d1.dh_class:12s}",
                        f"{('R' if d1.reverse else 'F')} {('B' if err else 'X')}",
                        f"{d1.angle:12.6f} {d2.angle:12.6f}",
                        # f"{d2.id:32s} {d2.aks[0]} {list(a4lst)}",
                        f"{d2.id32} {d2.aks[0]} {list(a4lst)}",
                        # list(f"{aai[d2.aks[0]]},{aai[a4]}" for a4 in a4lst),
                        list(f"{dpdiff[aai[a4], aai[d2.aks[0]]]:.4f}" for a4 in a4lst),
                        "   ",
                        f"{d2.aks[3]} {list(a4rlst)}",
                        list(f"{dpdiff[aai[a4], aai[d2.aks[0]]]:.4f}" for a4 in a4rlst),
                        # list(f"{dp2[aai[a4], aai[d2.aks[0]]]:.2f}" for a4 in a4lst),
                        # list(f"{distplot[aai[a4], aai[d2.aks[0]]]:.2f}" for a4 in a4lst),
                    )

                # plt.imshow(rdpdiff, cmap="hot", interpolation="nearest")
                # plt.show()
                # exit()

            # """
        if False:
            plt.imshow(dpdiff, cmap="hot", interpolation="nearest")
            plt.show()

        # write_PDB(pdb2, "foo.pdb")
        # write_PIC(pdb2, "foo.pic")
        pass

    if args.dpc:
        # print(args.dpc)
        s2 = pic_load(args.dpc)
        if s2 is not None:
            s2.internal_to_atom_coordinates()
        else:
            s2 = pdb_load(args.dpc)
            if s2 is None:
                s2 = cif_load(args.dpc)
            if s2 is not None:
                s2.atom_to_internal_coordinates()
        if s2 is not None:
            if isinstance(pdb_structure, Chain):
                chn = pdb_structure
            else:
                for chn in pdb_structure.get_chains():
                    break
            if not hasattr(chn, "internal_coord"):
                chn.atom_to_internal_coordinates()
            elif pic_input:
                chn.internal_to_atom_coordinates()
            cic = chn.internal_coord
            distplot = cic.distance_plot()
            dihedra_signs = cic.dihedral_signs()
            if isinstance(s2, Chain):
                chn2 = s2
            else:
                for chn2 in s2.get_chains():
                    break
            cic2 = chn2.internal_coord
            dp2 = cic2.distance_plot()
            dpdiff = np.abs(distplot - dp2)
            dpdiff[dpdiff < 0.0001] = 0.0
            # plt.imshow(distplot, cmap="hot", interpolation="nearest")
            # plt.show()
            # plt.imshow(dp2, cmap="hot", interpolation="nearest")
            # plt.show()

            plt.imshow(dpdiff, cmap="hot", interpolation="nearest")
            plt.show()

        pass

    if args.tdup:
        pdb_copy = copy.deepcopy(pdb_structure)
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            pdb_copy.atom_to_internal_coordinates(verbose=args.v)
        elif pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            pdb_copy.internal_to_atom_coordinates(verbose=args.v)
        else:
            print("tdup: no input")
            exit()
        pdb_copy = copy.deepcopy(pdb_structure)
        r = compare_residues(pdb_structure, pdb_copy, True)
        print(target, fileNo, r["report"])
        print(r)

    if args.tip:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            for chn in pdb_structure.get_chains():
                cic = chn.internal_coord
                cic.atomArrayValid[:] = False
                cic.atomArrayValid[0:3] = True
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
        elif pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)

    if args.wi:
        if pdb_input:
            # add_PIC(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
        write_PIC(
            pdb_structure, outfile + ".pic", picFlags=picFlags, hCut=hCut, pCut=pCut
        )
        print("wrote pic output for", outfile)

    if args.t or args.tv or args.tok:
        sp = StringIO()
        if pdb_input:
            with warnings.catch_warnings(record=True) as w:
                # warnings.simplefilter("error")
                r = {}
                if True:  # try:
                    if args.cp or args.ti:
                        doCpti(args, srt_string)
                    else:
                        if old_srt:
                            r = structure_rebuild_test(pdb_structure, args.tv)
                        else:
                            r = structure_rebuild_test(pdb_structure, args.tv, args.q)
                    warns = len(w) > 0
                    if args.tv and warns:
                        for wrn in w:
                            print(wrn.message)
                    if not (
                        args.q or args.ti or args.cp or (args.tok and not r["pass"])
                    ):
                        print(
                            target,
                            fileNo,
                            r["report"],
                            ("WARNINGS" if warns else ""),
                        )
                    elif args.q and not (args.tok and not r["pass"]):
                        print(
                            target,
                            fileNo,
                            f"OK: all {r['aCoordMatchCount']} atom coords match",
                        )

                # except Exception as e:
                #    print(prot_id, fileNo, "EXCEPTION:", type(e), e)

        elif pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            write_PDB(pdb_structure, sp)
            sp.seek(0)
            pdb2 = PDB_parser.get_structure(prot_id, sp)
            pdb2.atom_to_internal_coordinates(verbose=args.v)
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
            doCpti(
                args,
                "pdb_structure.atom_to_internal_coordinates(verbose=args.v)",
                pdb_structure,
            )

    if args.t6p or args.t6a or args.t6s or args.t6c or args.tt:
        pdb_copy = None
        if pic_input:
            # pdb_copy = read_PIC(
            #    gzip.open(filename, mode="rt") if filename.endswith(".gz") else filename
            # )
            # pdb_copy = copy.deepcopy(pdb_structure)
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            # pdb_copy.internal_to_atom_coordinates(verbose=args.v)
            pdb_copy = copy.deepcopy(pdb_structure)
        elif pdb_input:
            pdb_copy = copy.deepcopy(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            pdb_copy.atom_to_internal_coordinates(verbose=args.v)
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
        print(target, fileNo, r["report"])

    if args.t33:
        resTarg = 3  # 13  # 32  # account for 0-start child_list
        delta = 33  # degrees to change
        tdelta = delta  # / 10.0  # more realistic for bond angle
        pdb_copy = None
        if pic_input:
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            # pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            pdb_copy = copy.deepcopy(pdb_structure)
        elif pdb_input:
            pdb_copy = copy.deepcopy(pdb_structure)
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
            pdb_structure.internal_to_atom_coordinates(verbose=args.v)
            pdb_copy.atom_to_internal_coordinates(verbose=args.v)
        alist = [
            "omg",
            "phi",
            "psi",
            "N:CA:C:O",
            "chi1",
            "chi2",
            "chi3",
            "chi4",
            "chi5",
            "tau",
            "O:C:CA:CB",
            "-1C:N:CA:CB",
        ]
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
                    print(f"{ang if ang is not None else '  -  ':8.4}", end=" ")
                print()
                print("target:  ", end="")

                setter = True  # use edra.setter or direct array manipulation
                ang0 = None
                andx = None
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
                            cic0.atomArrayValid[cic0.h2aa[andx][2]] = False
                        else:
                            cic0.dihedraAngle[andx] += delta
                            if cic0.dihedraAngle[andx] > 180.0:
                                cic0.dihedraAngle[andx] -= 360.0
                            cic0.dihedraAngleRads[andx] = np.deg2rad(
                                cic0.dihedraAngle[andx]
                            )
                            cic0.dAtoms_needs_update[andx] = True
                            cic0.atomArrayValid[cic0.d2aa[andx][3]] = False
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
                cic0.internal_to_atom_coordinates(verbose=args.v)
                cic0.atom_to_internal_coordinates(verbose=args.v)

                # if a == "omg":
                #    write_PIC(pdb_structure, "omg1.pic")

                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")
                checkAllDihedra(a, ric, cic0, cic1, setter)
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
                cic0.internal_to_atom_coordinates(verbose=args.v)
                cic0.atom_to_internal_coordinates(verbose=args.v)

                print("finish:  ", end="")
                for ar in alist:
                    ang = ric.get_angle(ar)
                    tang = ric1.get_angle(ar)
                    print(color(ang, tang, (a == ar), 8.4, 0.001), end=" ")
                checkAllDihedra(a, ric, cic0, cic1, setter)
                print()

                print()

        r = compare_residues(pdb_structure, pdb_copy, True)
        # print(prot_id, fileNo, r["report"])
        pass

    if args.rama:
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
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
        # rtm update for bond_rotate / bond_set
        if pdb_input:
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
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
                    # ric.bond_set("tau", nv)  # test to raise exception
                    nvtau[str(r)] = nv

        pdb2 = IC_duplicate(pdb_structure)
        # write_PIC(pdb2, "foo.pic")
        pdb2.internal_to_atom_coordinates(verbose=args.v)
        sf = StringIO()
        write_PDB(pdb2, sf)
        sf.seek(0)
        new_pdb = PDB_parser.get_structure("1NEW", sf)
        new_pdb.atom_to_internal_coordinates(verbose=args.v)
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
            pdb_structure.atom_to_internal_coordinates(verbose=args.v)
        if args.flex:
            if args.gcb and pic_input:
                pdb_structure.internal_to_atom_coordinates(
                    verbose=args.v, start=start_offset, fin=fin_offset
                )
                pdb_structure.atom_to_internal_coordinates(
                    verbose=args.v
                )  # build C-beta's
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
            includeCode=(not args.nsc),
            maxPeptideBond=args.maxp,
            start=start_offset,
            fin=fin_offset,
        )
        # write_PIC(pdb_structure, "foo.pic")

    fileNo += 1

if args.hstat:
    print(f"# {len(hedra_stats)} total entries")
    print("hedra_defaults = {")
    print("    # [(averages), angle_sd]  # len1, len3 std dev, [total count]")
    for k, hdat in sorted(hedra_stats.items()):
        tot = len(hdat[0])
        avg = [0, 0, 0]
        sd = [0, 0, 0]
        for i in range(3):
            l12 = hdat[0]
            x = np.average(l12)
            avg[i] = np.average(hdat[i])
            sd[i] = np.std(hdat[i])
        print(
            f"    '{str(k)}' : [({avg[0]:7.5f}, {avg[1]:9.5f}, {avg[2]:7.5f}), {sd[1]:.5f}],  # {sd[0]:.5f}  {sd[2]:.5f} [{tot}]"
        )
    print("}\n")

if args.dstat:
    print(f"# {len(dihedra_stats)} total primary entries")
    print("dihedra_primary_defaults = {")

    for k, ddat in sorted(dihedra_stats.items()):
        tot = len(ddat[0])
        avg = Dihedron.angle_avg(ddat[0])
        sd = Dihedron.angle_pop_sd(ddat[0], avg)
        maxndx = np.where(ddat[1] == np.amax(ddat[1]))
        maxn = maxndx[0][0]
        maxn = maxn - 180
        print(f'    "{k}": [{maxn}, {sd:9.5f}],  # {tot} {avg:9.5f}')
    print("}\n")

    c = 0
    print("dihedra_secondary_defaults = {")
    print("    # primary angle to rotate from, average rotation # count std dev")
    for k, ddat in sorted(dihedra_stats.items()):
        for k2, ddat2 in ddat[2].items():
            if ddat[3][3] == "OXT" or k2[0] == "X":
                continue
            c += 1
            tot2 = len(ddat2[0])
            avg2 = Dihedron.angle_avg(ddat2[0])
            sd2 = Dihedron.angle_pop_sd(ddat2[0], avg2)
            print(f'    "{k2}": [{ddat[3]}, {avg2:9.5f}],  # {tot2} {sd2:9.5f} ')
    print("}" + f"  # {c} total secondary default entries\n")

    c = 0
    print("dihedra_secondary_xoxt_defaults = {")
    print("    # primary angle to rotate from, average rotation # count std dev")
    for k, ddat in sorted(dihedra_stats.items()):
        for k2, ddat2 in ddat[2].items():
            if ddat[3][3] != "OXT" and k2[0] != "X":
                continue
            c += 1
            tot2 = len(ddat2[0])
            avg2 = Dihedron.angle_avg(ddat2[0])
            sd2 = Dihedron.angle_pop_sd(ddat2[0], avg2)
            print(f'    "{k2}": [{ddat[3]}, {avg2:9.5f}],  # {tot2}  {sd2:9.5f} ')
    print("}" + f"  # {c} total secondary X / OXT default entries\n")


if args.pri:
    print(primary_angles)
    # for a in sorted(primary_angles):
    #    print(a)

if args.v:
    print("normal termination")
# print(AtomKey.icd)
