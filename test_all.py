#!/usr/bin/env python3

import sys
import gzip

# import warnings
import os
import re
import argparse

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

PDB_parser = PDBParser(PERMISSIVE=False, QUIET=False)
PDB_repository_base = None
pdbDirs = [
    "/Users/rob/data/pdb/",
    "/media/data/pdb/",
    "/Volumes/data/pdb/",
    "/media/data/structures/divided/pdb/",
]
for d in pdbDirs:
    if os.path.isdir(d):
        PDB_repository_base = d
        break

if PDB_repository_base is None:
    print("failed to set PDB_repoistory_base.")

CIF_parser = MMCIFParser(QUIET=True)
mmCIF_repository_base = None
cifDirs = [
    "/Users/rob/data/mmCIF/",
    "/media/data/mmCIF/",
    "/Volumes/data/mmCIF/",
    "/media/data/structures/divided/mmCIF/",
]
for d in cifDirs:
    if os.path.isdir(d):
        mmCIF_repository_base = d
        break

if mmCIF_repository_base is None:
    print("failed to set mmCIF_repoistory_base.")

pdbidre = re.compile(r"(^\d(\w\w)\w)(\w)?$")

arg_parser = argparse.ArgumentParser(description="Test parsing sructure files.")

arg_parser.add_argument(
    "-skip", dest="skip_count", help="count of pdb ID list entries to skip"
)
arg_parser.add_argument(
    "-limit",
    dest="limit_count",
    help="stop after processing this many pdb ID list entries",
)
arg_parser.add_argument(
    "-cif", help="get mmCIF file for PDB id (default pdb file)", action="store_true"
)

args = arg_parser.parse_args()
if args.skip_count:
    args.skip_count = int(args.skip_count)
if args.limit_count:
    args.limit_count = int(args.limit_count)

fileNo = 0

if not sys.stdin.isatty():
    pass
else:
    print("no data on stdin, waiting for input...")

for line in sys.stdin:  # ["7rsa", "2aty"]:  #
    line = line.rstrip()
    if "q" == line:
        break
    fileNo += 1
    if args.skip_count and fileNo <= args.skip_count:
        continue
    if args.limit_count is not None:
        if args.limit_count <= 0:
            # sys.exit(0)
            break
        args.limit_count -= 1

    fields = line.split()
    pid = ""
    m = pdbidre.match(fields[0])
    if m:
        if args.cif:
            fil = mmCIF_repository_base + m.group(2) + "/" + m.group(1) + ".cif.gz"
        else:
            fil = PDB_repository_base + m.group(2) + "/pdb" + m.group(1) + ".ent.gz"
    else:
        fil = fields[0]
        # if not fil.endswith(".gz"):  # hack
        #    fil += ".gz"

    ofil = fil

    if fil.endswith(".gz"):
        try:
            f = gzip.open(fil, mode="rt")
        except Exception as e:
            print(ofil, fileNo, "gzip open EXCEPTION:", type(e), e)
            continue
        fil = fil[:-3]
    else:
        f = fields[0]
    pid = fil[-8:-4]
    ptype = ""
    if fil.endswith(".ent"):
        ptype = "pdb"
        try:
            struct = PDB_parser.get_structure(pid, f)
        except Exception as e:
            print(ofil, pid, fileNo, ptype, "EXCEPTION:", type(e), e)
            continue
    elif fil.endswith(".cif"):
        ptype = "cif"
        try:
            struct = CIF_parser.get_structure(pid, f)
        except Exception as e:
            print(ofil, pid, fileNo, ptype, "EXCEPTION:", type(e), e)
            continue
    else:
        print("skip:", line)
        continue

    if struct.header["resolution"]:
        res = struct.header["resolution"]
    else:
        res = -1.0

    chnCount = 0
    for c in struct.get_chains():
        chnCount += 1
        rcount = 0  # residue count
        acount = 0  # atom count
        drcount = 0  # disordered residue count
        dacount = 0  # disordered atom count
        mrcount = 0  # missing residue count
        nscount = 0  # non-sequenctial count
        seqpos = None
        for r in c.get_residues():
            rcount += 1
            prev = int(r.id[1]) - 1
            if (
                seqpos is not None
                and seqpos != prev
                and seqpos < r.id[1]
                and " " == r.id[0]
                and " " == r.id[2]
            ):
                mrcount += r.id[1] - seqpos
            elif (
                seqpos is not None
                and seqpos > r.id[1]
                and " " == r.id[0]
                and " " == r.id[2]
            ):
                print(r)
                nscount += 1

            seqpos = int(r.id[1])

            if r.is_disordered():
                drcount += 1
        for a in c.get_atoms():
            acount += 1
            if a.is_disordered():
                dacount += 1

        print(
            fileNo,
            ptype,
            pid,
            c.id,
            chnCount,
            rcount,
            drcount,
            acount,
            dacount,
            res,
            mrcount,
            nscount,
        )
