#!/usr/local/bin/python3
# Copyright (c) 2021 Robert T. Miller.  All rights reserved.

# -*- coding: latin-1 -*-
"""Collect potential interacting PHEs using internal coordinates module.

Proof of concept only.
"""

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.internal_coords import AtomKey  # IC_Chain

from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBIO import PDBIO

pdbSet = [
    "3pbl",
    "3rze",
]

pdbl = PDBList(verbose=False)
parser = MMCIFParser()
io = PDBIO()

chainPerGroup = True  # groups have same chain ID else all center have same ID
chnId = "A"

for pdb in pdbSet:
    # get first chain and generate coordinate space transforms
    pdbl.retrieve_pdb_file(pdb, file_format="mmCif")
    fil = pdb[1:3] + "/" + pdb + ".cif"
    chnA = parser.get_structure(pdb, fil)[0]["A"]
    chnA.atom_to_internal_coordinates()
    chnA.internal_to_atom_coordinates()

    # variables to save typing
    cic = chnA.internal_coord
    atomArrayNdx = cic.atomArrayIndex
    resNameNdx = AtomKey.fields.resname

    # numpy selector for atoms in PHE residues
    pheAtomSelect = [
        atomArrayNdx.get(k) for k in atomArrayNdx.keys() if k.akl[resNameNdx] == "F"
    ]
    # get just the Phe atoms (fancy indexing creates new array)
    aaF = cic.atomArray[pheAtomSelect]

    lastPDB = ""
    io.set_structure(chnA)
    io.atom_number = 1
    done = []

    for res0 in chnA.get_residues():
        # outer loop, transform res0 to origin
        if res0.resname in ["PHE"]:
            groupList = [res0]  # to avoid duplicating interacting groups
            ric0 = res0.internal_coord

            # chi1 = ric0.pick_angle("chi1")  # chi1 space defined with CA at origin
            chi1 = ric0.pick_angle("N:CA:C:O")  # chi1 space defined with CA at origin
            cst = np.transpose(chi1.cst)  # transform TO chi1 space
            rcst = np.transpose(chi1.rcst)  # transform FROM chi1 space

            cic.atomArray[pheAtomSelect] = aaF.dot(cst)  # transform just the PHEs

            # announce source PDB file
            if pdb != lastPDB:
                print("REMARK 500 SOURCE " + pdb)
                lastPDB = pdb

            partners = ""
            for res1 in chnA.get_residues():
                # inner loop, test other PHEs for proximity to res0 at origin
                if res1.resname in ["PHE"]:
                    ric1 = res1.internal_coord
                    if ric1 != ric0:
                        caDistance = np.linalg.norm(res1["CA"].coord, 2)
                        if caDistance < 6.0:
                            # res1 CA is < 6 ang from res0 CA (res0 CA at 0,0,0)
                            partners += (
                                f"REMARK 500 CA DISTANCE {caDistance:.2f} ANGSTROMS\n"
                            )
                            cid = chnId if chainPerGroup else chr(ord(chnId[0]) + 1)
                            # partners += io.pdb_residue_string(res1, cid)
                            # partners += io.get_ter_str()
                            partners += ric1.pdb_residue_string()
                            groupList.append(res1)

            if partners != "":  # if found some close residues
                groupSet = set(groupList)
                if groupSet not in done:  # skip same group wih different order
                    # now print group in PDB format
                    print("REMARK 500 CENTER " + str(res0.full_id) + " " + res0.resname)
                    # print(io.pdb_residue_string(res0, chnId) + io.get_ter_str(), end="")
                    print(ric0.pdb_residue_string())
                    print(partners, end="")

                    if chainPerGroup:
                        # each group has same chain ID so get new one for next group
                        chnId = chr(ord(chnId[0]) + 1)
                    # remember did this group
                    done.append(groupSet)
                    # improvement: 3rze 432-434-435 should center on 435

            # return PHEs to original coordinate space
            cic.atomArray[pheAtomSelect] = aaF.dot(rcst)
