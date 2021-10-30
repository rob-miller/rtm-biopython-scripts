#!/usr/local/bin/python3
# Copyright (c) 2021 Robert T. Miller.  All rights reserved.

# -*- coding: latin-1 -*-
"""Generate C-alpha plot with internal coordinates module.

Proof of concept only.
"""
import numpy as np
import matplotlib.pyplot as plt

from Bio.PDB import MMCIFParser
from Bio.PDB import internal_coords

# get the file
# from Bio.PDB.PDBList import PDBList
# pdbl = PDBList()
# pdbl.retrieve_pdb_file("7RSA")

# read structure of ribonuclease A
parser = MMCIFParser()
pdb7rsaA = parser.get_structure("7rsa", "./rs/7rsa.cif")[0]["A"]

# compute internal coordinates
pdb7rsaA.atom_to_internal_coordinates()

# references to data structures to save typing:
atomArrayAll = pdb7rsaA.internal_coord.atomArray
atomArrayNdx = pdb7rsaA.internal_coord.atomArrayIndex
atmNameNdx = internal_coords.AtomKey.fields.atm

# generate list of atomArray indexes for only C-alpha atoms
CaSelection = [
    atomArrayNdx.get(k) for k in atomArrayNdx.keys() if k.akl[atmNameNdx] == "CA"
]

# numpy fancy indexing to get just the C-alpha atoms
atomArrayCa = atomArrayAll[CaSelection]

# create the distance matrix (see also scipy)
# from https://jbencook.com/pairwise-distance-in-numpy/
caDist = np.linalg.norm(atomArrayCa[:, None, :] - atomArrayCa[None, :, :], axis=-1)

# show the heatmap
# https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib
plt.imshow(caDist, cmap="hot", interpolation="nearest")
plt.show()
