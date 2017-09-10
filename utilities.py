# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:32:04 2017

@author: Benjamin
"""

import math
import numpy as np
from collections import Counter

# -----------------------------------------------------------------------------

def remove_all(input_list, elem_to_remove):
    """
    Returns a copy of "input_list" in which all instances of "elem_to_remove"
    have been removed.
    """
    idx_to_keep = []
    for idx, elem in enumerate(input_list):
        if not elem == elem_to_remove:
            idx_to_keep.append(idx)
    # Build new list "return_list":
    return_list = []
    for idx in idx_to_keep:
        return_list = return_list + list([input_list[idx]])

    return return_list

def count_atoms(chains):
    num_atoms = 0
    sorted_atom_idxs = []
    for gen in chains:
        for bn in gen:
            for atom_idx in bn:
                if atom_idx not in sorted_atom_idxs:
                    sorted_atom_idxs.append(atom_idx)
                    num_atoms += 1
                else:
                    print("In count_atoms: duplicate atom found.")
    return num_atoms

def print_generations(bin_generations):
    for gen_idx, gen in enumerate(bin_generations):
        #print("gen_idx: " + str(gen_idx))
        #print(gen)
        line_str = ""
        for bn in gen:
            #print("bn: " + str(bn))
            line_str = line_str + str([x for x in bn])
            line_str = line_str + " -- "
        print(line_str)





