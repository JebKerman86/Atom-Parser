# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:39:13 2017

@author: Benjamin
"""


import numpy as np
from copy import deepcopy




def merge(chains, col_gen_idx, chain1_idx, chain2_idx):
    """
    Merge chain2 into chain1, starting at col_gen_idx
    """
    merged_chains = []
    for gen in chains:
        merged_chains.append(deepcopy(gen))

    # merged_chains.append(col_gen)

    for gen_idx in range(len(merged_chains)):
        # print("gen_idx " + str(gen_idx))
        merged_bin = np.concatenate(
            (deepcopy(merged_chains[gen_idx][chain1_idx]),
                 deepcopy(merged_chains[gen_idx][chain2_idx])),
                 axis=0)
        merged_bin = np.unique(np.array(merged_bin))
        # merged_bin = [int(x) for x in merged_bin]
        print("merged_bin: ")
        print([x+1 for x in merged_bin])
        merged_chains[gen_idx][chain1_idx] = merged_bin
        merged_chains[gen_idx][chain2_idx] = np.array([])

    return merged_chains