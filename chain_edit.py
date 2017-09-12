# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:39:13 2017

@author: Benjamin
"""


import numpy as np
from copy import deepcopy

from bin_sort import bins_are_neighbours

def remove_duplicates_from_all_tips(chains):
    tips = deepcopy(chains[-1])

    existing_atoms_idxs = []

    for chain_idx in range(len(tips)):
        for atom_idx in chains[-1][chain_idx]:
            if atom_idx not in existing_atoms_idxs:
                existing_atoms_idxs.append(atom_idx)
            else:
                bn = chains[-1][chain_idx]
                index = np.argwhere(bn==atom_idx)[0][0]
                print("atom_idx to remove = " + str(atom_idx+1))
                print("index to remove = " + str(index))
                print("bn before remove: " + str([x+1 for x in bn]))
                chains[-1][chain_idx] = np.delete(bn, index)
                print("bn after remove: " + str([x+1 for x in chains[-1][chain_idx]]))

    print("existing_atoms_idxs: " + str(existing_atoms_idxs))


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


def glue_chains(chain1, chain2, interact_mtrx):
    """
    Glue dangling ends of two chains together:
    [0 --> chain1  -->  -1]  +  [-1  -->  chain2  -->  0]
    """
    # print("In glue_chains: \n")
    new_chain = []

    for gen_idx, bn in enumerate(chain1):
        new_chain.append(deepcopy(bn))

    # print("new_chain: "+ str([x+1 for x in new_chain]))

    for ii in range(1, len(chain2)+1):
        idx = len(chain2) - ii
        bn = deepcopy(chain2[idx])

        # print("bn" + str([x+1 for x in bn]))
        if ii == 1:
            print("chain1[-2]: " + str([x+1 for x in chain1[-2]]))
            print("chain2[idx-1]: " + str([x+1 for x in chain2[idx-1]]))

            if bins_are_neighbours(chain1[-1], chain2[idx-1], interact_mtrx) \
                or bins_are_neighbours(chain1[-2], chain2[idx], interact_mtrx):
                print("Merge ends.")
                bn = np.append(bn, deepcopy(new_chain[-1]))
                # print("bn after merge" + str([x+1 for x in bn]))
                new_chain = new_chain[0:-1]
            else:
                print("Don't merge ends.")

        new_chain.append(bn)

    return new_chain