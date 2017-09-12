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
                # print("atom_idx to remove = " + str(atom_idx+1))
                # print("index to remove = " + str(index))
                # print("bn before remove: " + str([x+1 for x in bn]))
                chains[-1][chain_idx] = np.delete(bn, index)
                # print("bn after remove: " + str([x+1 for x in chains[-1][chain_idx]]))

    # print("existing_atoms_idxs: " + str(existing_atoms_idxs))


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
        # print("merged_bin: ")
        # print([x+1 for x in merged_bin])
        merged_chains[gen_idx][chain1_idx] = merged_bin
        merged_chains[gen_idx][chain2_idx] = np.array([])

    return merged_chains


def shorten_dead_ends(dead_ends, chain_length_until_last_collision):
    shortened_dead_ends = []

    # print("Shorten dead ends:")
    for dead_end in dead_ends:
        # print("dead_end: " + str([x+1 for x in dead_end]))
        shortened_dead_end = dead_end
        dead_end_length = len(dead_end)
        if dead_end_length > 0:
            # print("dead_end_length > 0")
            # subtract "1", because otherwise we would merge into the contact
            while dead_end_length > chain_length_until_last_collision-1:
                # print("chain_length_until_last_collision = " + str(chain_length_until_last_collision))
                shortened_dead_end = []
                # Shorten dead end to make it fit
                # print("dead_end_length: " + str(dead_end_length))
                for bn_idx, bn in enumerate(dead_end):
                    merged_bn = np.array([])
                    if dead_end_length%2 == 1 and bn == dead_end[-1]:
                        # print("dead_end_length%2 == 1 and bn == dead_end[-1]")
                        merged_bn = np.append(merged_bn, dead_end[-1])
                        merged_bn = [int(x) for x in merged_bn]
                        shortened_dead_end.append(deepcopy(merged_bn))
                    if bn_idx%2 == 1:
                        # print("bn_idx%2 == 1")
                        merged_bn = np.append(merged_bn, np.array([dead_end[bn_idx-1], dead_end[bn_idx]]))
                        merged_bn = [int(x) for x in merged_bn]
                        # print("merged_bn: " + str(merged_bn))
                        shortened_dead_end.append(deepcopy(merged_bn))
                dead_end_length = len(shortened_dead_end)
                dead_end = shortened_dead_end
                # print("dead_end_length: " + str(dead_end_length))

        shortened_dead_ends.append(shortened_dead_end)
        
    return shortened_dead_ends


def merge_dead_ends_into_final_chains(chains, shortened_dead_ends, final_chain_idxs, gen_idx_of_last_collision):
    
    num_gen = len(chains)
    chain_length_until_last_collision = gen_idx_of_last_collision+1
    # print("\nMerge dead ends: \n")
    for idx, dead_end in enumerate(shortened_dead_ends):
        chain_idx = final_chain_idxs[idx]
        # print("chain_idx: " + str(chain_idx))
        # print("dead_end: " + str(dead_end))
        other_chain_idx = final_chain_idxs[(chain_idx+1)%2]
        # print("other_chain_idx: " + str(other_chain_idx))
        dead_end_length = num_gen - chain_length_until_last_collision
        # print("dead_end_length = " + str(dead_end_length))
        # print("deleting end of chain: " + str(chain_idx))
        for gen_idx in range(dead_end_length):
            # print("delete bin of generation index: " + str(chain_length_until_last_collision+gen_idx))
            chains[chain_length_until_last_collision+gen_idx][chain_idx] = np.array([])

        for gen_idx, bn in enumerate(dead_end):
            # print("bn: " + str(bn))
            # print("chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]: " + str([x+1 for x in chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]]))
            merged_bn = np.append(chains[gen_idx_of_last_collision-gen_idx][other_chain_idx],bn)
            chains[gen_idx_of_last_collision-gen_idx][other_chain_idx] = merged_bn
            # print("chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]: " + str([x+1 for x in chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]]))



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
            # print("chain1[-2]: " + str([x+1 for x in chain1[-2]]))
            # print("chain2[idx-1]: " + str([x+1 for x in chain2[idx-1]]))

            if bins_are_neighbours(chain1[-1], chain2[idx-1], interact_mtrx) \
                or bins_are_neighbours(chain1[-2], chain2[idx], interact_mtrx):
                # print("Merge ends.")
                bn = np.append(bn, deepcopy(new_chain[-1]))
                # print("bn after merge" + str([x+1 for x in bn]))
                new_chain = new_chain[0:-1]
            else:
                # print("Don't merge ends.")
                pass

        new_chain.append(bn)

    return new_chain




def build_final_chain(chains, contacts, final_chain_idxs, interact_mtrx):
    
    final_chain1 = []
    for gen in chains:
        add_bin = gen[final_chain_idxs[0]]
        if not len(add_bin) == 0:
            final_chain1.append(add_bin)
        else:
            break

    cntct_bin1 = deepcopy(final_chain1[0])
    for atom_idx in contacts[final_chain_idxs[0]]:
        if atom_idx not in final_chain1[0]:
            cntct_bin1 = np.append(cntct_bin1, atom_idx)
    # brauche ich hier deepcopy???????
    final_chain1[0] = deepcopy(cntct_bin1)

    final_chain2 = []
    for gen in chains:
        add_bin = gen[final_chain_idxs[1]]
        if not len(add_bin) == 0:
            final_chain2.append(add_bin)
        else:
            break

    cntct_bin2 = deepcopy(final_chain2[0])
    for atom_idx in contacts[final_chain_idxs[1]]:
        if atom_idx not in final_chain2[0]:
            cntct_bin2 = np.append(cntct_bin2, atom_idx)
    final_chain2[0] = deepcopy(cntct_bin2)

    final_chain = glue_chains(final_chain1, final_chain2, interact_mtrx)
    
    return final_chain
    