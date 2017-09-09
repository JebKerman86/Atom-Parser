# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import numpy as np
from copy import deepcopy

from prep_input import prep_data
from file_io import chache_data, load_data, write_bins, read_xyz_file, read_transport_file
from bin_sort import get_contact_bins, get_next_bins, remove_common_elems, \
                     find_all_collisions, merge_chain, glue_chains, bins_are_neighbours, merge, remove_duplicates_from_tips, test_solution
from utilities import remove_all, print_generations

LOAD_CACHE_DATA = False
# Name without file ending:
# INPUT_FILE_NAME = "caffeine"
INPUT_FILE_NAME = "caffeine_no_simultaneous_collision"
# INPUT_FILE_NAME = "1d_kette"
# INPUT_FILE_NAME = "t-kreuzung_sackgasse"
# INPUT_FILE_NAME = "t-kreuzung_dick"
# INPUT_FILE_NAME = "zno2wire"
# INPUT_FILE_NAME = "SiNW"
OPEN_JMOL = True

MAX_GENERATIONS = 100

def main():

    """
    Entry Point for Program
    """

    # The loaded data are NOT numpy arrays (change later?)
    (atom_types, atom_positions) = read_xyz_file(str(INPUT_FILE_NAME))
    (region_list, interaction_distances) = read_transport_file(str(INPUT_FILE_NAME))

    if LOAD_CACHE_DATA:
        data = load_data(INPUT_FILE_NAME)
    else:
        data = prep_data(atom_types, atom_positions, region_list, interaction_distances)
        chache_data(INPUT_FILE_NAME, data)

    dist_mtrx = data["dist_mtrx"]
    interact_mtrx = data["interact_mtrx"]
    ordered_idx_mtrx = data["ordered_idx_mtrx"]
    ordered_dist_mtrx = data["ordered_dist_mtrx"]
    ordered_interact_mtrx = data["ordered_interact_mtrx"]

    num_atoms = np.size(dist_mtrx, axis=0)

    # Don't convert to np array in file_io because this complictes caching
    # (since numpy arrays non json-serialzable)

    np_region_list = []
    for region in region_list:
        np_region_list.append(np.array(region))


    device = np_region_list[0]
    contacts = np_region_list[1:]
    print("contacts: " + str(contacts))
    contact_bins = get_contact_bins(device, contacts, interact_mtrx)
    prev_bins = list.copy(contacts)
    #print("prev_bins" + str(prev_bins))

    #Each element in "bin_generations" is a list of bins. Each of these lists
    #contains the bins of a specific generation. The bins are sorted in the
    #order of ascending contact indices.
    #All bins that are the same number of steps away from the contacts are
    #assigned to the same "generation", the atoms in "contact_bins" are
    #in generation zero.
    #"contact_bins": All contact atoms that are interacting with the device
    # are assigned to this bin.

    print("contact_bins: " + str(contact_bins))
    chains = []
    chains.append(contact_bins)
    # print("bin_generations" + str(bin_generations))
    num_chains = len(contacts)
    """
    active_chains = []
    for c in range(num_chains):
        active_chains.append(True)
    print("active_chains: " + str(active_chains))
    """
    print("num_chains: " + str(num_chains))
    curr_gen_idx = 1

    # contact atoms count as sorted from the beginning
    num_sorted_atoms = 0
    for c in contacts:
        num_sorted_atoms += len(c)
        
    final_collision_found = False
    final_chain_idxs = []
    gen_idx_of_last_collision = -1

    while curr_gen_idx < MAX_GENERATIONS:
        print("curr_gen_idx: " + str(curr_gen_idx))
        curr_gen = get_next_bins(chains[-1], prev_bins, interact_mtrx)

        chains.append(curr_gen)
        prev_bins = prev_bins + curr_gen
        
        print("\n Chains before merge step ")
        print_generations(chains)

        # PROBLEM: RIGHT NOW WE CAN ONLY HANDLE ONE COLLISION PER GENERATION
        if not final_collision_found:
            collision_found = False
            for chain1_idx, bn1 in enumerate(curr_gen):
                for chain2_idx, bn2 in enumerate(curr_gen):
                    if chain2_idx > chain1_idx:
                        if bins_are_neighbours(bn1, bn2, interact_mtrx):
                            if num_chains > 2:
                                # Merge chains
                                print("bn1: " +str([x+1 for x in bn1]))
                                print("bn2: " +str([x+1 for x in bn2]))
                                chains = merge(chains, curr_gen_idx, chain1_idx, chain2_idx)
                                num_chains -= 1
                                
                                remove_duplicates_from_tips(chains)
                                
                                # CHECK IF OTHER CHAIN TIP CONTAINS ATOMS FROM MERGED CHAIN TIP
                                # IF YES: DELETE ATOMS TO AVOID DUPLICATES
                                """
                                all_chain_idxs = list(range(len(curr_gen)))
                                print("all_chain_idxs" + str(all_chain_idxs))
                                print("chain1_idx = " + str(chain1_idx))
                                print("chain2_idx = " + str(chain2_idx))
                                merged_chain_idx = chain1_idx
                                not_merged_chain_idxs = [x for x in all_chain_idxs if not x in [chain1_idx, chain2_idx]]
                                print("merged_chain_idx: " + str(merged_chain_idx))
                                merged_chain_tip = chains[curr_gen_idx][merged_chain_idx]
                                print("not_merged_chain_idxs" + str(not_merged_chain_idxs))
                                for not_merged_chain_idx in not_merged_chain_idxs:
                                    bn = chains[curr_gen_idx][not_merged_chain_idx]
                                    bn = np.array([x for x in bn if not x in merged_chain_tip])
                                    chains[curr_gen_idx][not_merged_chain_idx] = bn
                                """
                                print("\n Chains after merge step: ")
                                print_generations(chains)
                                collision_found = True
                            else:
                                if not num_chains == 2:
                                    print("WEIRD PROBLEM: num_chains should be 2")
                                print("\n ---- final_collision_found! ---- \n")
                                final_collision_found = True
                                final_chain_idxs = [chain1_idx, chain2_idx]
                                gen_idx_of_last_collision = curr_gen_idx
                                remove_duplicates_from_tips(chains)

                    if collision_found or final_collision_found:
                        break
                if collision_found or final_collision_found:
                    break

        for bn in chains[-1]:
            num_sorted_atoms += len(np.unique(bn))

        if num_sorted_atoms >= num_atoms:
            print("num_sorted_atoms = " + str(num_sorted_atoms))
            break
        curr_gen_idx += 1

    # print("chains" + str(chains))
    print("\n Chain before culling dead ends: ")
    print_generations(chains)
    step = 0
    num_gen = len(chains)
    for gen_idx in range(gen_idx_of_last_collision+1, num_gen):
        print("gen_idx = " + str(gen_idx))

        for idx, src_chain_idx in enumerate(final_chain_idxs):
            target_chain_idx = final_chain_idxs[(idx+1)%2]
            print("chain_idx = " + str(src_chain_idx))
            print("target_chain_idx = " + str(target_chain_idx))
            target_bin = chains[gen_idx_of_last_collision-step][target_chain_idx]
            src_bin = chains[gen_idx][src_chain_idx]
            print("target_bin = " + str(target_bin))
            print("src_bin = " + str(src_bin))
            target_bin = np.append(target_bin, deepcopy(src_bin))
            target_bin = np.array([int(x) for x in target_bin])
            chains[gen_idx_of_last_collision-step][target_chain_idx] = target_bin
            chains[gen_idx][src_chain_idx] = np.array([])
        step += 1
        


    remove_duplicates_from_tips(chains)
    
    print("\n chain after removing duplicates from tips, and before glueing: ")
    print_generations(chains)

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


    print("\nfinal_chain: ")
    
    line_str = ""
    for bn in final_chain:
        line_str = line_str + str([x+1 for x in bn])
        line_str = line_str + "\n"

    print(line_str)

    is_solution = test_solution(final_chain, interact_mtrx)
    print("is_solution: " + str(is_solution))
    write_bins(final_chain, atom_positions, INPUT_FILE_NAME, OPEN_JMOL)

    """
    ToDo:
        Fix caffeine (problem: program breaks if there are two simultaneous collisions, since no final collision recognized)

    """



    """
          Algorithm: Simultaneously, from all contacts move into device.
                     Partition such that all new partition atoms are in
                     contact with previous partition
                     when partitions starting from two different contacts
                     collide: Merge partitions from the two branches into one.
                     
                     If more than two branches collide: Merge branches such
                     that the largest branch is as small as possible
                     
                     If there is still ambigouity, use contact indices of
                     branches to establich merge order
    """


if __name__ == "__main__":
    main()
