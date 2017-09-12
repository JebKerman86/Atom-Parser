# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import sys
import numpy as np
from copy import deepcopy

from prep_input import prep_data

from file_io import chache_data, load_data, write_bins, read_xyz_file, \
                    read_transport_file

from bin_sort import get_contact_bins, get_next_bins, bins_are_neighbours, \
                     get_chain_length, get_dead_ends
                     
from chain_edit import remove_duplicates_from_all_tips, merge, shorten_dead_ends, \
                       merge_dead_ends_into_final_chains, build_final_chain

from utilities import print_generations, count_atoms, \
                      print_final_chain, test_solution

# Name without file ending:
ALL_TEST_FILE_NAMES = ["1d_kette", "zno2wire", "t-kreuzung_dick", "t-kreuzung_langer_arm", "t-kreuzung_sackgasse", "caffeine_no_simultaneous_collision", "caffeine", "kompliziert"]
TEST_FILE_NAMES = ALL_TEST_FILE_NAMES[0:]

LOAD_CACHE_DATA = False
DISPLAY_FILE_NAME = TEST_FILE_NAMES[0:0]

# Maximal number of Generations, this is a maximum value for safety, to
# protect the program from getting stuck in an infinite loop.
# Increase if necessary.
MAX_GENERATIONS = 100


def main():

    """
    Entry Point for Program
    """
    is_solution_list = []
    final_chains_list = []
    atom_positions_list = []
    
    for INPUT_FILE_NAME in TEST_FILE_NAMES:
        # The loaded data are NOT numpy arrays (change later?)
        (atom_types, atom_positions) = read_xyz_file(str(INPUT_FILE_NAME))
        (region_list, interaction_distances) = \
            read_transport_file(str(INPUT_FILE_NAME))
    
        if LOAD_CACHE_DATA:
            data = load_data(INPUT_FILE_NAME)
        else:
            data = prep_data(atom_types, atom_positions, region_list,
                             interaction_distances)
            chache_data(INPUT_FILE_NAME, data)
    
        atom_positions_list.append(atom_positions)
        dist_mtrx = data["dist_mtrx"]
        interact_mtrx = data["interact_mtrx"]
        # ordered_idx_mtrx = data["ordered_idx_mtrx"]
        # ordered_dist_mtrx = data["ordered_dist_mtrx"]
        # ordered_interact_mtrx = data["ordered_interact_mtrx"]
    
        num_atoms = np.size(dist_mtrx, axis=0)
    
        # Didn't convert to np array in file_io because this complicates caching
        # (since numpy arrays non json-serialzable)
        np_region_list = []
        for region in region_list:
            np_region_list.append(np.array(region))
    
        device = np_region_list[0]
        contacts = np_region_list[1:]
        print("contacts: " + str(contacts))
        contact_bins = get_contact_bins(device, contacts, interact_mtrx)
        num_unlisted_contact_atoms = \
            count_atoms([contacts]) - count_atoms([contact_bins])
        prev_bins = list.copy(contacts)

        # Each element in "chains" is a list of bins. Each of these lists
        # contains the bins of a specific generation. The bins are sorted in the
        # order of ascending contact indices.
        # All bins that are the same number of steps away from the contacts are
        # assigned to the same "generation", the atoms in "contact_bins" are
        # in generation zero.
        # "contact_bins": All contact atoms that are interacting with the device
        # are assigned to this bin.

        print("contact_bins: " + str(contact_bins))
        chains = []
        chains.append(contact_bins)
        # print("bin_generations" + str(bin_generations))
        num_chains = len(contacts)
        print("num_chains: " + str(num_chains))
        curr_gen_idx = 1

        final_collision_found = False
        final_chain_idxs = []
        gen_idx_of_last_collision = -1

        # This condition is a failsafe, to avoid infinite loops
        while curr_gen_idx < MAX_GENERATIONS:
            collisions_found = []
            print("curr_gen_idx: " + str(curr_gen_idx))
            curr_gen = get_next_bins(chains[-1], prev_bins, interact_mtrx)

            chains.append(curr_gen)
            prev_bins = prev_bins + curr_gen

            print("\n Chains before merge step ")
            print_generations(chains)

            if not final_collision_found:
                for chain1_idx, bn1 in enumerate(curr_gen):
                    for chain2_idx, bn2 in enumerate(curr_gen):
                        if chain2_idx > chain1_idx:
                            if bins_are_neighbours(bn1, bn2, interact_mtrx):
                                if num_chains > 2:
                                    collisions_found.append((chain1_idx,chain2_idx))
                                    print("collisions_found: " + str(collisions_found))
                                    num_chains -= 1
                                    print("num_chains = " + str(num_chains))

                                else:
                                    if num_chains < 2:
                                        sys.exit("FATAL ERROR: num_chains < 2")

                                    print("\n ---- final_collision_found! ---- \n")
                                    final_collision_found = True
                                    final_chain_idxs = [chain1_idx, chain2_idx]
                                    print("final_collision: " + str(final_chain_idxs))
                                    gen_idx_of_last_collision = curr_gen_idx
                                    print("gen_idx_of_last_collision = " + str(gen_idx_of_last_collision))
                                    print("remove_duplicates_from_ALL_tips")
                                    remove_duplicates_from_all_tips(chains)
                                    print("\n Chains after removing duplicates:")
                                    print_generations(chains)
    
                        if final_collision_found:
                            break
                    if final_collision_found:
                        break
    
            for col_tuple in collisions_found:
                # Merge from src_chain_idx into target_chain_idx
                src_chain_idx = col_tuple[0]
                target_chain_idx = col_tuple[1]
    
                if col_tuple[0] in final_chain_idxs:
                    if col_tuple[1] in final_chain_idxs:
                        sys.exit("FATAL ERROR: Should never merge the two \
                                 final chains into eachother.")
    
                # Make sure we are merging into the final chain.
                # If not, swap src_chain_idx with target_chain_idx.
                if col_tuple[0] in final_chain_idxs:
                    src_chain_idx = col_tuple[1]
                    target_chain_idx = col_tuple[0]
    
                # Merge chains
                print("Merge chain_idxs:" + str(col_tuple))
                print("src_chain bin: " +str([x+1 for x in curr_gen[src_chain_idx]]))
                print("target_chain bin: " +str([x+1 for x in curr_gen[target_chain_idx]]))
    
                ########################################
                # CONSIDER: FOR MULTIPLE COLLISIONS, TRY TO MERGE SMALLER CHAINS TOGETHER FIRST
                ########################################
                # Duplicates have to be removed AFTER collision recognition,
                # since otherwise this could prevent finding collisions
                print("remove_duplicates_from_ALL_tips")
                remove_duplicates_from_all_tips(chains)
                chains = merge(chains, curr_gen_idx, target_chain_idx, src_chain_idx)
                print("\n Chains after merge step: ")
                print_generations(chains)
    
            print("remove_duplicates_from_ALL_tips")
            remove_duplicates_from_all_tips(chains)
            num_sorted_atoms = count_atoms(chains) + num_unlisted_contact_atoms
            print("num_sorted_atoms: " + str(num_sorted_atoms))
            if num_sorted_atoms >= num_atoms:
                if num_sorted_atoms > num_atoms:
                    sys.exit("FATAL ERROR: num_sorted_atoms > num_atoms")
                print("All atoms sorted.")
                break
            curr_gen_idx += 1
    
        if curr_gen_idx >= MAX_GENERATIONS:
            sys.exit("FATAL ERROR: MAX_GENERATIONS exceeded. (Increase MAX_GENERATIONS?)")
    
        if not final_collision_found:
            sys.exit("FATAL ERROR: No final collision found, don't know which chains to keep")
        
        print("\n Chain before culling dead ends: ")
        print_generations(chains)
        
        #Find dead ends in the two final chains
        dead_ends = get_dead_ends(chains, final_chain_idxs, gen_idx_of_last_collision)
    
        # for dead_end in dead_ends:
        #     print("dead_end: " + str([x+1 for x in dead_end]))
    
        # Before Merging dead ends, we have to make sure the dead end isn't longer
        # than the final chain we are attempting to merge it into
        chain_length_until_last_collision = gen_idx_of_last_collision+1
        shortened_dead_ends = shorten_dead_ends(dead_ends, chain_length_until_last_collision)

        merge_dead_ends_into_final_chains(chains, shortened_dead_ends, final_chain_idxs, gen_idx_of_last_collision)

        remove_duplicates_from_all_tips(chains)

        print("\n chain after removing duplicates from tips, and before glueing: ")
        print_generations(chains)

        final_chain = build_final_chain(chains, contacts, final_chain_idxs, interact_mtrx)
        
        final_chains_list.append(final_chain)

        print("\nfinal_chain: ")
        print_final_chain(final_chain)
    
        is_solution = test_solution(final_chain, interact_mtrx)

        
        is_solution_list.append(is_solution)


    print("- INPUT_FILE_NAME ---------------- solution found:")
    for idx, INPUT_FILE_NAME in enumerate(TEST_FILE_NAMES):
        # print(INPUT_FILE_NAME + ": " + str(is_solution_list[idx]))
        print("%-*s  success: %s" % (35,INPUT_FILE_NAME,str(is_solution_list[idx])))


    OPEN_JMOL = []
    for INPUT_FILE_NAME in TEST_FILE_NAMES:
        
        if INPUT_FILE_NAME is DISPLAY_FILE_NAME:
            OPEN_JMOL.append(True)
        else:
            OPEN_JMOL.append(False)


    for idx, INPUT_FILE_NAME in enumerate(TEST_FILE_NAMES):
        write_bins(final_chains_list[idx], atom_positions_list[idx], INPUT_FILE_NAME, OPEN_JMOL[idx])
    
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
