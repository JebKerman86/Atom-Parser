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
                     glue_chains, merge, remove_duplicates_from_all_tips, \
                     remove_duplicates_from_tips, get_chain_length

from utilities import print_generations, count_atoms, \
                      print_final_chain, test_solution

# Name without file ending:
# INPUT_FILE_NAME = "caffeine"
# INPUT_FILE_NAME = "caffeine_no_simultaneous_collision"
# INPUT_FILE_NAME = "1d_kette"
# INPUT_FILE_NAME = "t-kreuzung_sackgasse"
# INPUT_FILE_NAME = "t-kreuzung_dick"
# INPUT_FILE_NAME = "t-kreuzung_langer_arm"
# INPUT_FILE_NAME = "kompliziert"
# INPUT_FILE_NAME = "zno2wire"
# INPUT_FILE_NAME = "SiNW"

ALL_TEST_FILE_NAMES = ["1d_kette", "zno2wire", "t-kreuzung_dick", "t-kreuzung_langer_arm", "t-kreuzung_sackgasse", "caffeine_no_simultaneous_collision", "caffeine", "kompliziert"]
TEST_FILE_NAMES = ALL_TEST_FILE_NAMES[0:]

LOAD_CACHE_DATA = False
OPEN_JMOL = False

# Maximal number of Generations, this is a maximum value for safety, to
# protect the program from getting stuck in an infinite loop.
# Increase if necessary.
MAX_GENERATIONS = 100


def main():

    """
    Entry Point for Program
    """
    is_solution_list = []
    
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
        gen_idx_of_last_generation = -1
    
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
                                    if not num_chains == 2:
                                        print("WEIRD PROBLEM: num_chains should be 2")
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
                # remove_duplicates_from_tips(chains, target_chain_idx, src_chain_idx)
                remove_duplicates_from_all_tips(chains)
                chains = merge(chains, curr_gen_idx, target_chain_idx, src_chain_idx)
                print("\n Chains after merge step: ")
                print_generations(chains)
    
            print("remove_duplicates_from_ALL_tips")
            remove_duplicates_from_all_tips(chains)
            num_sorted_atoms = count_atoms(chains) + num_unlisted_contact_atoms
            print("num_sorted_atoms: " + str(num_sorted_atoms))
            if num_sorted_atoms >= num_atoms:
                # Fehler ausgeben falls größer!!
                print("All atoms sorted.")
                gen_idx_of_last_generation = curr_gen_idx
                break
            curr_gen_idx += 1
    
        if curr_gen_idx >= MAX_GENERATIONS:
            sys.exit("FATAL ERROR: MAX_GENERATIONS exceeded. (Increase MAX_GENERATIONS?)")
    
        if not final_collision_found:
            sys.exit("FATAL ERROR: No final collision found, don't know which chains to keep")
        
        print("\n Chain before culling dead ends: ")
        print_generations(chains)
        
        # length = get_chain_length(chains, final_chain_idxs[0], gen_idx_of_last_collision+1)
        # print("length = " + str(length))
        
        num_gen = len(chains)
        
        #Find dead ends in the two final chains
        dead_end_tuples = []
        dead_ends = []
        dead_end_start_idx = gen_idx_of_last_collision+1
        for chain_idx in final_chain_idxs:
            print("chain_idx = " + str(chain_idx))
            print("dead_end_start_idx = " + str(dead_end_start_idx))
            length = get_chain_length(chains, chain_idx, dead_end_start_idx)
            print("length = " + str(length))
            if length == 0:
                dead_end_tuples.append((-1,-1))
                dead_ends.append([])
            else:
                dead_end_tip_idx = dead_end_start_idx + length
                dead_end_tuples.append((dead_end_start_idx, dead_end_tip_idx-1))
                dead_end = []
                for gen in chains[dead_end_start_idx:]:
                    print("gen: " + str([x+1 for x in gen]))
                    dead_end.append(gen[chain_idx])
                dead_ends.append(deepcopy(dead_end))
    
        print("dead_end_tuples: " + str(dead_end_tuples))
    
        for dead_end in dead_ends:
            print("dead_end: " + str([x+1 for x in dead_end]))
    
        # Before Merging dead ends, we have to make sure the dead end isn't longer
        # than the final chain we are attempting to merge it into
    
        chain_length_until_last_collision = gen_idx_of_last_collision+1
        shortened_dead_ends = []
        print("Shorten dead ends:")
        for dead_end in dead_ends:
            print("dead_end: " + str([x+1 for x in dead_end]))
            shortened_dead_end = dead_end
            dead_end_length = len(dead_end)
            if dead_end_length > 0:
                print("dead_end_length > 0")
                # subtract "1", because otherwise we would merge into the contact
                while dead_end_length > chain_length_until_last_collision-1:
                    print("chain_length_until_last_collision = " + str(chain_length_until_last_collision))
                    shortened_dead_end = []
                    # Shorten dead end to make it fit
                    print("dead_end_length: " + str(dead_end_length))
                    for bn_idx, bn in enumerate(dead_end):
                        merged_bn = np.array([])
                        if dead_end_length%2 == 1 and bn == dead_end[-1]:
                            print("dead_end_length%2 == 1 and bn == dead_end[-1]")
                            merged_bn = np.append(merged_bn, dead_end[-1])
                            merged_bn = [int(x) for x in merged_bn]
                            shortened_dead_end.append(deepcopy(merged_bn))
                        if bn_idx%2 == 1:
                            print("bn_idx%2 == 1")
                            merged_bn = np.append(merged_bn, np.array([dead_end[bn_idx-1], dead_end[bn_idx]]))
                            merged_bn = [int(x) for x in merged_bn]
                            print("merged_bn: " + str(merged_bn))
                            shortened_dead_end.append(deepcopy(merged_bn))
                    dead_end_length = len(shortened_dead_end)
                    dead_end = shortened_dead_end
                    print("dead_end_length: " + str(dead_end_length))
    
            shortened_dead_ends.append(shortened_dead_end)
            
    
    
        print("\nMerge dead ends: \n")
        for idx, dead_end in enumerate(shortened_dead_ends):
            chain_idx = final_chain_idxs[idx]
            print("chain_idx: " + str(chain_idx))
            print("dead_end: " + str(dead_end))
            other_chain_idx = final_chain_idxs[(chain_idx+1)%2]
            print("other_chain_idx: " + str(other_chain_idx))
            dead_end_length = num_gen - chain_length_until_last_collision
            print("dead_end_length = " + str(dead_end_length))
            print("deleting end of chain: " + str(chain_idx))
            for gen_idx in range(dead_end_length):
                print("delete bin of generation index: " + str(chain_length_until_last_collision+gen_idx))
                chains[chain_length_until_last_collision+gen_idx][chain_idx] = np.array([])
    
            for gen_idx, bn in enumerate(dead_end):
                print("bn: " + str(bn))
                print("chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]: " + str([x+1 for x in chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]]))
                merged_bn = np.append(chains[gen_idx_of_last_collision-gen_idx][other_chain_idx],bn)
                chains[gen_idx_of_last_collision-gen_idx][other_chain_idx] = merged_bn
                print("chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]: " + str([x+1 for x in chains[gen_idx_of_last_collision-gen_idx][other_chain_idx]]))
    
    
        """
    
        step = 0
        for gen_idx in range(gen_idx_of_last_collision+1, num_gen):
            # print("gen_idx = " + str(gen_idx))
    
            for idx, src_chain_idx in enumerate(final_chain_idxs):
                # This should always yield the other index
                target_chain_idx = final_chain_idxs[(idx+1)%2]
                # print("chain_idx = " + str(src_chain_idx))
                # print("target_chain_idx = " + str(target_chain_idx))
                src_bin = chains[gen_idx][src_chain_idx]
                target_bin = chains[gen_idx_of_last_collision-step][target_chain_idx]
                # print("target_bin = " + str(target_bin))
                # print("src_bin = " + str(src_bin))
                target_bin = np.append(target_bin, deepcopy(src_bin))
                target_bin = np.array([int(x) for x in target_bin])
                chains[gen_idx_of_last_collision-step][target_chain_idx] = target_bin
                chains[gen_idx][src_chain_idx] = np.array([])
            step += 1
            
        """
    
        remove_duplicates_from_all_tips(chains)
    
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
    
        print("\nfinal_chain: ")
        print_final_chain(final_chain)
    
        is_solution = test_solution(final_chain, interact_mtrx)
        print("is_solution: " + str(is_solution))
        write_bins(final_chain, atom_positions, INPUT_FILE_NAME, OPEN_JMOL)
        
        is_solution_list.append(is_solution)
        
        
    print("- INPUT_FILE_NAME ---------------- solution found:")
    for idx, INPUT_FILE_NAME in enumerate(TEST_FILE_NAMES):
        # print(INPUT_FILE_NAME + ": " + str(is_solution_list[idx]))
        print("%-*s  success: %s" % (35,INPUT_FILE_NAME,str(is_solution_list[idx])))
    
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
