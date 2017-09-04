# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import numpy as np

from prep_input import prep_data
from file_io import chache_data, load_data, write_bins, read_xyz_file, read_transport_file
from bin_sort import get_contact_bins, get_next_bins, remove_common_elems, \
                     find_all_collisions, merge_chain, glue_chains
from utilities import find_duplicates, remove_all, print_generations

LOAD_CACHE_DATA = False
# Name without file ending:
# INPUT_FILE_NAME = "caffeine"
INPUT_FILE_NAME = "t-kreuzung"
# INPUT_FILE_NAME = "t-kreuzung_sackgasse"
# INPUT_FILE_NAME = "zno2wire"
#INPUT_FILE_NAME = "SiNW"
OPEN_JMOL = False

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
    # print("contacts: " + str(contacts))
    contact_bins = get_contact_bins(device, contacts, interact_mtrx)
    added_bins = list.copy(contacts)
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
    bin_generations = []
    bin_generations.append(contact_bins)
    # print("bin_generations" + str(bin_generations))

    curr_generation = 1

    # contact atoms count as sorted from the beginning
    num_sorted_atoms = 0
    for c in contacts:
        num_sorted_atoms += len(c)

    while curr_generation < MAX_GENERATIONS:
        (curr_bins, prev_bins) = get_next_bins(bin_generations[-1], added_bins, interact_mtrx)
        # print("curr_generation = " + str(curr_generation))
        for b in curr_bins:
            num_sorted_atoms += len(b)
            # print("bin: " + str(b))
        added_bins = added_bins + curr_bins
        # print("prev_bins: " + str(prev_bins))
        bin_generations.append(curr_bins)
        (duplicates, common_elems) = find_duplicates(curr_bins)
        if duplicates > 0:
            # print("common_elems: " + str(common_elems))
            remove_common_elems(bin_generations[-1], common_elems)
            num_sorted_atoms -= duplicates

        if num_sorted_atoms >= num_atoms:
            break
        curr_generation += 1

    print("num_sorted_atoms: " + str(num_sorted_atoms))

    # contiguous_bin_generations = []
    # contiguous_bin_generations = create_subchain_tree(bin_generations, interact_mtrx)
    # print("contiguous_bin_generations: ")
    # print_generations(contiguous_bin_generations)
    print("bin_generations: ")
    print_generations(bin_generations)


    collision_list = find_all_collisions(bin_generations, interact_mtrx)
    print("collision_list: ")
    print(collision_list)

    
    # The last collision between chains is the collision between the longest
    # two chains. These are the chains we want to keep.
    last_col_list = collision_list[-1][1]
    if last_col_list == []:
        print("No collisions between chains detected (since last_col_list == [])")
        print("Something is weird. Exiting...")

    if len(last_col_list) > 1:
        print("Last collision is not unique. Don't know how to handle this yet.")
        print("Something is weird. Exiting...")


    last_col = last_col_list[0]
    keep_chains = (min(last_col),max(last_col))
    print("keep_chains: ")
    print(keep_chains)


    # Merge all colliding chains that are not collisions between
    # the two keep_chains (ie. the last collision)
    bin_generations_merged = []
    # This leaves out the last collision
    for collisions in collision_list[0:-1]:
        for col_tuple in collisions[1]:
            col_gen_idx = collisions[0]
            bin_generations_merged = \
                merge_chain(bin_generations, col_gen_idx, col_tuple)

    print("bin_generations_merged: ")
    print_generations(bin_generations_merged)

    gen_idx_of_last_collision = collision_list[-1][0]

    print("gen_idx_of_last_collision: " + str(gen_idx_of_last_collision))

    # In the "else" branch, we are merging dead ends
    # (without merging, we would have a bin with three neighbours)
    final_chain1 = []
    for gen_idx, gen in enumerate(bin_generations_merged):
        if gen_idx <= gen_idx_of_last_collision:
            final_chain1.append(gen[keep_chains[0]])
        else:
            print("gen: " + str(gen))
            for atom_idx in gen[keep_chains[0]]:
                final_chain1[-1] = np.append(final_chain1[-1], [atom_idx])
                
                
                
    final_chain2 = []
    for gen_idx, gen in enumerate(bin_generations_merged):
        if gen_idx <= gen_idx_of_last_collision:
            final_chain2.append(gen[keep_chains[1]])
        else:
            print("gen: " + str(gen))
            for atom_idx in gen[keep_chains[1]]:
                final_chain2[-1] = np.append(final_chain2[-1], [atom_idx])



    final_chain = glue_chains(final_chain1, final_chain2)
    
    # print(final_chain)
    
    print("\nfinal_chain: ")
    
    line_str = ""
    for bn in final_chain:
        line_str = line_str + str(bn)
        line_str = line_str + "\n"

    print(line_str)



    write_bins(final_chain, atom_positions, INPUT_FILE_NAME, OPEN_JMOL)


    """
    Problem:
        When chains collide, the tips of the chains can end up in
        different generations (because one of the chains snatches up
        the last atoms between the chains)
        In this case, collisions are not detected at the moment.
        Possible solutions:
            - Make sure that tips of chains at collision point are always
             in the same generation?
             - Or maybe check for collisions between generations?
    """
    """
    ToDo:
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
