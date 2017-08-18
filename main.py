# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import numpy as np

from prep_input import prep_data
from file_io import chache_data, load_data, write_bins, read_xyz_file, read_transport_file
from bin_sort import get_contact_bins, get_next_bins, remove_common_elems, \
                     create_subdomains, find_all_collisions, merge_chain, glue_chains
from utilities import find_duplicates, remove_all, print_generations

LOAD_CACHE_DATA = False
# Name without file ending:
INPUT_FILE_NAME = "t-kreuzung_dick"
#INPUT_FILE_NAME = "t-kreuzung_sackgasse"
# INPUT_FILE_NAME = "zno2wire"
#INPUT_FILE_NAME = "SiNW"
OPEN_JMOL = True

MAX_GENERATIONS = 20

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
    bin_generations = []
    bin_generations.append(contact_bins)
    # print("bin_generations" + str(bin_generations))

    curr_generation = 1

    # contact atoms count as sorted from the beginning
    num_sorted_atoms = 0
    for c in contacts:
        num_sorted_atoms += len(c)

    while curr_generation < MAX_GENERATIONS:
        curr_bins = get_next_bins(bin_generations[-1], prev_bins, interact_mtrx)
        print("curr_generation = " + str(curr_generation))
        for b in curr_bins:
            num_sorted_atoms += len(b)
            print("bin: " + str(b))
        prev_bins = prev_bins + curr_bins
        # print("prev_bins: " + str(prev_bins))
        bin_generations.append(curr_bins)
        (duplicates, common_elems) = find_duplicates(curr_bins)
        if duplicates > 0:
            print("common_elems: " + str(common_elems))
            remove_common_elems(bin_generations[-1], common_elems)
            num_sorted_atoms -= duplicates

        if num_sorted_atoms >= num_atoms:
            break
        curr_generation += 1

    print("num_sorted_atoms: " + str(num_sorted_atoms))

    contiguous_bin_generations = []
    for gen in bin_generations:
        contiguous_gen = []
        for bn in gen:
            divided_bin = create_subdomains(bn, interact_mtrx)
            contiguous_gen.append(divided_bin)
        contiguous_bin_generations.append(contiguous_gen)



    print("contiguous_bin_generations: ")
    print_generations(contiguous_bin_generations)


    collision_list = find_all_collisions(contiguous_bin_generations, interact_mtrx)
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
    
    
    
    
    for collisions in collision_list[0:-1]:
        for col_tuple in collisions[1]:
            col_gen_idx = collisions[0]
            contiguous_bin_generations = \
                merge_chain(contiguous_bin_generations, col_gen_idx, col_tuple)

    print("contiguous_bin_generations: ")
    print_generations(contiguous_bin_generations)
    
    final_chain1 = []
    for gen in contiguous_bin_generations:
        final_chain1.append(gen[keep_chains[0]])
        
    final_chain2 = []
    for gen in contiguous_bin_generations:
        final_chain2.append(gen[keep_chains[1]])
    

    
    final_chain = glue_chains(final_chain1, final_chain2)
    
    print("\nfinal_chain: ")
    
    for bn in final_chain:
        line_str = ""
        for idx, sd in enumerate(bn):
            line_str = line_str + str(sd)
            if idx+1 < len(bn):
                line_str = line_str + " -- "
        print(line_str)
        
        
    
    

    write_bins(bin_generations, atom_positions, INPUT_FILE_NAME, OPEN_JMOL)
    
    
    
    # rewrite code to work with subchains instead of subdomains:
    # When merging bins, merge the two subchains that are next to eachother
    # Right now, when merging bins, subdomains are kept separate, even if
    # they are right next to eachother
    
    # BLind alley (Sackgasse) detection:
    # A subchain that doesn't end in a collision is blind alley.
    # These need to be trimmed off before any chain merging is attempted.

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
