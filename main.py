# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import numpy as np

from prep_input import prep_data
from file_io import chache_data, load_data, write_bins, read_xyz_file, read_transport_file
from bin_sort import get_contact_bins, get_next_bins
from utilities import common_elements

LOAD_CACHE_DATA = False
# Name without file ending:
INPUT_FILE_NAME = "t-kreuzung"
# INPUT_FILE_NAME = "zno2wire"
#INPUT_FILE_NAME = "SiNW"
OPEN_JMOL = False

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
    #print("bin_generations" + str(bin_generations))

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
        #print("prev_bins: " + str(prev_bins))
        bin_generations.append(curr_bins)
        if common_elements(curr_bins)[0]:
            break
        if num_sorted_atoms >= num_atoms:
            break
        curr_generation += 1

    print("num_sorted_atoms: " + str(num_sorted_atoms))



    print("bin_generations: ")
    for gen in bin_generations:
        line_str = ""
        for b in gen:
            line_str = line_str + str(b)
        print(line_str)

    #for now, leave off last generation, because it contains duplicate elements
    #write_bins(bin_generations[:-1], atom_positions, INPUT_FILE_NAME)

    write_bins(bin_generations, atom_positions, INPUT_FILE_NAME, OPEN_JMOL)

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

    # Why no index "9" in contact? Interaction matrix is wrong!


if __name__ == "__main__":
    main()
