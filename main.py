# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""

import numpy as np

from prep_input import prep_data
from file_io import chache_data, load_data
from bin_sort import get_contact_bins, get_next_bins

LOAD_CACHE_DATA = False
# Name without file ending:
INPUT_FILE_NAME = "1d_kette"


def main():

    """
    Entry Point for Program
    """

    if LOAD_CACHE_DATA:
        data = load_data(INPUT_FILE_NAME)
    else:
        data = prep_data(INPUT_FILE_NAME)
        chache_data(INPUT_FILE_NAME, data)

    dist_mtrx = data["dist_mtrx"]
    interact_mtrx = data["interact_mtrx"]
    ordered_idx_mtrx = data["ordered_idx_mtrx"]
    ordered_dist_mtrx = data["ordered_dist_mtrx"]
    ordered_interact_mtrx = data["ordered_interact_mtrx"]

    num_atoms = np.size(dist_mtrx, axis=0)
    region_list = data["region_list"]

    # Don't convert to np array in file_io because this complictes caching
    # (since numpy arrays non json-serialzable)
    np_region_list = []
    for region in region_list:
        np_region_list.append(np.array(region))

    # print(np_region_list)

    # List of np-arrays. Each list element corresponds to the contact with the
    # same index. The list element contains the indices of atoms in that
    # contact that are interacting with device atoms
    # contact_bins = []

    device = np_region_list[0]
    contacts = np_region_list[1:]
    
    contact_bins = get_contact_bins(device, contacts, interact_mtrx)


    #print(contact_bins)


    #List of bins. Each element contains the bins corresponding to a contact
    bins = []
    atoms = np.array(list(range(num_atoms)))

    for contact in contact_bins:
        bin_candidates = np.array([])
        for atom_idx in contact:
            bin_add_candidates = atoms[interact_mtrx[atom_idx, :]]
            for cntct in contacts:
                bin_add_candidates = [x for x in bin_add_candidates if x not in cntct]
            bin_candidates = np.r_[bin_candidates, bin_add_candidates]
            #print(bin_add_candidates)
        bin_atoms = np.unique(bin_candidates)
        #print("bin_atoms: " + str(bin_atoms))
        bins.append(bin_atoms)

    print(get_next_bins(contact_bins, contacts, interact_mtrx))
    print(bins)


    """
    ToDo: use NumPy instead of lists
          use version controll
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
