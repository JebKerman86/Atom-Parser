# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:35:15 2017

@author: Benjamin
"""


import numpy as np

from prep_input import prep_data
from file_io import chache_data, load_data
from utilities import print_matrix

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

    dist_matrix = data["dist_mtrx"]
    interaction_matrix = data["interact_mtrx"]
    ordered_index_matrix = data["ordered_idx_mtrx"]
    ordered_dist_matrix = data["ordered_dist_mtrx"]
    ordered_interaction_matrix = data["ordered_interact_mtrx"]

    region_list = data["region_list"]

    device = region_list[0]
    contact1 = region_list[1]  # contacts must start at index "1"
    contact2 = region_list[2]
    #print(device)
    #print(contact1)
    #print(contact2)
    # Each row contains the indices of atoms in a contact
    # that interact with device atoms.
    contact_edge_matrix = []
    """
    for contact in region_list[1:]:
        # print("contact: " + str(contact))
        contact_edge_list = []
        for index_contact in contact:
            for index_device in device:
                # print(index_contact)
                # print(index_device)
                if interaction_matrix[index_contact][index_device]:
                    if not (index_contact in contact_edge_list):
                        contact_edge_list.append(index_contact)
        contact_edge_matrix.append(contact_edge_list)

    # print_matrix(dist_matrix)
    print_matrix(interaction_matrix)
    print(contact_edge_matrix)
    """

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
