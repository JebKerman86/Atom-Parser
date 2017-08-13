# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:38:06 2017

@author: Benjamin
"""
import numpy as np
from numpy import linalg as LA

from file_io import read_xyz_file, read_transport_file
from utilities import vec3_dist, order_index_list, print_matrix, np_order


# -----------------------------------------------------------------------------

#Vectorize this function ???
def is_within_interaction_distance(index1, index2, atom_types,
                                   interaction_distances, dist_matrix):
    """
    Check whether atoms specified by indices "index1" and "index2" are close
    enough to interact with eachother. If yes, return true.
    """

    interaction_key = str(atom_types[index1]) + str(atom_types[index2])
    inter_dist = interaction_distances.get(interaction_key, -1.0)
    if inter_dist == -1:
        interaction_key = str(atom_types[index2]) + str(atom_types[index1])
        inter_dist = interaction_distances.get(interaction_key, -1.0)

    return bool(dist_matrix[index1,index2] < inter_dist)


# -----------------------------------------------------------------------------


def prep_data(input_file_name):
    """
    Determine the following (N x N) Matrices ( 0 <= n, m < N, N = numAtoms ):
    <*> dist_matrix          ------- Distance between atom n and atom m
    <*> interaction_matrix   ------- Entry n,m is "True", if atoms n,m are
                                     within interaction distance
    <*> ordered_index_matrix ------- Each line corresponds to a fixed atom. The
                                     The elements of the line are the indices
                                     of the fixed atoms neighbours, in order
                                     of increasing distance.
    <*> ordered_dist_matrix  ------- Take ordered_index_matrix, and replace the
                                     indices of the atom with their distances
                                     to the fixed atom.
    <*> ordered_interaction_matrix - Take ordered_dist_matrix, and set entry
                                     to "True" if distance is lower than the
                                     interaction distance for that pair of atoms
    """

    print("Input File: " + str(input_file_name) + ".xyz")

    (atom_types, atom_positions) = read_xyz_file(str(input_file_name))

    (region_list, interaction_distances) = read_transport_file()
    
    num_atoms = len(atom_positions)
    np_atom_positions = np.array(atom_positions)

    #dist_matrix = []
    np_dist_matrix = np.zeros((num_atoms,num_atoms))
    

    for idx1 in range(0, num_atoms):
        for idx2 in range(0, num_atoms):
            #print("ixd1 = " + str(idx1) + " / " + "idx2 = " + str(idx2))
            dist = np_atom_positions[idx2] - np_atom_positions[idx1]
            #print(dist)
            np_dist_matrix[idx1,idx2] = LA.norm(dist)

    np_interaction_matrix = np.zeros((num_atoms,num_atoms))

    for idx1 in range(num_atoms):
        for idx2 in range(num_atoms):
            is_interacting = is_within_interaction_distance(
                    idx1, idx2, atom_types, interaction_distances, np_dist_matrix)
            np_interaction_matrix[idx1,idx2] = is_interacting

    # ordered_index_matrix:
    # Each row gives a list of atoms by their index, order by increasing
    # distance the column index for each row corresponds to the one atom
    # relative to which all the distances in a row are measured
    
    np_ordered_index_matrix = np.zeros((num_atoms,num_atoms))
    np_ordered_dist_matrix = np.zeros((num_atoms,num_atoms))
    np_ordered_interaction_matrix = np.zeros((num_atoms,num_atoms))

    np_ordered_index_matrix = np.argsort(np_dist_matrix, axis = 1)

    np_ordered_dist_matrix = np.sort(np_dist_matrix, axis = 1)
    
    for idx1 in range(0, num_atoms):
        #print("---------------------------------------------------------")
        #print(np_ordered_index_matrix[idx1,:])
        #print("---------------------------------------------------------")
        for idx2, elem in enumerate(np_ordered_index_matrix[idx1,:]):
            #print("ixd1 = " + str(idx1) + " / " + "idx2 = " + str(idx2))
            is_interacting =  \
                is_within_interaction_distance(
                    idx1, elem, atom_types,
                    interaction_distances, np_dist_matrix)
            np_ordered_interaction_matrix[idx1, idx2] = is_interacting
            #print(is_interacting)
            #print(np_ordered_interaction_matrix)



    #print(np.array(ordered_index_matrix) - np_ordered_index_matrix)
    #print(np.array(ordered_dist_matrix) - np_ordered_dist_matrix)

    print("np_ordered_interaction_matrix")
    print(np_ordered_interaction_matrix)


    data = {"region_list"                :  region_list,
            "dist_matrix"                :  np_dist_matrix,
            "interaction_matrix"         :  np_interaction_matrix,
            "ordered_index_matrix"       :  np_ordered_index_matrix,
            "ordered_dist_matrix"        :  np_ordered_dist_matrix,
            "ordered_interaction_matrix" :  np_ordered_interaction_matrix}

    return data
