# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:38:06 2017

@author: Benjamin
"""


from file_io import read_xyz_file, read_transport_file
from utilities import vec3_dist, order_index_list


# -----------------------------------------------------------------------------


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

    # print(dist_matrix[index1][index2])

    return bool(dist_matrix[index1][index2] < inter_dist)


# -----------------------------------------------------------------------------


# Nmpy arrays verwenden!!
# variablen i, j länger machen (zB ii, jj)

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

    (atom_types, atom_positions) = read_xyz_file("input_files/" + str(input_file_name) + ".xyz")

    (region_list, interaction_distances) = read_transport_file("input_files/" + "transport")

    num_atoms = len(atom_positions)

    dist_matrix = []

    for i in range(0, num_atoms):
        dist_matrix.append([])
        for j in range(0, num_atoms):
            dist_matrix[i].append(vec3_dist(atom_positions[i],
                                            atom_positions[j]))

    interaction_matrix = []

    for i in range(num_atoms):
        interaction_list = []
        for j in range(num_atoms):
            interaction_list.append(
                is_within_interaction_distance(
                    i, j, atom_types, interaction_distances, dist_matrix))
        interaction_matrix.append(interaction_list)

    ordered_index_matrix = []
    ordered_dist_matrix = []
    ordered_interaction_matrix = []

    # ordered_index_matrix:
    # Each row gives a list of atoms by their index, order by increasing
    # distance the column index for each row corresponds to the one atom
    # relative to which all the distances in a row are measured

    print("---------------------------------------------------------------")
    for fixed_atom_index in range(num_atoms):
        dist_list = dist_matrix[fixed_atom_index][:]
        ordered_index_list = order_index_list(dist_list)

        ordered_dist_list = []
        for indices in ordered_index_list:
            ordered_dist_list.append(dist_list[indices])

        ordered_interaction_list = []
        for elem in ordered_index_list:
            ordered_interaction_list.append(is_within_interaction_distance(
                fixed_atom_index, elem, atom_types,
                interaction_distances, dist_matrix))

        ordered_index_matrix.append(ordered_index_list)
        ordered_dist_matrix.append(ordered_dist_list)
        ordered_interaction_matrix.append(ordered_interaction_list)

    """
    print(atom_types)

    for i in range(num_atoms):
        print(ordered_index_matrix[i])
        print(ordered_dist_matrix[i])
    """

    data = {"region_list"                :  region_list,
            "dist_matrix"                :  dist_matrix,
            "interaction_matrix"         :  interaction_matrix,
            "ordered_index_matrix"       :  ordered_index_matrix,
            "ordered_dist_matrix"        :  ordered_dist_matrix,
            "ordered_interaction_matrix" :  ordered_interaction_matrix}

    return data
