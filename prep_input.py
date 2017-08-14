# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:38:06 2017

@author: Benjamin
"""
import numpy as np
from numpy import linalg as LA

from file_io import read_xyz_file, read_transport_file


# -----------------------------------------------------------------------------

# Vectorize this function ???
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

    return bool(dist_matrix[index1, index2] < inter_dist)


# -----------------------------------------------------------------------------


def prep_data(atom_types, atom_positions, region_list, interaction_distances):
    """
    Determine the following (N x N) Matrices ( 0 <= n, m < N, N = numAtoms ):
    <*> dist_mtrx          ------- Distance between atom n and atom m
    <*> interact_mtrx      ------- Entry n,m is "True", if atoms n,m are
                                    within interaction distance
    <*> ordered_idx_mtrx   ------- Each line corresponds to a fixed atom. The
                                    The elements of the line are the indices
                                    of the fixed atoms neighbours, in order
                                    of increasing distance.
    <*> ordered_dist_mtrx  ------- Take ordered_index_matrix, and replace the
                                    indices of the atom with their distances
                                    to the fixed atom.
    <*> ordered_interact_mtrx  --- Take ordered_dist_matrix, and set entry
                                    to "True" if distance is lower than the
                                    interaction distance for that pair of atoms
    """

    num_atoms = len(atom_positions)
    atom_positions = np.array(atom_positions)

    dist_mtrx = np.zeros((num_atoms, num_atoms))
    for idx1 in range(0, num_atoms):
        for idx2 in range(0, num_atoms):
            # print("ixd1 = " + str(idx1) + " / " + "idx2 = " + str(idx2))
            dist = atom_positions[idx2] - atom_positions[idx1]
            # print(dist)
            dist_mtrx[idx1, idx2] = LA.norm(dist)

    interact_mtrx = np.empty((num_atoms, num_atoms), dtype=bool)
    for idx1 in range(num_atoms):
        for idx2 in range(num_atoms):
            is_interacting = is_within_interaction_distance(
                    idx1, idx2, atom_types, interaction_distances, dist_mtrx)
            interact_mtrx[idx1, idx2] = is_interacting

    # ordered_index_matrix:
    # Each row gives a list of atoms by their index, order by increasing
    # distance the column index for each row corresponds to the one atom
    # relative to which all the distances in a row are measured

    ordered_idx_mtrx = np.zeros((num_atoms, num_atoms))
    ordered_dist_mtrx = np.zeros((num_atoms, num_atoms))
    ordered_interact_mtrx = np.zeros((num_atoms, num_atoms))

    ordered_idx_mtrx = np.argsort(dist_mtrx, axis=1)
    ordered_dist_mtrx = np.sort(dist_mtrx, axis=1)

    for idx1 in range(0, num_atoms):
        # print("---------------------------------------------------------")
        # print(ordered_idx_mtrx[idx1,:])
        # print("---------------------------------------------------------")
        for idx2, elem in enumerate(ordered_idx_mtrx[idx1, :]):
            # print("ixd1 = " + str(idx1) + " / " + "idx2 = " + str(idx2))
            is_interacting =  \
                is_within_interaction_distance(
                    idx1, elem, atom_types,
                    interaction_distances, dist_mtrx)
            ordered_interact_mtrx[idx1, idx2] = is_interacting
            # print(is_interacting)
            # print(ordered_interact_mtrx)

    data =  {"dist_mtrx": dist_mtrx,
             "interact_mtrx": interact_mtrx,
             "ordered_idx_mtrx": ordered_idx_mtrx,
             "ordered_dist_mtrx": ordered_dist_mtrx,
             "ordered_interact_mtrx": ordered_interact_mtrx}

    return data
