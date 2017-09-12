# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:11:13 2017

@author: Benjamin
"""
import sys
import numpy as np
from copy import deepcopy

def get_contact_bins(device, contacts, interact_mtrx):
    """
    Returns the nested list "contact_starter_atom_list".
    For each contact, this nested list contains another list of all the
    atoms in said contact that are interacting with device atoms.
    """
    contact_starter_atom_list = []

    for contact in contacts:
        contact_edge_list = []
        for contact_idx in contact:
            for device_idx in device:
                if interact_mtrx[contact_idx, device_idx]:
                    if not contact_idx in contact_edge_list:
                        contact_edge_list.append(contact_idx)
        contact_starter_atom_list.append(np.array(contact_edge_list))

    return contact_starter_atom_list


def get_next_bins(curr_bins, prev_bins, interact_mtrx):
    """
    List of bins. Each element contains the bins corresponding to a contact.
    """
    next_bins = []
    num_atoms = np.size(interact_mtrx, axis=0)
    atoms = np.array(list(range(num_atoms)))

    for curr_bin in curr_bins:
        bin_candidates = np.array([])
        for atom_idx in curr_bin:
            # print(atom_idx)
            # print(interact_mtrx[atom_idx, :])
            bin_add_candidates = atoms[interact_mtrx[atom_idx, :]]
            for prev_bin in prev_bins:
                # print("prev_bin" + str(prev_bin))
                bin_add_candidates = [x for x in bin_add_candidates if x not in prev_bin]
                bin_add_candidates = np.array(bin_add_candidates)
            bin_candidates = np.append(bin_candidates, bin_add_candidates)
            bin_candidates = [int(x) for x in bin_candidates]
            # print("bin_add_candidates" + str(bin_add_candidates))
            # print("bin_candidates" + str(bin_candidates))
        bin_atoms = np.unique(bin_candidates)
        # print("bin_atoms: " + str(bin_atoms))
        next_bins.append(bin_atoms)

    return next_bins


def bins_are_neighbours(bin1, bin2, interact_mtrx):
    """
    Check if bin1 and bin2 contain interacting atoms, or have atoms in common.
    """
    max_atom_idx = len(interact_mtrx[:, 0])-1
    for atom_idx1 in bin1:
        if atom_idx1 > max_atom_idx:
            sys.exit("FATAL ERROR: atom_idx out of range!")
            continue
        for atom_idx2 in bin2:
            if atom_idx2 > max_atom_idx:
                sys.exit("FATAL ERROR: atom_idx out of range!")
                continue
            if interact_mtrx[atom_idx1, atom_idx2]:
                return True
    return False


def get_chain_length(chains, chain_idx, start_gen_idx):
    """
    Get length of the chain with index "chain_idx", starting from (and including)
    generation "start_gen_idx" to end of chain, or until first
    empty bin (while excluding empty bin).
    """
    
    length = 0
    for gen_idx, gen in enumerate(chains[start_gen_idx:]):
        bn = gen[chain_idx]
        if len(bn) == 0:
            break
        length += 1
        # print("\nbn: " + str([x+1 for x in bn]))
    
    return length


def get_dead_ends(chains, final_chain_idxs, gen_idx_of_last_collision):
    vrb_local = False
    dead_ends = []

    dead_end_start_idx = gen_idx_of_last_collision+1
    for chain_idx in final_chain_idxs:
        if vrb_local: print("chain_idx = " + str(chain_idx))
        if vrb_local: print("dead_end_start_idx = " + str(dead_end_start_idx))
        length = get_chain_length(chains, chain_idx, dead_end_start_idx)
        if vrb_local: print("length = " + str(length))
        if length == 0:
            dead_ends.append([])
        else:
            dead_end = []
            for gen in chains[dead_end_start_idx:]:
                if vrb_local: print("gen: " + str([x+1 for x in gen]))
                dead_end.append(gen[chain_idx])
            dead_ends.append(deepcopy(dead_end))

    return dead_ends
