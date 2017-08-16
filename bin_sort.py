# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:11:13 2017

@author: Benjamin
"""

import numpy as np

from utilities import remove_all

def get_contact_bins(device, contacts, interact_mtrx):

    contact_starter_atom_list = []

    for contact in contacts:
        contact_edge_list = []
        for index_contact in contact:
            for index_device in device:
                if interact_mtrx[index_contact, index_device]:
                    if not (index_contact in contact_edge_list):
                        contact_edge_list.append(index_contact)
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
            #print(atom_idx)
            #print(interact_mtrx[atom_idx, :])
            bin_add_candidates = atoms[interact_mtrx[atom_idx, :]]
            for prev_bin in prev_bins:
                #print("prev_bin" + str(prev_bin))
                bin_add_candidates = [x for x in bin_add_candidates if x not in prev_bin]
            bin_candidates = np.r_[bin_candidates, bin_add_candidates]
            bin_candidates = [int(x) for x in bin_candidates]
            #print("bin_candidates" + str(bin_candidates))
            #print(bin_add_candidates)
        bin_atoms = np.unique(bin_candidates)
        #print("bin_atoms: " + str(bin_atoms))
        next_bins.append(bin_atoms)
        
    return(next_bins)


def remove_common_elems(curr_gen_cntct_chains, common_elems):
    """
    Takes the bins at tips of each contact chain and ensures that there are no
    duplicate atoms in these bins. The atom in the contact chain with
    the smallest index is kept.
    """
    elems = []
    for cntct_idx, contact_list in enumerate(common_elems):
        elems_to_delete = []
        for idx, atom_idx in enumerate(contact_list):
            if not atom_idx in elems:
                elems.append(atom_idx)
            else:
                elems_to_delete.append(atom_idx)
        for atom_idx in elems_to_delete:
            print("elems_to_delete: " + str(elems_to_delete))
            curr_gen_cntct_chains[cntct_idx] =  \
                remove_all(curr_gen_cntct_chains[cntct_idx], atom_idx)


def contiguity_check(bin_to_check, starter_atom_idx, interact_mtrx):
    # num_atoms = np.size(interact_mtrx, axis=0)
    reached_atoms = []
    # print("num_atoms = " + str(num_atoms))

    if starter_atom_idx not in bin_to_check:
        print("Bad input in contiguity_check: " +
              "The provided ""starter_atom_idx"" was not in bin_to_check.")
        return []

    atoms_to_add_to_reached_atoms = [starter_atom_idx]
    while not atoms_to_add_to_reached_atoms == []:
        reached_atoms = reached_atoms + atoms_to_add_to_reached_atoms
        atoms_to_add_to_reached_atoms = []
        for atom_idx1 in reached_atoms:
            for atom_idx2 in bin_to_check:
                # print("atom_idx2: " + str(atom_idx2))
                if interact_mtrx[atom_idx1, atom_idx2]:
                        if atom_idx2 not in reached_atoms:
                            if atom_idx2 not in atoms_to_add_to_reached_atoms:
                                # print("added atom_idx2: " + str(atom_idx2))
                                atoms_to_add_to_reached_atoms.append(atom_idx2)

    return np.array(reached_atoms)


def create_subdomains(bin_to_divide, interact_mtrx):
    subdomains = []
    # num_atoms = np.size(interact_mtrx, axis=0)

    unreached_atoms = np.copy(bin_to_divide).tolist()
    iterations = 0
    while (not unreached_atoms == []) and (iterations < 1000):
        iterations += 1
        starter_atom_idx = unreached_atoms[0]
        subdomains.append(contiguity_check(bin_to_divide, starter_atom_idx, interact_mtrx))
    
        unreached_atoms = np.copy(bin_to_divide).tolist()
        for sd in subdomains:
            for atom_idx in sd:
                unreached_atoms.remove(atom_idx)

    return list.copy(subdomains)

