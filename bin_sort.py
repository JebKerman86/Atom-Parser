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

    #List of bins. Each element contains the bins corresponding to a contact
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



