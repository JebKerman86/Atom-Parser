# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:11:13 2017

@author: Benjamin
"""

import numpy as np
from copy import deepcopy

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


def remove_common_elems(last_gen_cntct_chains, curr_gen_cntct_chains, common_elems):
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
                print("add atom_idx: " + str(atom_idx))
                # last_gen_cntct_chains[cntct_idx] = \
                #     np.r_[ last_gen_cntct_chains[cntct_idx], np.array([atom_idx]) ]
                # last_gen_cntct_chains[cntct_idx] = \
                #     [int(x) for x in last_gen_cntct_chains[cntct_idx]]
                # print(last_gen_cntct_chains[cntct_idx])
                # elems_to_delete.append(atom_idx)
            else:
                elems_to_delete.append(atom_idx)
        for atom_idx in elems_to_delete:
            # print("elems_to_delete: " + str(elems_to_delete))
            curr_gen_cntct_chains[cntct_idx] =  \
                remove_all(curr_gen_cntct_chains[cntct_idx], atom_idx)



def bins_are_neighbours(bin1, bin2, interact_mtrx):

    for atom_idx1 in bin1:
        for atom_idx2 in bin2:
            if atom_idx2 > atom_idx1:
                if interact_mtrx[atom_idx1, atom_idx2]:
                    return True
    return False


def check_generation_for_collisions(generation, interact_mtrx):
    
    # Return list of collision tuples, since there can be multiple collisions
    # in a generation
    collisions = []
    for bin_idx1, bin1 in enumerate(generation):
        for bin_idx2, bin2 in enumerate(generation):
            if bin_idx2 > bin_idx1:
                col = bins_are_neighbours(bin1, bin2,interact_mtrx)
                if col:
                    collisions.append((bin_idx1, bin_idx2))
                    
    if collisions == []:
        return (False, [(-1,-1)])
    else:
        return (True, collisions)


def find_all_collisions(bin_generations, interact_mtrx):

    # col_list is list of tuples:
    # [ (gen_idx1, [(chain_idx1,chain_idx2), (chain_idx1,chain_idx2)]),
    #   (gen_idx2, [(chain_idx1, chain_idx2)]) ]
    col_list = []

    for gen_idx, gen in enumerate(bin_generations):
        # PROBLEM: WHAT HAPPENS IF THERE ARE COLLISIONS BETWEEN BINS IN DIFFERENT
        # GENERATIONS? 
        col = check_generation_for_collisions(gen, interact_mtrx)
        if col[0]:
            col_list.append((gen_idx, col[1]))

    return col_list


def merge_chain(merged_bin_generations, bin_generations, col_gen_idx, col_tuple):

    # Handle cases when more than two chains collide at same time at same place
    # col_tuple[0] always contains the smaller chain index
    # print("col_tuple")
    # print(col_tuple)
    # print("merged_bin_generations[0][2]: ")
    # print(merged_bin_generations[0][2] + merged_bin_generations[0][0])
    for gen_idx in range(col_gen_idx+1):
        # print("gen_idx " + str(gen_idx))
        merged_bin = np.concatenate(
                (deepcopy(merged_bin_generations[gen_idx][col_tuple[0]]),
                 deepcopy(merged_bin_generations[gen_idx][col_tuple[1]])),
                 axis=0)
        merged_bin = [int(x) for x in merged_bin]
        print("merged_bin: ")
        print(merged_bin)
        merged_bin_generations[gen_idx][col_tuple[0]] = merged_bin
        merged_bin_generations[gen_idx][col_tuple[1]] = []
        
    # print("merged_bin_generations: ")
    # print(merged_bin_generations)



def glue_chains(chain1, chain2):
    """
    Glue dangling ends of two chains together:
    [0 --> chain1  -->  -1]  +  [-1  -->  chain2  -->  0]
    """
    
    new_chain = []
    
    for bn in chain1:
        new_chain.append(deepcopy(bn))
        
    for ii in range(1, len(chain2)+1):
        idx = len(chain2) - ii
        new_chain.append(deepcopy(chain2[idx]))
        
    return new_chain



"""
def find_contacts_to_keep(bin_generations):
    
    Find the two contacts that should be left over after merging.
    These are the two contacts that have the most bins (principal layers)
    between them.
    
    num_ctcts = len(bin_generations[0])
    for cntct_idx1 in num_ctcts:
        for cntc_idx2 in num_ctcts:
            if cntc_idx2 > cntct_idx1:
                
"""


