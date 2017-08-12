# -*- coding: utf-8 -*-
"""
Einlesen von .XYZ Dateien
"""


import json
import os
from pathlib import Path


# -----------------------------------------------------------------------------


def read_xyz_file(input_file_name):
    """
    Reads File "input_file_name".xyz, and returns two lists,
    specifying Atom types and Atom positions, respectively.
    """

    file = open(input_file_name, 'r')
    max_file_lines = 1000

    num_atoms = int(file.readline())
    # file_comment = file.readline()
    atom_positions = []
    atom_types = []

    # skip comment section in file
    line = file.readline()

    iterations = 0
    line = file.readline()
    stripped_line = line.replace(" ", "").replace("\n", "")

    while stripped_line != "" and iterations < max_file_lines:

        line_string = line.split()
        atom_type_string = line_string[0]
        atom_types.append(atom_type_string)
        coord_string = line_string[1:]
        coord_xyz = [float(i) for i in coord_string]
        atom_positions.append(coord_xyz)

        iterations += 1
        line = file.readline()
        stripped_line = line.replace(" ", "").replace("\n", "")

    print("Length of Atom List: " + str(len(atom_positions)))

    if iterations >= max_file_lines:
        print("Input file too long!")

    if num_atoms != len(atom_positions):
        print("falsche Atomanzahl: num_atoms != Atoms.count")

    return (atom_types, atom_positions)


# -----------------------------------------------------------------------------


def read_transport_file(input_file_name):
    """
    Reads File "input_file_name".dat, and returns lists containing the atom
    indices of the device atoms, as well as the atom indices of
    the contact atoms. Also, a dictionary "interaction_distances" is generated,
    which spcifies the maximum interaction distance between each type of atom.
    """

    file = open(input_file_name + ".dat", 'r')
    max_file_lines = 1000

    iterations = 0

    # IMPORTANT: In file, first atom has index is one, but in my program,
    # first atom has index is zero

    region_list = []  # List of regions, starting with device region

    line = file.readline()
    stripped_line = line.replace(" ", "").replace("\n", "")
    entries = line.split()

    iterations = 0
    while stripped_line != ''             \
        and iterations < max_file_lines   \
        and ("Device" in entries[0] or "Contact" in entries[0]):

        region_list.append(tuple(range(int(entries[1]) - 1, int(entries[2]))))
        line = file.readline()
        stripped_line = line.replace(" ", "").replace("\n", "")
        entries = line.split()
        iterations += 1

    interaction_distances = {}

    line = file.readline()
    stripped_line = line.replace(" ", "").replace("\n", "")
    entries = line.split()

    # loop terminates at first empty line, or at end of file
    # (since readline() returns empty string at end of file)
    iterations = 0
    while stripped_line != '' and iterations < max_file_lines:
        key = entries[0] + entries[1]
        interaction_distances[key] = float(entries[2])

        line = file.readline()
        entries = line.split()
        iterations += 1
        stripped_line = line.replace(" ", "").replace("\n", "")

    return (region_list, interaction_distances)


# -----------------------------------------------------------------------------


def chache_data(input_file_name, data):
    """
    Caches matrices from function PrepInput() to lower execution time when
    program called repeatedly with the same input files.
    """
    with open("./cached_data/" + input_file_name + ".json", 'w') as outfile:
        json.dump(data, outfile)


# -----------------------------------------------------------------------------


def load_data(input_file_name):
    """
    Loads cached data for reuse.
    """
    file_name = input_file_name + ".json"
    my_file = Path(file_name)
    if my_file.is_file():
        with open("./cached_data/" + str(file_name), 'r') as file:
            return json.load(file)
    return "No File Found"
