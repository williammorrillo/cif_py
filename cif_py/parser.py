"""
Parser for cif files.
"""

import os
import re
import numpy as np


def parse_loop(block, dict):
    """
    Parse a single block of a cif file.
    """
    # retreive all data after loop_
    try:
        loop = block.split('loop_')[1]
        loop = loop.strip()

        # split into lines
        lines = loop.split('\n')

        titles = [x.strip() for x in lines if x.strip().startswith('_')]
        data = [x.strip().split(' ') for x in lines
                if not x.strip().startswith('_')]
        data = np.array(data, dtype=object)

        dict.update({title: data[:, i] for i, title in enumerate(titles)
                     if not (title.startswith('_sym')
                     or title.startswith('_space'))})

        dict.update({title: data for title in titles
                     if title.startswith('_sym')
                     or title.startswith('_space')})
    except:
        pass


def parse_vectors(block, dict):
    """
    looks for lattice constants in the block
    """

    # prefix
    prefix = '_cell_'

    lines = block.split('\n')

    for line in lines:
        if line.startswith(prefix):
            line = line.strip()
            line = line.split(' ')
            line = [x for x in line if x != '']
            line[1] = line[1].split('(')[0]
            dict.update({line[0]: line[1:]})


def parse_formula(block, dict):
    """
    looks for chemical formula in the block
    """

    prefix = '_chemical_'

    lines = block.split('\n')

    for line in lines:
        if line.startswith(prefix):
            line = line.strip()
            line = line.split(' ')
            line = [x for x in line if x != '']
            dict.update({line[0]: line[1:]})


def parse_cif(cif_file):
    """
    Parse a cif file
    """

    with open(cif_file, 'r') as f:
        cif = f.read()

    # separate into blocks by empty lines
    blocks = re.split(r'\n\s*\n', cif)

    # parse blocks
    data = {}

    for block in blocks:
        parse_loop(block, data)
        parse_vectors(block, data)
        parse_formula(block, data)

    return data
