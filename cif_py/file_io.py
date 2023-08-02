"""
a module to deal with file io
"""

import os
import h5py


def h5_write(dict, filename):
    """
    writes a h5 file from a dictionary
    """
    filename = os.path.splitext(filename)[0] + '.h5'

    with h5py.File(filename, 'w') as f:
        for key, value in dict.items():
            f.create_dataset(key, data=value)


def h5_read(h5_file_name):
    """
    reads a h5 file and returns a dictionary
    """

    # check if file exists
    if not os.path.isfile(h5_file_name):
        raise FileNotFoundError('File {} not found'.format(h5_file_name))

    with h5py.File(h5_file_name, 'r') as f:
        data = {key: value[:]
                for key, value in f.items()}

    return data
