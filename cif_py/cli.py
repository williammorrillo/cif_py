"""
Command line interface for cif_py
"""
import argparse
import os

from h5py._hl import files

from . import parser
from . import output
from . import file_io


def main(args):
    """
    Main function for parsing cif files.
    """

    cif_file_name = args.cif_file_name
    h5_file_name = args.h5_file_name
    vasp = args.vasp

    if cif_file_name is None:
        # check if file exists
        if not os.path.isfile(cif_file_name):
            raise FileNotFoundError('File {} not found'.format(cif_file_name))

        data = parser.parse_cif(cif_file_name)
        file_io.h5_write(data, cif_file_name)

    elif h5_file_name is None:
        # check if file exists
        if not os.path.isfile(h5_file_name):
            raise FileNotFoundError('File {} not found'.format(h5_file_name))

        data = file_io.h5_read(h5_file_name)

    else:
        raise NotImplementedError(
                'Please select the input format (cif or h5)')

    # write output
    output_file = output.OutputWriter(data)

    if vasp:
        output_file.write_vasp()
    else:
        raise NotImplementedError(
                'Please select the output format')


def read_args():
    """
    Read command line arguments.
    """
    parser = argparse.ArgumentParser(
            prog='cif_py',
            description=dedent('''
                               Parser for cif files

                               Input: cif or h5 file
                               -- cif_file_name: cif file to parse
                               -- h5_file_name: h5 file to parse

                               output:
                               -- vasp: POSCAR
                               '''),
            epilog='Author: William Morrillo',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
            '--cif_file_name',
            type=str,
            help='cif file to parse',
            default=None)

    parser.add_argument(
            '--vasp',
            action='store_true',
            help='write vasp output')

    parser.add_argument(
            '--h5_file_name',
            type=str,
            help='h5 file to parse',
            default=None)

    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    args = read_args()
