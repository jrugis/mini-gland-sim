#!/usr/bin/env python
# coding: utf-8

import os
import argparse

import h5py
import numpy as np


def main(results_file, restart_file):
    with h5py.File(results_file, 'r') as results:
        print(f"Loading results from {results_file}")

        # load the x array
        x = np.asarray(results['_duct/x'])
        print("  x shape", x.shape)

        # just the final time step
        xfinal = x[-1, :]
        print("  x final shape", xfinal.shape)

        # attributes for x array
        for key in results['_duct'].attrs.keys():
            print(f"{key}: {results['_duct'].attrs[key]}")
        v_c = results['_duct'].attrs['number of cellular variables']
        v_l = results['_duct'].attrs['number of lumenal variables']
        n_c = results['_duct'].attrs['number of cells']
        n_l = results['_duct'].attrs['number of lumen discs']

    # create restart HDF5 file
    print(f"Writing initial conditions file to {restart_file}")
    with h5py.File(restart_file, 'w') as fh:
        fh.create_dataset("/x", data=xfinal)
        fh['/x'].attrs['number of cellular variables'] = v_c
        fh['/x'].attrs['number of lumenal variables'] = v_l
        fh['/x'].attrs['number of cells'] = n_c
        fh['/x'].attrs['number of lumen discs'] = n_l


def get_args():
    parser = argparse.ArgumentParser(description='Create a restart file containing initial conditions from a results file.')
    parser.add_argument("-i", "--results-file", default="_duct_results.h5", help="The result file to read the initial conditions from")
    parser.add_argument("-o", "--restart-file", default="initial_conditions.h5", help="The restart file to create")
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_args()
    main(args.results_file, args.restart_file)
