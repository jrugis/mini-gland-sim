#!/usr/bin/env python
# coding: utf-8

import os
import math

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


VARS = [
    "V_A",
    "V_B",
    "Cell volume",
    "Na_C",
    "K_C",
    "Cl_C",
    "NCO_C",
    "Cellular pH",
    "CO_C"
]

UNITS = [
    "(mV)",
    "(mV)",
    "(um3)",
    "(mM)",
    "(mM)",
    "(mM)",
    "(mM)",
    "",
    "mM"
]

SCALE = {
    "Cellular pH": lambda x: -np.log10(x * 1e-3),
}


def main(args):
    # index of the chosen variable
    varidx = VARS.index(args.variable)

    # load results
    results = h5py.File(args.results_file, 'r')
    print("file", results)

    x = np.asarray(results['_duct/x'])
    print("x shape", x.shape)

    # cell position
    cellpos = np.asarray(results['_duct/CellPos'])
    print("CellPos shape", cellpos.shape)

    # sorting cell positions
#    sortedidx = np.argsort(cellpos)
#    cellpos = np.sort(cellpos)
#    sortedidx = np.flip(sortedidx)
#    cellpos = np.flip(cellpos)
#    print(cellpos)

    # lumen segment positions
    intpos = np.asarray(results['_duct/IntPos'])
    print("IntPos shape", intpos.shape)

    # attributes for x array
    for key in results['_duct'].attrs.keys():
        print(f"{key}: {results['_duct'].attrs[key]}")
    v_c = results['_duct'].attrs['number of cellular variables']
    v_l = results['_duct'].attrs['number of lumenal variables']
    n_c = results['_duct'].attrs['number of cells']
    n_l = results['_duct'].attrs['number of lumen discs']

    # time values
    dt = results.attrs["output time interval"]
    print("output time interval: %f s" % dt)
    timevals = np.arange(0.0, x.shape[0] * dt - 0.5 * dt, dt)
    assert timevals.shape[0] == x.shape[0]

    # limiting time
    mintime = args.min_time if args.min_time is not None else timevals[0] - 1
    maxtime = args.max_time if args.max_time is not None else timevals[-1] + 1
    print("min time", mintime)
    print("max time", maxtime)
    indx = np.logical_and(timevals >= mintime, timevals <= maxtime)
    tbeg = timevals[indx]
    xbeg = x[indx, :]
    print("xbeg shape", xbeg.shape)
    xcbeg = xbeg[:, :n_c * v_c].reshape(len(tbeg), n_c, v_c)
    print('xcbeg shape', xcbeg.shape)

    # now plot

    # cell plot, n_c is number of cells
    ncol = math.ceil(math.sqrt(n_c))
    nrow = math.ceil(n_c / ncol)

    # create the plot
    fig, plots = plt.subplots(nrow, ncol, squeeze=False)
    plt.subplots_adjust(wspace=0.5)
    fig.set_size_inches(ncol * 7.6, nrow * 5.0)
#    fig.suptitle("t = %.3f s (%s)" % (timevals[-1], os.getcwd()))

    # loop over cells
    for cell0 in range(n_c):
        # subplot coord
        row = cell0 // ncol
        col = cell0 % ncol

        # yvals, scale if necessary
        yvals = xcbeg[:, cell0, varidx]
        if args.variable in SCALE:
            yvals = SCALE[args.variable](yvals)

        # plot
        plots[row, col].plot(tbeg, yvals)
        plots[row, col].set_title(f"Cell {cell0+1} ({cellpos[cell0]:.1f} um)")
        plots[row, col].set_xlabel("time (s)")
        plots[row, col].set_ylabel(f"{args.variable} {UNITS[varidx]}")
        plots[row, col].ticklabel_format(useOffset=False)

    # save to file
    print(f"saving result to: {args.output_file}")
    fig.savefig(args.output_file)


def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot a variable over time for all cells")
    parser.add_argument("-v", "--variable", default="Cell volume", help=f"Variable to plot ({', '.join(VARS)}) (default: Cell volume)")
    parser.add_argument("-m", "--min-time", type=float, default=None, help="Minimum time to plot (default: plot all times)")
    parser.add_argument("-M", "--max-time", type=float, default=None, help="Maximum time to plot (default: plot all times)")
    parser.add_argument("-f", "--results-file", default="_duct_results.h5", help="Results file (default: _duct_results.h5)")
    parser.add_argument("-o", "--output-file", default=None, help="Name of file to save figure to (default: plot_cell_time_<variable>.pdf)")
    args = parser.parse_args()

    if args.variable not in VARS:
        raise ValueError(f"variable must be one of: {', '.join(VARS)}")

    if not os.path.exists(args.results_file):
        raise ValueError(f"results file does not exist: {args.results_file}")

    if args.output_file is None:
        args.output_file = f"plot_cell_time_{args.variable.replace(' ', '-')}.pdf"

    return args


if __name__ == "__main__":
    args = get_args()
    main(args)
