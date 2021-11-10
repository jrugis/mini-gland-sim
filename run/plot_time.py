#!/usr/bin/env python
# coding: utf-8

import os

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    # load results
    results = h5py.File(args.results_file, 'r')
    print("file", results)

    x = np.asarray(results['_duct/x'])
    print("x shape", x.shape)

    # cell position
    cellpos = np.asarray(results['_duct/CellPos'])
    print("CellPos shape", cellpos.shape)

    # sorting cell positions
    sortedidx = np.argsort(cellpos)
    cellpos = np.sort(cellpos)
    sortedidx = np.flip(sortedidx)
    cellpos = np.flip(cellpos)

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

    if args.cell > n_c or args.cell < 1:
        raise ValueError(f'Cell {args.cell} is out of range!')
    cell_no = args.cell - 1

    # which discs this cell interfaces with
    disc_mask = np.asarray(results['_duct/CellDiscMask'])
    cell_disc_mask = disc_mask[cell_no][:]
    disc_index = np.where(cell_disc_mask == 1)[0][0]  # just take the first disc this cell interacts with
    print("disc index", disc_index)

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
    xlbeg = xbeg[:, n_c * v_c:].reshape(len(tbeg), n_l, v_l)
    print('xlbeg shape', xlbeg.shape)

    # just the cell and disc we want to plot
    yy_c = xcbeg[:, cell_no, :]
    yy_l = xlbeg[:, disc_index, :]

    # now plot
    nrow = 3
    ncol = 2
    fig, plots = plt.subplots(nrow, ncol, squeeze=False)
    plt.subplots_adjust(wspace=0.5)
    fig.set_size_inches(ncol * 7.6, nrow * 5.0)
    fig.suptitle(f"Cell {cell_no} ({os.getcwd()})")

    plots[0, 0].plot(tbeg, yy_c[:, 0], label="V_A")
    plots[0, 0].plot(tbeg, yy_c[:, 1], label="V_B")
    plots[0, 0].legend(loc='best')
    plots[0, 0].set_title("Membrane Potentials")
    plots[0, 0].set_xlabel("time (s)")
    plots[0, 0].set_ylabel("mV")
    plots[0, 0].ticklabel_format(useOffset=False)

    plots[0, 1].plot(tbeg, yy_c[:, 2], label="Cell volume")
    plots[0, 1].set_title("Cell volume")
    plots[0, 1].set_xlabel("time (s)")
    plots[0, 1].set_ylabel("\\mu m^3")
    plots[0, 1].ticklabel_format(useOffset=False)

    plots[1, 0].plot(tbeg, yy_c[:, 3], label="Na_C")
    plots[1, 0].plot(tbeg, yy_c[:, 4], label="K_C")
    plots[1, 0].plot(tbeg, yy_c[:, 5], label="Cl_C")
    plots[1, 0].plot(tbeg, yy_c[:, 6], label="HCO_C")
    plots[1, 0].legend(loc='best')
    plots[1, 0].set_title("Cellular concentrations")
    plots[1, 0].set_xlabel("time (s)")
    plots[1, 0].set_ylabel("mM")
    plots[1, 0].ticklabel_format(useOffset=False)

    plots[1, 1].plot(tbeg, -np.log10(yy_c[:, 7] * 1e-3), label="Cellular pH")
    plots[1, 1].set_title("Cellular pH")
    plots[1, 1].set_xlabel("time (s)")
    plots[1, 1].ticklabel_format(useOffset=False)

    plots[2, 0].plot(tbeg, yy_l[:, 0], label="Na_A")
    plots[2, 0].plot(tbeg, yy_l[:, 1], label="K_A")
    plots[2, 0].plot(tbeg, yy_l[:, 2], label="Cl_A")
    plots[2, 0].plot(tbeg, yy_l[:, 3], label="HCO_A")
    plots[2, 0].legend(loc='best')
    plots[2, 0].set_title("Lumenal concentrations")
    plots[2, 0].set_xlabel("time (s)")
    plots[2, 0].set_ylabel("mM")
    plots[2, 0].ticklabel_format(useOffset=False)

    plots[2, 1].plot(tbeg, -np.log10(yy_l[:, 4] * 1e-3), label="Lumenal pH")
    plots[2, 1].set_title("Lumenal pH")
    plots[2, 1].set_xlabel("time (s)")
    plots[2, 1].ticklabel_format(useOffset=False)

    # save to file
    print(f"saving result to: {args.output_file}")
    fig.savefig(args.output_file)


def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot dynamic simulation result")
    parser.add_argument('-c', '--cell', type=int, default=1, help="The cell to plot (1-based index, default: 1)")
    parser.add_argument("-m", "--min-time", type=float, default=None, help="Minimum time to plot (default: plot all times)")
    parser.add_argument("-M", "--max-time", type=float, default=None, help="Maximum time to plot (default: plot all times)")
    parser.add_argument("-f", "--results-file", default="_duct_results.h5", help="Results file (default: _duct_results.h5)")
    parser.add_argument("-o", "--output-file", default=None, help="Name of file to save figure to (default: plot_time_cell<cell_number>.pdf)")
    args = parser.parse_args()

    if not os.path.exists(args.results_file):
        raise ValueError(f"results file does not exist: {args.results_file}")

    if args.output_file is None:
        args.output_file = f"plot_time_cell{args.cell}.pdf"

    return args


if __name__ == "__main__":
    args = get_args()
    main(args)
