#!/usr/bin/env python
# coding: utf-8

import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    results = h5py.File(args.results_file, 'r')
    print(results)

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
    print(cellpos)

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

    # which time step to plot
    plot_time = args.time
    if plot_time is None:
        plot_ind = -1
        plot_time = timevals[-1]
    else:
        plot_ind = np.where(timevals == plot_time)[0][0]

    # full x array for last step
    xfin = x[plot_ind, :]
    print("xfin shape", xfin.shape)

    # extract cellular array
    x_c = np.transpose(xfin[:n_c * v_c].reshape(n_c, v_c))
    print("x_c shape", x_c.shape)

    # extract lumenal array
    x_l = np.transpose(xfin[n_c * v_c:].reshape(n_l, v_l))
    print("x_l shape", x_l.shape)

    # create the plot
    nrow = 3
    ncol = 2
    fig, plots = plt.subplots(nrow, ncol, squeeze=False)
    plt.subplots_adjust(wspace=0.5)
    fig.set_size_inches(ncol * 7.6, nrow * 5.0)
    fig.suptitle("t = %.3f s (%s)" % (plot_time, os.getcwd()))

    plots[0, 0].plot(cellpos, x_c[0, sortedidx], label="V_A")
    plots[0, 0].plot(cellpos, x_c[1, sortedidx], label="V_B")
    plots[0, 0].legend(loc='best')
    plots[0, 0].set_title("Membrane Potentials")
    plots[0, 0].set_xlabel(u"Duct Length (\u03bcm)")
    plots[0, 0].set_ylabel("mV")
    plots[0, 0].ticklabel_format(useOffset=False)

    plots[0, 1].plot(cellpos, x_c[2, sortedidx], label="Cell volume")
    plots[0, 1].set_title("Cell volume")
    plots[0, 1].set_xlabel(u"Duct Length (\u03bcm)")
    plots[0, 1].set_ylabel(u"\u03bcm^3")
    plots[0, 1].ticklabel_format(useOffset=False)

    plots[1, 0].plot(cellpos, x_c[3, sortedidx], label="Na_C")
    plots[1, 0].plot(cellpos, x_c[4, sortedidx], label="K_C")
    plots[1, 0].plot(cellpos, x_c[5, sortedidx], label="Cl_C")
    plots[1, 0].plot(cellpos, x_c[6, sortedidx], label="HCO_C")
    plots[1, 0].legend(loc='best')
    plots[1, 0].set_title("Cellular concentrations")
    plots[1, 0].set_xlabel(u"Duct Length (\u03bcm)")
    plots[1, 0].set_ylabel("mM")
    plots[1, 0].ticklabel_format(useOffset=False)

    plots[1, 1].plot(cellpos, -np.log10(x_c[7, sortedidx] * 1e-3), label="Cellular pH")
    plots[1, 1].set_title("Cellular pH")
    plots[1, 1].set_xlabel(u"Duct Length (\u03bcm)")
    plots[1, 1].ticklabel_format(useOffset=False)

    plots[2, 0].plot(intpos, x_l[0, :], label="Na_A")
    plots[2, 0].plot(intpos, x_l[1, :], label="K_A")
    plots[2, 0].plot(intpos, x_l[2, :], label="Cl_A")
    plots[2, 0].plot(intpos, x_l[3, :], label="HCO_A")
    plots[2, 0].legend(loc='best')
    plots[2, 0].set_title("Lumenal concentrations")
    plots[2, 0].set_xlabel(u"Duct Length (\u03bcm)")
    plots[2, 0].set_ylabel("mM")
    plots[2, 0].ticklabel_format(useOffset=False)

    plots[2, 1].plot(intpos, -np.log10(x_l[4, :] * 1e-3), label="Lumenal pH")
    plots[2, 1].set_title("Lumenal pH")
    plots[2, 1].set_xlabel(u"Duct Length (\u03bcm)")
    plots[2, 1].ticklabel_format(useOffset=False)

    if args.output_file is None:
        args.output_file = f"plot_state_{plot_time:.1f}s.pdf"
    print(f"writing {args.output_file}")
    fig.savefig(args.output_file)


def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot dynamic simulation result at given time point")
    parser.add_argument('-t', '--time', type=int, default=None, help="The time in seconds to plot (default: end)")
    parser.add_argument("-f", "--results-file", default="_duct_results.h5", help="Results file (default: _duct_results.h5)")
    parser.add_argument("-o", "--output-file", default=None, help="Name of file to save figure to (default: plot_state_<time>s.pdf)")
    args = parser.parse_args()

    if not os.path.exists(args.results_file):
        raise ValueError(f"results file does not exist: {args.results_file}")

    return args


if __name__ == "__main__":
    args = get_args()
    main(args)
