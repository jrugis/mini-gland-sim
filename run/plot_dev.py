#!/usr/bin/env python
# coding: utf-8

import os

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


results = h5py.File('_duct_results.h5', 'r')
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
timevals = np.arange(0.0, x.shape[0]*dt-0.5*dt, dt)
assert timevals.shape[0] == x.shape[0]

# Plotting last step

# full x array for last step
xfin = x[-1, :]
print("xfin shape", xfin.shape)

# extract cellular array
x_c = np.transpose(xfin[:n_c*v_c].reshape(n_c, v_c))
print("x_c shape", x_c.shape)

# extract lumenal array
x_l = np.transpose(xfin[n_c*v_c:].reshape(n_l, v_l))
print("x_l shape", x_l.shape)

# create the plot
fig, plots = plt.subplots(2, 5, squeeze=False)
plt.subplots_adjust(wspace = 0.5)
fig.set_size_inches(5 * 7.6, 2 * 5.0)
fig.suptitle("t = %.3f s (%s)" % (timevals[-1], os.getcwd()))

plots[0,0].plot(cellpos, x_c[0, sortedidx], label="V_A")
plots[0,0].plot(cellpos, x_c[1, sortedidx], label="V_B")
plots[0,0].legend(loc='best')
plots[0,0].set_title("Membrane Potential")
plots[0,0].set_ylabel("mV")

plots[0,1].plot(cellpos, x_c[3, sortedidx], label="Na_C")
plots[0,1].plot(cellpos, x_c[4, sortedidx], label="K_C")
plots[0,1].legend(loc='best')
plots[0,1].set_title("Cellular Concentration")
plots[0,1].set_ylabel("mM")

plots[0,2].plot(cellpos, x_c[5, sortedidx], label="Cl_C")
plots[0,2].plot(cellpos, x_c[6, sortedidx], label="NCO_C")
plots[0,2].legend(loc='best')
plots[0,2].set_title("Cellular Concentration")
plots[0,2].set_ylabel("mM")

plots[0,3].plot(cellpos, -np.log10(x_c[7, sortedidx] * 1e-3))
plots[0,3].set_title("Cellular pH")

plots[0,4].plot(cellpos, x_c[8, sortedidx], label="CO_C")
plots[0,4].legend(loc='best')
plots[0,4].set_title("Cellular Concentration")
plots[0,4].set_ylabel("mM")
plots[0,4].ticklabel_format(useOffset=False)

plots[1,0].plot(cellpos, x_c[2, sortedidx])
plots[1,0].set_title("Cell Volume")
plots[1,0].set_ylabel("um3")

plots[1,1].plot(intpos, x_l[0, :], label="Na_A")
plots[1,1].plot(intpos, x_l[1, :], label="K_A")
plots[1,1].legend(loc='best')
plots[1,1].set_title("Lumenal Concentration")
plots[1,1].set_ylabel("mM")

plots[1,2].plot(intpos, x_l[2, :], label="Cl_A")
plots[1,2].plot(intpos, x_l[3, :], label="HCO_A")
plots[1,2].legend(loc='best')
plots[1,2].set_title("Lumenal Concentration")
plots[1,2].set_ylabel("mM")

plots[1,3].plot(intpos, -np.log10(x_l[4, :] * 1e-3))
plots[1,3].set_title("Lumenal pH")

plots[1,4].plot(intpos, x_l[5, :], label="CO_A")
plots[1,4].legend(loc='best')
plots[1,4].set_title("Lumenal Concentration")
plots[1,4].set_ylabel("mM")
plots[1,4].ticklabel_format(useOffset=False)

print("writing finalstep.pdf", )
fig.savefig("finalstep.pdf")
