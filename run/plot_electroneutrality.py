#!/usr/bin/env python
# coding: utf-8

import os

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


results = h5py.File('d1s1_results.h5', 'r')
print(results)

# attributes for x array
for key in results['d1s1'].attrs.keys():
    print(f"{key}: {results['d1s1'].attrs[key]}")
n_c = results['d1s1'].attrs['number of cells']

# electroneutrality
electro = np.asarray(results['d1s1/electroneutrality'])
print(electro.shape)
assert n_c == electro.shape[1]

# time values
dt = results.attrs["output time interval"]
print("output time interval: %f s" % dt)
timevals = np.arange(0.0, electro.shape[0]*dt-0.5*dt, dt)
assert timevals.shape[0] == electro.shape[0]

cell = np.random.choice(n_c)
print("plotting electroneutrality for cell %d" % (cell+1,))

plt.plot(timevals, electro[:, cell])
plt.title("Electroneutrality check for cell {0}".format(cell+1))
plt.xlabel("Time (s)")
plt.ylabel("Electroneutrality (mol)")

print("writing electroneutrality.pdf", )
plt.savefig("electroneutrality.pdf")
