# This file is distributed under MIT license as part of the project flocksims.
# See LICENSE file for details.
# (c) Simon Grosse-Holz, 2019

# Quick script to sweep some parameters for the flocking simulations
import os, sys
import numpy as np

import mkl
mkl.set_num_threads(1)
from multiprocessing import Pool

def flocking(args):
    T = args[0]
    J = args[1]
    rep = args[2]
    base_folder = '/net/levsha/share/simongh/sims/flocking/sweep_Jordan'
    # base_folder = '/home/simongh/simulations/flocking/dense_sweep_dt0.00001'
    filename = os.path.join(base_folder, '{2}_T={0}_J={1}.h5'.format(T, J, rep))
    cmd = './c++/flocking -N 128 -o {2} -dt 1e-4 -v0 0.088 -doTheory -ae 1000 -runtime 500 -T {0} -J {1}'.format(T, J, filename)
    os.system(cmd)

if __name__ == '__main__':
    Ts = np.logspace(-3, 2, 50);
    Js = np.array([1]);#np.logspace(-3, 3, 30);
    rep = 0;
    paramlist = []
    for i in range(15):
        for T in Ts:
            for J in Js:
                paramlist.append((T, J, rep))
                rep += 1

    with Pool(32) as mypool:
        mypool.map(flocking, paramlist)
