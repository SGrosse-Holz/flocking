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
    # base_folder = '/net/levsha/share/simongh/sims/flocking/sweep'
    base_folder = '/home/simongh/simulations/flocking/dense_sweep_dt0.1'
    filename = os.path.join(base_folder, '{2}_T={0}_J={1}.h5'.format(T, J, rep))
    cmd = './c++/flocking -runtime 200 -ae 1 -T {0} -J {1} -o {2} -dt_save 0 -dt 0.1'.format(T, J, filename)
    os.system(cmd)

if __name__ == '__main__':
    Ts = np.logspace(-2, 4, 30);
    Js = np.logspace(-4, 4, 30);
    rep = 0
    paramlist = []
    for i in range(10):
        for T in Ts:
            for J in Js:
                paramlist.append((T, J, rep))
                rep += 1

    with Pool(32) as mypool:
        mypool.map(flocking, paramlist)
