import oct2py
import numpy as np
import blackbox as bb
from math import *
import tempfile
import shutil

def fun(p):
    tmpdir = tempfile.mkdtemp()
    print('evaluating objective function with p = %s' % (str(p),))
    oc=oct2py.Oct2Py()
    ret = oc.antenna_sim(p, tmpdir)
    print('f(%s) = %e' % (str(p),ret,))
    oc.exit()
    shutil.rmtree(tmpdir)
    return ret

if __name__ == '__main__':
    C0 = 299792458
    f_0 = 6.4896e9
    lambda0 = C0/f_0

    bounds = [
        [.66*lambda0/(2*pi), 1.5*lambda0/(2*pi)], # radius
        [1e-2, 2e-2], #length
        [lambda0*.05, lambda0*.5], #pitch2
        [0.,lambda0*.5] #pitch1
        ]

    bb.search(
        f=fun,  # given function
        box=bounds,  # range of values for each parameter
        n=16,  # number of function calls on initial stage (global search)
        m=8,  # number of function calls on subsequent stage (local search)
        batch=4,  # number of calls that will be evaluated in parallel
        resfile='output.csv') # text file where results will be saved
