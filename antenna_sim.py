import oct2py
import numpy as np
import blackbox as bb
from math import *
import tempfile

def fun(p):
    tmpdir = tempfile.mkdtemp()
    print('evaluating objective function with p = %s (%s)' % (str(p),tmpdir))
    oc=oct2py.Oct2Py()
    ret = oc.antenna_sim(p, tmpdir)
    print('f(%s) = %f (%s)' % (str(p),ret,tmpdir))
    oc.exit()
    return ret

if __name__ == '__main__':
    C0 = 299792458
    f_0 = 6.4896e9
    lambda0 = C0/f_0

    bounds = [
        [0.33*lambda0/(2*pi), 3.*lambda0/(2*pi)], # radius
        [3e-3, 2e-2], #length
        [lambda0*.05, lambda0*.5], #pitch
        [0.,2.] #taper ratio
        ]

    bb.search(
        f=fun,  # given function
        box=bounds,  # range of values for each parameter
        n=150,  # number of function calls on initial stage (global search)
        m=50,  # number of function calls on subsequent stage (local search)
        batch=4,  # number of calls that will be evaluated in parallel
        resfile='output.csv') # text file where results will be saved
