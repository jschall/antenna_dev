import numpy as np
import blackbox as bb
from math import *
from dask.distributed import Client
from contextlib import contextmanager

def install(package):
    import pip
    pip.main(['install', package])

def fun(p):
    global mcode
    try:
        install('oct2py')
    except:
        pass
    import oct2py
    import tempfile
    import shutil
    import os
    fd, matfile = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'wb') as f:
        f.write(mcode)
        f.flush()
    tmpdir = tempfile.mkdtemp()
    print('evaluating objective function with p = %s' % (str(p),))
    oc=oct2py.Oct2Py()
    ret = oc.feval(matfile, p, tmpdir)
    print('f(%s) = %f' % (str(p),ret,))

    try:
        oc.exit()
    except:
        pass

    try:
        f.close()
    except:
        pass

    try:
        os.close(fd)
    except:
        pass

    try:
        os.remove(matfile)
    except:
        pass

    try:
        shutil.rmtree(tmpdir)
    except:
        pass
    return ret

if __name__ == '__main__':
    with open('antenna_sim.m','rb') as f:
        mcode = f.read()

    client = Client('tcp://10.1.10.78:8786')

    C0 = 299792458
    f_0 = 6.4896e9
    lambda0 = C0/f_0

    bounds = [
        [6.0429e-03/1.15, 6.0429e-03*1.15], # radius
        [1.0871e-02/1.15, 1.0871e-02*1.15], #length
        [4.4621e-03/1.15, 4.4621e-03*1.15], #pitch2
        [4.2901e-03/1.5, 4.2901e-03*1.5] #pitch1
        ]

    class Mapper:
        def map(self, *args, **kwargs):
            global client
            return client.gather(client.map(*args, **kwargs))
    mapper = Mapper()

    @contextmanager
    def executor():
        yield mapper

    bb.search(
        f=fun,  # given function
        box=bounds,  # range of values for each parameter
        n=80,  # number of function calls on initial stage (global search)
        m=24,  # number of function calls on subsequent stage (local search)
        batch=8,  # number of calls that will be evaluated in parallel
        resfile='output.csv',  # text file where results will be saved
        executor = executor
        )
