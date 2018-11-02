from multiprocessing import Pool
import oct2py
import tempfile
import shutil
import time
import sys

def fun(p):
    import oct2py
    import tempfile
    import shutil
    tmpdir = tempfile.mkdtemp()
    oc=oct2py.Oct2Py()
    ret = oc.antenna_sim(p, tmpdir)
    oc.exit()
    shutil.rmtree(tmpdir)
    return ret

if __name__ == '__main__':
    cost_per_hour = float(sys.argv(1))
    pool = Pool()
    params = [[0.0054383338063546697, 0.015233551609548139, 0.011012782095881371, 0.012108818391510621] for _ in range(pool._processes)]
    t0 = time.time()
    pool.map(fun, params)
    t1 = time.time()
    print("score: %f" % ((t1-t0)/(pool._processes*cost_per_hour),))
