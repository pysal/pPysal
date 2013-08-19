import pysal
import multiprocessing as mp
import numpy as np
from pysal.weights.spatial_lag import lag_spatial as slag
from pysal.common import *
import ctypes

PERMUTATIONS = 999

class Moran:
    """Moran's I Global Autocorrelation Statistic

    Parameters
    ----------

    y               : array
                      variable measured across n spatial units
    w               : W
                      spatial weights instance
    transformation  : string
                      weights transformation,  default is row-standardized "r".
                      Other options include "B": binary,  "D":
                      doubly-standardized,  "U": untransformed (general weights),
                      "V": variance-stabilizing.
    permutations    : int
                      number of random permutations for calculation of pseudo-p_values


    Attributes
    ----------
    y            : array
                   original variable
    w            : W
                   original w object
    permutations : int
                   number of permutations
    I            : float
                   value of Moran's I
    EI           : float
                   expected value under normality assumption
    VI_norm      : float
                   variance of I under normality assumption
    seI_norm     : float
                   standard deviation of I under normality assumption
    z_norm       : float
                   z-value of I under normality assumption
    p_norm       : float
                   p-value of I under normality assumption (one-sided)
                   for two-sided tests, this value should be multiplied by 2
    VI_rand      : float
                   variance of I under randomization assumption
    seI_rand     : float
                   standard deviation of I under randomization assumption
    z_rand       : float
                   z-value of I under randomization assumption
    p_rand       : float
                   p-value of I under randomization assumption (1-tailed)
    sim          : array (if permutations>0)
                   vector of I values for permutated samples
    p_sim        : array (if permutations>0)
                   p-value based on permutations (one-sided)
                   null: spatial randomness
                   alternative: the observed I is extreme
                                it is either extremely greater or extremely lower
    EI_sim       : float (if permutations>0)
                   average value of I from permutations
    VI_sim       : float (if permutations>0)
                   variance of I from permutations
    seI_sim      : float (if permutations>0)
                   standard deviation of I under permutations.
    z_sim        : float (if permutations>0)
                   standardized I based on permutations
    p_z_sim      : float (if permutations>0)
                   p-value based on standard normal approximation from

    Examples
    --------
    >>> import pysal
    >>> w = pysal.open(pysal.examples.get_path("stl.gal")).read()
    >>> f = pysal.open(pysal.examples.get_path("stl_hom.txt"))
    >>> y = np.array(f.by_col['HR8893'])
    >>> mi = Moran(y,  w)
    >>> "%7.5f" % mi.I
    '0.24366'
    >>> mi.EI
    -0.012987012987012988
    >>> mi.p_norm
    0.00027147862770937614

    SIDS example replicating OpenGeoda

    >>> w = pysal.open(pysal.examples.get_path("sids2.gal")).read()
    >>> f = pysal.open(pysal.examples.get_path("sids2.dbf"))
    >>> SIDR = np.array(f.by_col("SIDR74"))
    >>> mi = pysal.Moran(SIDR,  w)
    >>> "%6.4f" % mi.I
    '0.2477'
    >>> mi.p_norm
    0.0001158330781489969

    """
    def __init__(self, y, w, transformation="r", permutations=PERMUTATIONS):
        self.y = y
        w.transform = transformation
        self.w = w
        self.permutations = permutations
        self.__moments()
        self.I = self.__calc(self.z)
        self.z_norm = (self.I - self.EI) / self.seI_norm
        self.p_norm = 2.0 * (1 - stats.norm.cdf(np.abs(self.z_norm)))
        self.z_rand = (self.I - self.EI) / self.seI_rand
        self.p_rand = 2.0 * (1 - stats.norm.cdf(np.abs(self.z_rand)))

        sim = [self.__calc(np.random.permutation(self.z))
                for i in xrange(permutations)]

        self.sim = sim = np.array(sim)
        above = sim >= self.I
        larger = sum(above)
        if (self.permutations - larger) < larger:
            larger = self.permutations - larger
        self.p_sim = (larger + 1.) / (permutations + 1.)
        self.EI_sim = sum(sim) / permutations
        self.seI_sim = np.array(sim).std()
        self.VI_sim = self.seI_sim ** 2
        self.z_sim = (self.I - self.EI_sim) / self.seI_sim
        self.p_z_sim = 2.0 * (1 - stats.norm.cdf(np.abs(self.z_sim)))

    def __moments(self):
        self.n = len(self.y)
        y = self.y
        #z = (y-y.mean())/y.std()
        z = y - y.mean()
        self.z = z
        self.z2ss = sum(z * z)
        self.EI = -1. / (self.n - 1)
        n = self.n
        s1 = self.w.s1
        s0 = self.w.s0
        s2 = self.w.s2
        s02 = s0 * s0
        v_num = n * n * s1 - n * s2 + 3 * s0 * s0
        v_den = (n - 1) * (n + 1) * s0 * s0
        self.VI_norm = v_num / v_den - (1.0 / (n - 1)) ** 2
        self.seI_norm = self.VI_norm ** (1 / 2.)

        k = (1 / (sum(z ** 4)) * ((sum(z ** 2)) ** 2))
        vi = (1 / (((n - 1) ** 3) * s02)) * ((n * ((n * n - 3 * n + 3) * s1 - n * s2 + 3 * s02))
                                             - (k * ((n * n - n) * s1 - 2 * n * s2 + 6 * s02)))
        self.VI_rand = vi
        self.seI_rand = vi ** (1 / 2.)

    def __calc(self, z):
        zl = slag(self.w, z)
        inum = sum(z * zl)
        return self.n / self.w.s0 * inum / self.z2ss


def moran_mp(y,w,transformation='r',permutations=PERMUTATIONS):
    w.transform = transformation
    n,z2ss,z,EI,seI_norm,seI_rand   = get_moments(y,w)
    I = _calc(z,w,n,z2ss)
    z_norm = (I - EI) / seI_norm
    p_norm = 2.0 * (1 - stats.norm.cdf(np.abs(z_norm)))
    z_rand = (I - EI) / seI_rand
    p_rand = 2.0 * (1 - stats.norm.cdf(np.abs(z_rand)))

    global c_perm
    c_perm = mp.RawArray(ctypes.c_double, permutations)

    cores = mp.cpu_count()
    pool = mp.Pool()
    step = permutations / cores
    for start in range(0,permutations,step):
        pool.apply_async(calc, args=(z,w,n,z2ss,start,step))
    pool.close()
    pool.join()
    sim = np.frombuffer(c_perm)

    above = sim >= I
    larger = sum(above)
    if (permutations - larger) < larger:
        larger = permutations - larger
    p_sim = (larger + 1.) / (permutations + 1.)
    EI_sim = sum(sim) / permutations
    seI_sim = np.array(sim).std()
    VI_sim = seI_sim ** 2
    z_sim = (I - EI_sim) / seI_sim
    p_z_sim = 2.0 * (1 - stats.norm.cdf(np.abs(z_sim)))

    return I

def get_moments(y,w):
    n = len(y)
    z = y - y.mean()
    z2ss = sum(z*z)
    EI = -1 / (n-1)
    s1 = w.s1
    s0 = w.s0
    s2 = w.s2
    s02 = s0 * s0

    v_num = n * n * s1 - n * s2 + 3 * s0 * s0
    v_den = (n - 1) * (n + 1) * s0 * s0

    VI_norm = v_num / v_den - (1.0 / (n - 1)) ** 2
    seI_norm = VI_norm ** (1 / 2.)

    k = (1 / (sum(z ** 4)) * ((sum(z ** 2)) ** 2))
    vi = (1 / (((n - 1) ** 3) * s02)) * ((n * ((n * n - 3 * n + 3) * s1 - n * s2 + 3 * s02))
            - (k * ((n * n - n) * s1 - 2 * n * s2 + 6 * s02)))

    VI_rand = vi
    seI_rand = vi ** (1 / 2.)

    return n,z2ss,z,EI,seI_norm,seI_rand

def _calc(z,w,n,z2ss):
    zl = slag(w, z)
    inum = sum(z * zl)
    return n / w.s0 * inum / z2ss


def calc(z,w,n,z2ss,start,stop):
    pid=mp.current_process()._identity[0]
    shared_sim = np.frombuffer(c_perm)
    for i in range(start,start+stop):
        r_num = np.random.RandomState(pid + i)
        z = np.random.permutation(z)
        zl = slag(w, z)
        inum = sum(z * zl)
        shared_sim[i] =  n / w.s0 * inum / z2ss

if __name__ == "__main__":
    global w
    w = pysal.open(pysal.examples.get_path("stl.gal")).read()
    f = pysal.open(pysal.examples.get_path("stl_hom.txt"))
    y = np.array(f.by_col['HR8893'])

    ###999###
    t1 = time.time()
    mi = Moran(y,  w, 'r', 999)
    t2 = time.time()
    print "PySAL serial (999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 999)
    t2 = time.time()
    print "PySAL MP (999): ", t2-t1
    assert(mi.I == mp_i)

    ###9999###
    t1 = time.time()
    mi = Moran(y,  w, 'r', 9999)
    t2 = time.time()
    print "PySAL serial (9999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 9999)
    t2 = time.time()
    print "PySAL MP (9999): ", t2-t1
    assert(mi.I == mp_i)

    ###19999###
    t1 = time.time()
    mi = Moran(y,  w, 'r', 19999)
    t2 = time.time()
    print "PySAL serial (19999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 19999)
    t2 = time.time()
    print "PySAL MP (19999): ", t2-t1
    assert(mi.I == mp_i)

    ###29999###
    t1 = time.time()
    mi = Moran(y,  w, 'r', 29999)
    t2 = time.time()
    print "PySAL serial (29999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 29999)
    t2 = time.time()
    print "PySAL MP (29999): ", t2-t1
    assert(mi.I == mp_i)

