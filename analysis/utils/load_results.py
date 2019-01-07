import numpy as np
import collections
import mpmath

mpmath.mp.dps = 500
def pval(x): return float(mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2)))))
#pval(8.24400983)

res_fields = ['rsid', 'qscore', 'pval', 'mu', 'sigma']
class TejaasResult(collections.namedtuple('_TejaasResult', res_fields)):
    __slots__ = ()

    @property
    def logp(self):
        if self.pval > 0:
            res = np.log10(self.pval)
        elif self.pval == 0:
            res = pval((self.qscore - self.mu) / self.sigma)
        elif self.pval < 0:
            res = -2000
        return res

    def __repr__(self):
        parent_string = super(TejaasResult, self).__repr__().strip(')')
        return '%s, logp=%f)' %(parent_string, self.logp)

def tejaas(filepath):
    res = dict()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            p     = float(arr[5])
            q     = float(arr[2])
            mu    = float(arr[3])
            sigma = float(arr[4])
            logp  = np.log10(p) if p != 0 else pval( (q - mu) / sigma)
            res[rsid] = -logp
    return res

def jpa(filepath):
    res = dict()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr = line.strip().split("\t")
            rsid = arr[0]
            res[rsid] = float(arr[1])
    return res


def matrixeqtl(filepath):
    res = dict()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            pval = float(arr[4])
            if rsid not in res:
                res[rsid] = -np.log10(pval)
    return res

def matrixeqtl_fdr(filepath):
    res = dict()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            fdr  = float(arr[5])
            if rsid not in res:
                res[rsid] = -np.log10(fdr)
    return res
