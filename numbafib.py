import numpy as np
from numba import jit

@jit(nopython=True)
def pfib(n):
    t0 = 0
    t1 = 1
     
    for i in range(2, n+1):
        t = t0 + t1
        t0 = t1
        t1 = t

    return t1
