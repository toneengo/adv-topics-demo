import numpy as np
from scipy.special import binom
import random
import math
#def f(x): return np.cos(20 * np.pi * (x**4)) + 1 - 2 * x
#def fprime(x): return (-80 * np.pi * (x**3) * np.sin(20 * np.pi * (x**4)) - 2)
DFS_LIMIT = 35
LEVEL_2D = 30
SEED = 3
DEGREE = 10

def bernstein_eqn(n, i, x):
    return binom(n, i) * (x**i) * (1 - x)**(n - i)


def bernstein_eqn_sum_value(n, m, x):
    y = 0
    for j in m:
        y += bernstein_eqn(n, j, x)
    return y
# ###################################################################


def Bernstein(n, i):
    def fxn(x):
        return binom(n, i) * (x**i) * (1 - x)**(n - i)
    return fxn


def Bernstein_prime(n, i):
    def fxn(x):
        try:
            return binom(n, i) * (i * (x**(i - 1)) * ((1 - x) **
                                                      (n - i)) - (x**i) * (n - i) * (1 - x)**(n - i - 1))
        except ZeroDivisionError:
            return 0
    return fxn


n = DEGREE
for i in range(n + 1):
    globals()['fxn' + '_{}'.format(i)] = Bernstein(n, i)
    globals()['fxn_prime' + '_{}'.format(i)] = Bernstein_prime(n, i)


def eqn_creator(to_sum):
    def fxn_sum(x):
        y = 0
        for i in to_sum:
            y += globals()['fxn' + '_{}'.format(i)](x)
        return y

    def fxn_prime_sum(x):
        y_prime = 0
        for i in to_sum:
            y_prime += globals()['fxn_prime' + '_{}'.format(i)](x)
        return y_prime
    return fxn_sum, fxn_prime_sum


a = np.arange(0, DEGREE + 1)
random.Random(SEED).shuffle(a)
newarr = np.array_split(a, 4)

eqn_1, eqn_1_prime = eqn_creator(newarr[0])
eqn_2, eqn_2_prime = eqn_creator(newarr[1])
eqn_3, eqn_3_prime = eqn_creator(newarr[2])
eqn_4, eqn_4_prime = eqn_creator(newarr[3])

fs = [eqn_1, eqn_2, eqn_3, eqn_4]
fprimes = [eqn_1_prime, eqn_2_prime, eqn_3_prime, eqn_4_prime]

def f(x): return np.cos(20 * np.pi * (x**4)) + 1 - 2 * x
def fprime(x): return (-80 * np.pi * (x**3) * np.sin(20 * np.pi * (x**4)) - 2)

def f2(x): return np.sin(x) + 2
def f3(x): return np.cos(x) + 2

def f2prime(x): return np.cos(x)
def f3prime(x): return np.sin(x) * -1.0
fs2 = [f2, f3]
fprimes2 = [f2prime, f3prime]
