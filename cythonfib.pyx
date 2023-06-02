# cythonfib.pyx

""" Example cython interface definition """


cdef extern from "cfib.hpp":
    int fib(int n)

def pfib(n):
    return fib(n) 

cdef extern from "cfib.hpp":
    int fibb(int n)

def pfibb(n):
    return fibb(n) 

