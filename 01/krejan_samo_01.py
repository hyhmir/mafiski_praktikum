# imports
import numpy as np
from mpmath import mp


# constants and presets

mp.dps = 60
alpha = mp.mpf('0.355028053887817239')
beta = mp.mpf('0.258819403792806798')


# funcs

def f(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        term *= x**3 / (3*i * (3*i - 1))
        out += term
    return out


def g(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf(x)
    term = mp.mpf(x)
    for i in range(1, n + 1):
        term *= x**3 / (3*i * (3*i + 1))
        out += term
    return out


def L(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        multiply = (3*i + mp.mpf('5/2')) * (3*i + mp.mpf('1/2')) / ((i + 1) * 18 * x)
        if np.abs(multiply) > 1:
            break
        term *= multiply
        out += term
    return out


def P(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        multiply = - (6*i - mp.mpf('5/2')) * (3*i - mp.mpf('1/2')) * (6*i - mp.mpf('7/2')) * (3*i - mp.mpf('11/2')) / ((2*i - 1) * 2*i * 18**2 * x**2)
        if np.abs(multiply) > 1:
            break
        term *= multiply
        out += term
    return out


def Q(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf('5/72') / mp.mpf(x)
    term = mp.mpf('5/72') / mp.mpf(x)
    for i in range(1, n + 1):
        multiply = - (6*i - mp.mpf('5/2')) * (3*i - mp.mpf('1/2')) * (6*i + mp.mpf('5/2')) * (3*i + mp.mpf('1/2')) / ((2*i + 1) * 2*i * 18**2 * x**2)
        if np.abs(multiply) > 1:
            break
        term *= multiply
        out += term
    return out

print(alpha*f(mp.mpf(5), 100000) - beta*g(mp.mpf(5), 100000))
print(mp.airyai(5))


# graphs
