# imports
import numpy as np
import mpmath
from mpmath import mp


# constants and presets

mp.dps = 300
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
            return out, i - 1
        term *= multiply
        out += term
    return out, n


def P(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        multiply = - (6*i - mp.mpf('5/2')) * (3*i - mp.mpf('1/2')) * (6*i - mp.mpf('7/2')) * (3*i - mp.mpf('11/2')) / ((2*i - 1) * 2*i * 18**2 * x**2)
        if np.abs(multiply) > 1:
            return out, i - 1
        term *= multiply
        out += term
    return out, n


def Q(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('5/72') / mp.mpf(x)
    term = mp.mpf('5/72') / mp.mpf(x)
    for i in range(1, n + 1):
        multiply = - (6*i - mp.mpf('5/2')) * (3*i - mp.mpf('1/2')) * (6*i + mp.mpf('5/2')) * (3*i + mp.mpf('1/2')) / ((2*i + 1) * 2*i * 18**2 * x**2)
        if np.abs(multiply) > 1:
            return out, i - 1
        term *= multiply
        out += term
    return out, n


def Ai_mac(x, n):
    return alpha * f(x, n) - beta * g(x, n), n


def Ai_neg(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    q = Q(x, n)
    p = P(x, n)
    return mpmath.re((1 / (mp.sqrt(mp.pi) * (- x) ** (mp.mpf('1/4')))) * (mp.sin(ksi - (mp.pi / mp.mpf(4))) * q[0] + mp.cos(ksi - (mp.pi / mp.mpf(4))) * p[0])), max(p[1], q[1])


def Ai_pos(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    l = L(- ksi, n)
    return mp.exp(- ksi) * l[0] / (mp.mpf(2) * mp.sqrt(mp.pi) * x ** mp.mpf('1/4')), l[1]


def Bi_mac(x, n):
    return mp.sqrt(3) * (alpha * f(x, n) + beta * g(x, n)), n


def Bi_neg(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    q = Q(x, n)
    p = P(x, n)
    return mpmath.re((1 / (mp.sqrt(mp.pi) * (- x) ** (mp.mpf('1/4')))) * (- mp.sin(ksi - (mp.pi / mp.mpf(4))) * p[0] + mp.cos(ksi - (mp.pi / mp.mpf(4))) * q[0])), max(p[1], q[1])


def Bi_pos(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    l = L(ksi, n)
    return mp.exp(ksi) * l[0] / (mp.mpf(2) * mp.sqrt(mp.pi) * x ** mp.mpf('1/4')), l[1]





print(Ai_mac(50, 100000), Ai_neg(50, 100000), Ai_pos(50, 100000), mp.airyai(50))


# graphs
