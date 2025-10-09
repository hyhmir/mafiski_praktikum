# imports
import numpy as np
import mpmath
from mpmath import mp


# constants and presets

mp.dps = 60
alpha = mp.mpf('0.355028053887817239')
beta = mp.mpf('0.258819403792806798')


# funcs
# tudi vsi int-i so v mp.mpf saj mi je nekaj delalo tezave z natancnostjo in jih kasneje nisem spremenil nazaj. na performance ne vplivajo skoraj nic

def f(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        term *= mp.mpf(x)**mp.mpf('3') / (mp.mpf('3')*mp.mpf(i) * (mp.mpf('3')*mp.mpf(i) - mp.mpf('1')))
        out += term
    return out


def g(x: mp.mpf, n: int) -> mp.mpf:
    out = mp.mpf(x)
    term = mp.mpf(x)
    for i in range(1, n + 1):
        term *= mp.mpf(x)**mp.mpf('3') / (mp.mpf('3')*mp.mpf(i) * (mp.mpf('3')*mp.mpf(i) + mp.mpf('1')))
        out += term
    return out


def L(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        multiply = (mp.mpf(3)*mp.mpf(i) + mp.mpf('5/2')) * (mp.mpf(3)*mp.mpf(i) + mp.mpf('1/2')) / ((mp.mpf(i) + mp.mpf(1)) * mp.mpf(18) * mp.mpf(x))
        if np.abs(multiply) > 1:
            return out, i - 1
        term *= multiply
        out += term
    return out, n


def P(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        multiply = - (mp.mpf(6)*mp.mpf(i) - mp.mpf('5/2')) * (mp.mpf(3)*mp.mpf(i) - mp.mpf('1/2')) * (mp.mpf(6)*mp.mpf(i) - mp.mpf('7/2')) * (mp.mpf(3)*mp.mpf(i) - mp.mpf('11/2')) / ((mp.mpf(2)*mp.mpf(i) - mp.mpf(1)) * mp.mpf(2)*mp.mpf(i) * mp.mpf(18)**mp.mpf(2) * mp.mpf(x)**mp.mpf(2))
        if np.abs(multiply) > 1:
            return out, i - 1
        term *= multiply
        out += term
    return out, n


def Q(x: mp.mpf, n: int=100) -> mp.mpf:
    out = mp.mpf('5/72') / mp.mpf(x)
    term = mp.mpf('5/72') / mp.mpf(x)
    for i in range(1, n + 1):
        multiply = - (mp.mpf(6)*mp.mpf(i) - mp.mpf('5/2')) * (mp.mpf(3)*mp.mpf(i) - mp.mpf('1/2')) * (mp.mpf(6)*mp.mpf(i) + mp.mpf('5/2')) * (mp.mpf(3)*mp.mpf(i) + mp.mpf('1/2')) / ((mp.mpf(2)*mp.mpf(i) + mp.mpf(1)) * mp.mpf(2)*mp.mpf(i) * mp.mpf(18)**mp.mpf(2) * mp.mpf(x)**mp.mpf(2))
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
    return mpmath.re((mp.mpf(1) / (mp.sqrt(mp.pi) * (- mp.mpf(x)) ** (mp.mpf('1/4')))) * (mp.sin(ksi - (mp.pi / mp.mpf(4))) * q[0] + mp.cos(ksi - (mp.pi / mp.mpf(4))) * p[0])), max(p[1], q[1])


def Ai_pos(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    l = L(- ksi, n)
    return mp.exp(- ksi) * l[0] / (mp.mpf(2) * mp.sqrt(mp.pi) * mp.mpf(x) ** mp.mpf('1/4')), l[1]


def Bi_mac(x, n):
    return mp.sqrt(3) * (alpha * f(x, n) + beta * g(x, n)), n


def Bi_neg(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    q = Q(x, n)
    p = P(x, n)
    return mpmath.re((mp.mpf(1) / (mp.sqrt(mp.pi) * (- mp.mpf(x)) ** (mp.mpf('1/4')))) * (- mp.sin(ksi - (mp.pi / mp.mpf(4))) * p[0] + mp.cos(ksi - (mp.pi / mp.mpf(4))) * q[0])), max(p[1], q[1])


def Bi_pos(x, n):
    ksi = mp.mpf('2/3') * mpmath.fabs(x) ** (mp.mpf('3/2'))
    l = L(ksi, n)
    return mp.exp(ksi) * l[0] / (mp.mpf(2) * mp.sqrt(mp.pi) * mp.mpf(x) ** mp.mpf('1/4')), l[1]





print(Bi_mac(30, 100), mp.airybi(30))
# print(mp.mpf(3))

# graphs
