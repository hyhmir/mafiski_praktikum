# imports
import numpy as np
import mpmath
from mpmath import mp
import matplotlib.pyplot as plt
import time


# constants and presets

mp.dps = 500
alpha = mp.mpf('0.35502805388781723926006318600418317639797917419917724058332651030081004245')
beta = mp.mpf('0.25881940379280679840518356018920396347909113835493458221000181385610277267')


# funcs


def f(x, n):
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        term *= mp.mpf(x)**3 / (3*i * (3*i - 1))
        out += term
    return out


def g(x, n):
    out = mp.mpf(x)
    term = mp.mpf(x)
    for i in range(1, n + 1):
        term *= mp.mpf(x)**3 / (3*i * (3*i + 1))
        out += term
    return out


def L(x, n):
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        new = term * (3 * i - mp.mpf('5/2')) * (3 * i - mp.mpf('1/2')) / (18 * x * i)
        if np.abs(new) > np.abs(term):
            return out, i - 1
        term = new
        out += term
    return out, n


def P(x, n):
    out = mp.mpf('1')
    term = mp.mpf('1')
    for i in range(1, n + 1):
        new = term * (6 * i - mp.mpf('1/2')) * (6 * i - mp.mpf('5/2')) * (6 * i - mp.mpf('7/2')) * (6 * i - mp.mpf('11/2')) / (mp.mpf(-324) * x**2 * 2 * i * (2 * i - 1))
        if np.abs(new) > np.abs(term):
            return out, i - 1
        term = new
        out += term
    return out, n


def Q(x, n):
    out = mp.mpf('5/72') / mp.mpf(x)
    term = out
    for i in range(1, n + 1):
        new = term * (6 * i - mp.mpf('1/2')) * (6 * i - mp.mpf('5/2')) * (6 * i + mp.mpf('1/2')) * (6 * i + mp.mpf('5/2')) / (mp.mpf(-324) * x**2 * 2 * i * (2 * i + 1))
        if np.abs(new) > np.abs(term):
            return out, i - 1
        term = new
        out += term
    return out, n


def Ai_mac(x, n):
    return alpha * f(x, n) - beta * g(x, n), n


def Ai_neg(x, n):
    ksi = mp.mpf('2/3') * mp.sqrt(mpmath.fabs(x) ** (mp.mpf('3')))
    q = Q(ksi, n)
    p = P(ksi, n)
    T = mp.mpf(1) / (mp.sqrt(mp.pi) * mp.sqrt(mp.sqrt(-1*x)))
    fi = mp.pi / mp.mpf(4)
    return mpmath.re(T * (mp.sin(ksi - fi) * q[0] + mp.cos(ksi - fi) * p[0])), max(p[1], q[1])


def Ai_pos(x, n):
    ksi = mp.mpf('2/3') * mp.sqrt(mpmath.fabs(x) ** (mp.mpf('3')))
    l = L(- ksi, n)
    return mp.exp(- ksi) * l[0] / (mp.mpf(2) * mp.sqrt(mp.pi) * mp.sqrt(mp.sqrt(x))), l[1]


def Bi_mac(x, n):
    return mp.sqrt(3) * (alpha * f(x, n) + beta * g(x, n)), n


def Bi_neg(x, n):
    ksi = mp.mpf('2/3') * mp.sqrt(mpmath.fabs(x) ** (mp.mpf('3')))
    q = Q(ksi, n)
    p = P(ksi, n)
    T = mp.mpf(1) / (mp.sqrt(mp.pi) * mp.sqrt(mp.sqrt(-1*x)))
    fi = mp.pi / mp.mpf(4)
    return mpmath.re(T * (mp.cos(ksi - fi) * q[0] - mp.sin(ksi - fi) * p[0])), max(p[1], q[1])


def Bi_pos(x, n):
    ksi = mp.mpf('2/3') * mp.sqrt(mpmath.fabs(x) ** (mp.mpf('3')))
    l = L(ksi, n)
    return mp.exp(ksi) * l[0] / (mp.sqrt(mp.pi) * mp.sqrt(mp.sqrt(x))), l[1]


def absolute(x, func, ref, abs_err, max=500):
    n = 1
    i = n
    while n < max:
        f = func(x, n)
        err = np.abs(f[0] - ref)
        i = f[1]
        if err <= abs_err:
            break
        n += 1
    start = time.perf_counter()
    f = func(x, i)
    end = time.perf_counter()
    timer = end - start
    return (f[0], f[1], timer, err)


def relative(x, func, ref, rel_err, max=500):
    n = 1
    i = n
    while n < max:
        f = func(x, n)
        err = np.abs(f[0] - ref) / np.abs(ref)
        i = f[1]
        if err <= rel_err:
            break
        n += 1
    start = time.perf_counter()
    f = func(x, i)
    end = time.perf_counter()
    timer = end - start
    return (f[0], f[1], timer, err)


# graphs

x_neg = np.arange(-45, -5, 0.1)
x_mac = np.arange(-20, 20, 0.1)
x_pos = np.arange(5, 45, 0.1)


y_neg_abs = np.zeros((len(x_neg), 4))
y_mac_abs = np.zeros((len(x_mac), 4))
y_pos_abs = np.zeros((len(x_pos), 4))

y_neg_abs = np.zeros((len(x_neg), 4))
y_mac_abs = np.zeros((len(x_mac), 4))
y_pos_abs = np.zeros((len(x_pos), 4))
