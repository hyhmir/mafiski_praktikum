# imports
import numpy as np
import mpmath
from mpmath import mp
import matplotlib.pyplot as plt


# constants and presets

mp.dps = 500
alpha = mp.mpf('0.355028053887817239')
beta = mp.mpf('0.258819403792806798')


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


# graphs

x = np.arange(-30, -10, 0.1)

y = [Ai_neg(i, 100)[0] for i in x]

y1 = [float(i) for i in y]
y2 = [float(mp.airyai(i)) for i in x]

y3 = np.abs([float(mp.airyai(i) - Ai_neg(mp.mpf(i), 100)[0]) for i in x])

# y4 = [float(Ai_neg(mp.mpf(i), 1000)[1]) for i in x]


plt.plot(x, y1)
plt.plot(x, y2)
# plt.yscale('log')
plt.show()
plt.plot(x, y3)
plt.yscale('log')
plt.show()
# plt.plot(x, y4)
# plt.yscale('log')
# plt.show()

