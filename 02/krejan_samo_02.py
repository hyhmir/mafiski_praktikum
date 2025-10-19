import numpy as np


def gen_fi(n, seed):
    np.random.seed(seed)
    return np.random.uniform(0, 2 * np.pi, n)


def gen_l(n, l, mu, seed):
    np.random.seed(seed)
    return (np.random.pareto(mu - 1, n) + 1) * l


def gen_flight(t, l, mu, seed1, seed2):
    Fi = gen_fi(t, seed1)
    L = gen_l(t, l, mu, seed2)
    time = np.arange(1, t + 1)
    return np.cumsum(np.column_stack((L * np.cos(Fi), L * np.sin(Fi))), axis=0), time


def gen_walk(t, l, mu, seed1, seed2):
    Fi = gen_fi(t, seed1)
    L = gen_l(t, l, mu, seed2)
    time = np.cumsum(L)
    return np.cumsum(np.column_stack((L * np.cos(Fi), L * np.sin(Fi))), axis=0), time

print(gen_walk(10, 1, 2, 12345, 15445))