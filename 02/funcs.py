import numpy as np
from scipy.optimize import curve_fit
from numba import njit, prange
import jax
import jax.numpy as jnp

np.random.seed(1298456)


@njit(parallel=True)
def fun_interp(times, dist, time):
    n_walks = times.shape[0]
    n_time = len(time)
    dists = np.empty((n_walks, n_time))
    for i in prange(n_walks):
        dists[i] = np.interp(time, times[i], dist[i])
    return dists


def fun_pareto(mu, a, b, n):
    podatki = (
        np.random.uniform(0, 1, n) * (b ** (1 - mu) - a ** (1 - mu)) + a ** (1 - mu)
    ) ** (1 / (1 - mu))

def fun_mad_std(arr):
    return 1.4826 * np.median(np.abs(arr - np.median(arr, axis=0)), axis=0)


def gen_fi(n, walks, seed):
    # np.random.seed(seed)
    return np.random.uniform(0, 2 * np.pi, (walks, n))


def gen_l(n, walks, l, mu, seed):
    # np.random.seed(seed)
    return (np.random.pareto(mu - 1, (walks, n)) + 1) * l


def gen_flight(t, walks, l, mu, seed1, seed2):
    Fi = gen_fi(t, walks, seed1)
    L = gen_l(t, walks, l, mu, seed2)
    time = np.arange(1, t + 1)
    return np.cumsum(np.stack((L * np.cos(Fi), L * np.sin(Fi))), axis=-1), time


def gen_walk(t, walks, l, mu, seed1, seed2):
    Fi = gen_fi(t, walks, seed1)
    L = gen_l(t, walks, l, mu, seed2)
    time = np.cumsum(L, axis=-1)
    return np.cumsum(np.stack((L * np.cos(Fi), L * np.sin(Fi))), axis=-1), time


def gen_dist(n, walks, mu, flight):
    if flight:
        fly = gen_flight(n, walks, 1, mu, 10000, 10000)[0]
        # print(fly.shape)
        dists = np.linalg.norm(fly, axis=0)[:, 100:]
        # print(dists.shape)
        time = np.arange(100, n)
    else:
        walk, times = gen_walk(n, walks, 1, mu, 10000, 10000)
        # print(walk.shape)
        # print(times.shape)
        time = np.linspace(100, np.min(times[:,-1]), 500)
        # print(walk.shape)
        dist = np.linalg.norm(walk, axis=0)
        # print(dist.shape)
        dists = np.empty((walks, len(time)))
        dists = fun_interp(times, dist, time)
    return dists, time


def calc_gamma(n, walks, mu, flight):
    distances, times = gen_dist(n, walks, mu, flight)
    log_mad = 2 * np.log(fun_mad_std(distances))
    # print(distances.shape)
    log_time = np.log(times)
    # print(log_time.shape)

    def lin(x, a, b):
        return a * x + b
    
    popt, pcov = curve_fit(lin, log_time, log_mad)
    a, b = popt
    a_err, b_err = np.sqrt(np.diag(pcov))
    return a, a_err


def gen_sequence(flight, repetitions):

    mus = np.arange(1.1, 4., 0.1)
    gammas = np.empty((len(mus), repetitions))
    for i, mu in enumerate(mus):
        for j in range(repetitions):
            gammas[i, j] = calc_gamma(10000, 1000, mu, flight)[0]

    if flight:
        gammas = gammas ** (-1)
    gamma_means = np.mean(gammas, axis=1)
    gamma_std = np.std(gammas, axis=1, ddof=1)

    return mus, gamma_means, gamma_std




if __name__ == '__main__':
    print(gen_sequence(False, 10))


