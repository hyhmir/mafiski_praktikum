import numpy as np
from scipy.optimize import curve_fit
from numba import njit, prange
import jax
import jax.numpy as jnp

# np.random.seed(1298456)


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


def gen_fi(n, walks, wait=None):
    if wait == None:
        return np.random.uniform(0, 2 * np.pi, (walks, n))
    return np.random.uniform(0, 2 * np.pi, (walks, 2*n))


def gen_l(n, walks, L, mu, wait=None):
    if wait == None:
        return (np.random.pareto(mu - 1, (walks, n)) + 1) * L
    l = np.zeros((walks, 2*n))
    l[:,::2] = (np.random.pareto(mu - 1, (walks, n)) + 1) * L
    return l


def gen_flight(t, walks, l, mu):
    Fi = gen_fi(t, walks)
    L = gen_l(t, walks, l, mu)
    time = np.arange(1, t + 1)
    return np.cumsum(np.stack((L * np.cos(Fi), L * np.sin(Fi))), axis=-1), time


def gen_walk(t, walks, l, mu, wait=None):
    Fi = gen_fi(t, walks, wait)
    L = gen_l(t, walks, l, mu, wait)
    time = np.copy(L)
    if wait != None:
        time[:,1::2] += (np.random.pareto(wait - 1, (walks, t)) + 1)
    time = np.cumsum(time, axis=-1)
    return np.cumsum(np.stack((L * np.cos(Fi), L * np.sin(Fi))), axis=-1), time


def gen_dist(n, walks, mu, flight, wait=None):
    if flight:
        fly = gen_flight(n, walks, 1, mu)[0]
        # print(fly.shape)
        dists = np.linalg.norm(fly, axis=0)[:, 100:]
        # print(dists.shape)
        time = np.arange(100, n)
    else:
        walk, times = gen_walk(n, walks, 1, mu, wait)
        # print(walk.shape)
        # print(times.shape)
        time = np.linspace(100, np.min(times[:,-1]), 1000)
        # print(walk.shape)
        dist = np.linalg.norm(walk, axis=0)
        # print(dist.shape)
        dists = np.empty((walks, len(time)))
        dists = fun_interp(times, dist, time)
    return dists, time


def calc_gamma(n, walks, mu, flight, wait=None):
    distances, times = gen_dist(n, walks, mu, flight, wait)
    log_mad = 2 * np.log(fun_mad_std(distances))
    # print(distances.shape)
    log_time = np.log(times)
    # print(log_time.shape)
    (a, b), cov = np.polyfit(log_time, log_mad, 1, cov=True)

    return a, np.sqrt(np.diag(cov))[0]

# @njit(parallel=True)
def gen_sequence(flight, repetitions, wait=None):
    if wait == None:
        mus = np.arange(1.1, 4., 0.1)
        gammas = np.empty((len(mus), repetitions))
        for i, mu in enumerate(mus):
            for j in range(repetitions):
                gammas[i, j] = calc_gamma(10000, 100, mu, flight, wait)[0]

        if flight:
            gammas = gammas ** (-1)
        gamma_means = np.mean(gammas, axis=1)
        gamma_std = np.std(gammas, axis=1, ddof=1)

        return mus, gamma_means, gamma_std, 1
    mus = np.arange(1.1, 4., 0.1)
    nis = np.arange(1.1, 4., 0.1)
    gammas = np.empty((len(mus), len(nis), repetitions))
    for i, mu in enumerate(mus):
        for j, ni in enumerate(nis):
            for k in range(repetitions):
                gammas[i, j, k] = calc_gamma(10000, 1000, mu, flight, ni)[0]

    gamma_means = np.mean(gammas, axis=2)
    gamma_std = np.std(gammas, axis=2, ddof=1)

    return mus, nis, gamma_means, gamma_std



def gen_phase():
    dist2, time2 = gen_dist(10000, 100, 2, False)
    dist3, time3 = gen_dist(10000, 100, 3, False)
    mads2 = fun_mad_std(dist2)
    mads3 = fun_mad_std(dist3)
    def fun2(t, a):
        return a*t**2/np.log(t)
    
    def fun3(t, a):
        return a*t/np.log(t)
    
    popt2, cov2 = curve_fit(fun2, time2, mads2**2)
    a2 = popt2
    a2_err = np.sqrt(np.diag(cov2))
    popt3, cov3 = curve_fit(fun3, time3, mads3**2)
    a3 = popt3
    a3_err = np.sqrt(np.diag(cov3))

    return (a2, a2_err, mads2**2, time2), (a3, a3_err, mads3**2, time3)



if __name__ == '__main__':
    print(gen_sequence(False, 2, 1))


