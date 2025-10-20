import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from funcs import *

np.random.seed(12345)
plt.style.use(['science'])



def save(name='', xlabel='x', ylabel='y', legend=True):
    plt.grid()
    if legend:
        plt.legend(loc='best',
            frameon=True,
            framealpha=0.9,
            facecolor='white',
            edgecolor='gray'
        )
    plt.title(name)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.savefig(f'02/graphs/{name}.pdf', dpi=512, bbox_inches='tight')
    plt.clf()

### ploting walks ###

# for n in [10, 100, 1000, 10000]:
#     walks10, t = gen_flight(n, 3, 1, 2)

#     for i in range(3):
#         x = walks10[0, i] - walks10[0, i, 0]
#         y = walks10[1, i] - walks10[1, i, 0]
#         plt.plot(x, y)
#     save(f'N = {n}, mu = 2', legend=False)

# for n in [10, 100, 1000, 10000]:
#     walks10, t = gen_flight(n, 3, 1, 3)

#     for i in range(3):
#         x = walks10[0, i] - walks10[0, i, 0]
#         y = walks10[1, i] - walks10[1, i, 0]
#         plt.plot(x, y)
#     save(f'N = {n}, mu = 3', legend=False)


### primer dolocitve gama ###
# dists, times = gen_dist(10000, 100, 2.5, False)
# print(dists.shape)
# log_mad = 2*np.log(fun_mad_std(dists))
# log_time = np.log(times)
# a, b = np.polyfit(log_time, log_mad, 1)
# plt.scatter(log_time, log_mad, s=1, label=f'Povprecje $sigma$')
# plt.plot(log_time, a*log_time+b, label=f'y = {round(a, 4)}x  {round(b, 4)}')
# save('Primer dolocitve $gama$', xlabel='lnt', ylabel='2ln sigma')


# mus, gamma, gamma_err, c = gen_sequence(True, 10)
# def expect(mu):
#     if mu<3:
#         return (mu - 1) / 2
#     return 1.
# e_gamma = []
# for mu in mus:
#     e_gamma.append(expect(mu))
# plt.errorbar(mus, gamma, gamma_err, fmt='.', capsize=1, label='Simulacija')
# plt.plot(mus, e_gamma, label='Teoreticna pricakovanja')
# save('poleti', 'mu', '1/gamma')


# mus, gamma, gamma_err, c = gen_sequence(False, 10)
# def expect(mu):
#     if mu<2:
#         return 2.
#     elif mu<3:
#         return 4 - mu
#     return 1.
# e_gamma = []
# for mu in mus:
#     e_gamma.append(expect(mu))
# plt.errorbar(mus, gamma, gamma_err, fmt='.', capsize=1, label='Simulacija')
# plt.plot(mus, e_gamma, label='Teoria')
# save('sprehodi', 'mu', 'gamma')


# mu2, mu3 = gen_phase()

# plt.plot(mu2[3], mu2[2], label='Povprecje sigma sq')
# plt.plot(mu2[3], mu2[0]*mu2[3]**2/np.log(mu2[3]), label='pricakovano')
# save('mu = 2', 't', 'sigma sq')


# mu2, mu3 = gen_phase()

# plt.plot(mu3[3], mu3[2], label='Povprecje sigma sq')
# plt.plot(mu3[3], mu3[0]*mu3[3]/np.log(mu3[3]), label='pricakovano')
# save('mu = 3', 't', 'sigma sq')


def expect(mu, ni):
    if ni<=2 and mu>=ni:
        return 2.
    elif 2<ni<=3 and mu>=2:
        return 4 - mu
    elif ni>3 and mu>=2:
        return 1.
    elif mu<2:
        if ni>=2+mu:
            return ni-1
        return 2 + ni - mu
    
mus, nis, gamma, c = gen_sequence(False, 10, 1)

# e_gamma = np.zeros((len(mus), len(nis)))
# for i,mu in enumerate(mus):
#     for j,ni in enumerate(nis):
#         e_gamma[i, j] = expect(mu, ni)

# plt.pcolormesh(mus, nis, e_gamma)
# plt.colorbar(label='Value')
# plt.title('Heatmap with pcolormesh')
# plt.xlabel('Y axis')
# plt.ylabel('X axis')
# plt.show()

plt.pcolormesh(mus, nis, gamma)
plt.colorbar(label='Value')
plt.title('Heatmap with pcolormesh')
plt.xlabel('Y axis')
plt.ylabel('X axis')
plt.show()
