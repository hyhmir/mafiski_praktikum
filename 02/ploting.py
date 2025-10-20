import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from funcs import *

np.random.seed(123456)
plt.style.use(['science'])



def save(name='', xlabel='x', ylabel='y', legend=True, grid=True):
    if grid:
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

    plt.savefig(f'02/graphs/{name.replace(' ', '').replace('$', '').replace('\\', '')}.pdf', dpi=512, bbox_inches='tight')
    plt.clf()

# ### ploting walks ###

for n in [10, 100, 1000, 10000]:
    walks10, t = gen_flight(n, 3, 1, 2)

    for i in range(3):
        x = walks10[0, i] - walks10[0, i, 0]
        y = walks10[1, i] - walks10[1, i, 0]
        plt.plot(x, y)
    save(f'$N = {n},\\ \\mu = 2$', legend=False)

for n in [10, 100, 1000, 10000]:
    walks10, t = gen_flight(n, 3, 1, 3)

    for i in range(3):
        x = walks10[0, i] - walks10[0, i, 0]
        y = walks10[1, i] - walks10[1, i, 0]
        plt.plot(x, y)
    save(f'$N = {n},\\ \\mu = 3$', legend=False)


### primer dolocitve gama ###
dists, times = gen_dist(100000, 100, 2.5, False)
print(dists.shape)
log_mad = 2*np.log(fun_mad_std(dists))
log_time = np.log(times)
a, b = np.polyfit(log_time, log_mad, 1)
plt.scatter(log_time, log_mad, s=1, label=f'Povprecje $\\sigma$')
plt.plot(log_time, a*log_time+b, label=f'y = {round(a, 4)}x  {round(b, 4)}')
save('Primer dolocitve $\\gamma$', xlabel='$\\ln t$', ylabel='$2 \\ln \\sigma$')


mus, gamma, gamma_err, c = gen_sequence(True, 20)
def expect3(mu):
    if mu<3:
        return (mu - 1) / 2
    return 1.
e_gamma = []
for mu in mus:
    e_gamma.append(expect3(mu))
plt.errorbar(mus, gamma, gamma_err, fmt='.', capsize=1, label='Simulacija')
plt.plot(mus, e_gamma, label='Teorija')
save('poleti', '$\\mu$', '$\\gamma^{-1}$')


mus, gamma, gamma_err, c = gen_sequence(False, 20)
def expect2(mu):
    if mu<2:
        return 2.
    elif mu<3:
        return 4 - mu
    return 1.
e_gamma = []
for mu in mus:
    e_gamma.append(expect2(mu))
plt.errorbar(mus, gamma, gamma_err, fmt='.', capsize=1, label='Simulacija')
plt.plot(mus, e_gamma, label='Teorjia')
save('sprehodi', '$\\mu$', '$\\gamma$')


mu2, mu3 = gen_phase()

plt.plot(mu2[3], mu2[2], label='Povprecje $\\sigma^2$')
plt.plot(mu2[3], mu2[0]*mu2[3]**2/np.log(mu2[3]), label='pricakovano')
save('$\\mu$ = 2', '$t$', '$\\sigma^2$')


plt.plot(mu3[3], mu3[2], label='Povprecje $\\sigma^2$')
plt.plot(mu3[3], mu3[0]*mu3[3]/np.log(mu3[3]), label='pricakovano')
save('$\\mu$ = 3', '$t$', '$\\sigma^2$')


def expect1(mu, ni):
    if mu<=2 and mu<=ni:
        return 2.
    elif 2<mu<=3 and ni>=2:
        return 4 - mu
    elif mu>3 and ni>=2:
        return 1.
    elif ni<2:
        if ni>-1 + mu:
            return 2 + ni - mu
        return ni - 1
    
mus, nis, gamma, c = gen_sequence(False, 2, 1)

e_gamma = np.zeros((len(mus), len(nis)))
for i,mu in enumerate(mus):
    for j,ni in enumerate(nis):
        e_gamma[i, j] = expect1(mu, ni)

plt.pcolormesh(mus, nis, e_gamma)
plt.colorbar()
save('Pričakovani režimi', '$\\nu$', '$\\mu$', False, False)


plt.pcolormesh(mus, nis, gamma)#, vmax=3, vmin=0)
plt.colorbar()
save('Dejanski režimi', '$\\nu$', '$\\mu$', False, False)


def e_15(x):
    if x <= 1.5:
        return 2.
    if x >= 2.5:
        return 0.5
    return -1.5*x + 4.25

exp = []

mus = np.arange(1.1, 4., 0.1)
for mu in mus:
    exp.append(e_15(mu))
mi = []
mi_err = []
for mu in mus:
    m, merr = calc_gamma(10000, 100, mu, False, 1.5)
    mi.append(m)
    mi_err.append(merr)

plt.errorbar(mus, mi, mi_err, fmt='.', capsize=1, label='Simulacija')
plt.plot(mus, exp, label='Napoved')
save('$\\nu$ = 1.5', '$\\mu$', '$\\gamma$')

def e_23(x):
    if x <= 2:
        return 2.
    if x >= 3:
        return 1
    return -x + 4

exp = []

mus = np.arange(1.1, 4., 0.1)
for mu in mus:
    exp.append(e_23(mu))
mi = []
mi_err = []
for mu in mus:
    m, merr = calc_gamma(10000, 100, mu, False, 2)
    mi.append(m)
    mi_err.append(merr)

plt.errorbar(mus, mi, mi_err, fmt='.', capsize=1, label='Simulacija')
plt.plot(mus, exp, label='Napoved')
save('$\\nu$ = 2', '$\\mu$', '$\\gamma$')

exp = []

mus = np.arange(1.1, 4., 0.1)
for mu in mus:
    exp.append(e_23(mu))
mi = []
mi_err = []
for mu in mus:
    m, merr = calc_gamma(10000, 100, mu, False, 3)
    mi.append(m)
    mi_err.append(merr)

plt.errorbar(mus, mi, mi_err, fmt='.', capsize=1, label='Simulacija')
plt.plot(mus, exp, label='Napoved')
save('$\\nu$ = 3', '$\\mu$', '$\\gamma$')

