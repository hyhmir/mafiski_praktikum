import rust as rs
import numpy as np
import matplotlib.pyplot as plt
import timeit
# import pandas as pd
# import time
from scipy.signal import correlate
import scienceplots


plt.style.use('science')


#### checking da rust koda dela ####
# t = np.linspace(0, 10, 1024)
# h = np.sin(t) + np.cos(2*t)

# plt.plot(rs.fft(h, True))
# plt.plot(np.fft.ifft(h))
# plt.show()

# plt.plot(np.fft.ifft(np.fft.fft(h)) - rs.fft(rs.fft(h, False), True))
# plt.show()

# plt.plot(rs.fft(rs.fft(h, False), True))
# plt.plot(h)
# plt.show()

# plt.plot(rs.fft(rs.fft(h, False), True) - h)
# plt.show()

#### funcs ####

def fft(x, fs, N=0, shift=False):
    if N == 0:
        try:
            N = x.shape[0]
        except:
            N = len(x)
    X = np.fft.fft(x, N)
    f = np.fft.fftfreq(N, 1/fs)
    if shift:
        X = np.fft.fftshift(X)
        f = np.fft.fftshift(f)
    return X, f


def ifft(X, fs, N=0):
    if N == 0:
        try:
            N = X.shape[0]
        except:
            N = len(X)
    x = np.fft.ifft(X, N)
    t = np.fft.fftfreq(N, 1/fs)
    return x, t


def pad_data(data):
    try:
        N = len(data)
    except:
        N = data.shape[0]
    pad = np.zeros(N)
    return np.hstack((data, pad))


def auto_corr(h):
    times = []
    t0 = timeit.default_timer()
    N = h.shape[0]
    ns = np.arange(N)
    fis = []
    for n in ns:
        A = 1 / (N - n)
        sum = 0
        for i in range(N - n):
            sum += h[i] * h[i + n]
        fis.append(A * sum)
        t1 = timeit.default_timer()
        times.append(t1 - t0)
        t0 = t1
    return fis, ns, times


def auto_corr_np(h):
    times = []
    t0 = timeit.default_timer()
    N = h.shape[0]
    ns = np.arange(N)
    fis = []
    h_pad = pad_data(h)
    for n in ns:
        A = 1 / (N - n)
        h_pad = pad_data(h)
        h_pad_roll = np.roll(h_pad, n)
        out = A * np.dot(h_pad, h_pad_roll)
        fis.append(out)
        t1 = timeit.default_timer()
        times.append(t1 - t0)
        t0 = t1

    return fis, ns, times

def auto_corr_fft(h, minus_n = True):
    times = []
    t0 = timeit.default_timer()
    N = h.shape[0]
    ns = np.arange(N)
    fis = []
    h_pad = pad_data(h)
    h_fft, h_freqs = fft(h_pad, 1)
    h_fft, h_freqs = ifft(np.abs(h_fft)**2, 1)
    for n in ns:
        A = 1 / (N)
        if minus_n:
            A = 1 / (N - n)
        out = A * h_fft[n]
        fis.append(out)
        t1 = timeit.default_timer()
        times.append(t1 - t0)
        t0 = t1

    return fis, ns, times

def auto_corr_fft_rel(h):
    ac_fft, ns_fft, t = auto_corr_fft(h, minus_n=False)
    h_mean = np.mean(h)
    return (ac_fft - h_mean**2) / (ac_fft[0] - h_mean**2), ns_fft



#### speeeeed #####

loglengs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
lengs = 2**loglengs
# print(lengs)

numpt = []
numperr = []

rustt = []
rusterr = []

for len in lengs:
    h = np.random.random(len)
    timesr = []
    timesn = []
    for i in range(10):
        timesn.append(timeit.timeit(lambda: np.fft.fft(h), number=1))
        timesr.append(timeit.timeit(lambda: rs.fft(h, False), number=1)) # type: ignore
    numpt.append(np.mean(timesn))
    numperr.append(np.std(timesn))
    rustt.append(np.mean(timesr))
    rusterr.append(np.std(timesr))

plt.errorbar(lengs, numpt, numperr, fmt='-', label='Numpy')
plt.errorbar(lengs, rustt, rusterr, fmt='-', label='Rust')
plt.legend()
plt.title('Hitrost FFT')
plt.xlabel('Dolžina signala')
plt.ylabel('t [s]')
plt.savefig('05/graphs/speed.pdf', dpi=512)
plt.clf()
# plt.show()


#### spektrogram ####

from scipy.io import wavfile

p1, bubo = wavfile.read('05/samples/bubomono.wav')
p2, bubo2 = wavfile.read('05/samples/bubo2mono.wav')
p3, mix = wavfile.read('05/samples/mix.wav')
p4, mix1 = wavfile.read('05/samples/mix1.wav')
p5, mix2 = wavfile.read('05/samples/mix2.wav')
p6, mix22 = wavfile.read('05/samples/mix22.wav')

data = [bubo, bubo2, mix, mix1, mix2, mix22]
sampling = [p1, p2, p3, p4, p5, p6]
names = ['bubo', 'bubo2', 'mix', 'mix1', 'mix2', 'mix22']



fig, axes = plt.subplots(3, 2, figsize=(10, 15))
axes = axes.flatten()

norm = 1

for i, ax in enumerate(axes):
    d = data[i] / norm
    sample_rate = sampling[i]
    pxx, freqs, bins, im = ax.specgram(d, Fs=sample_rate, NFFT=512)
    ax.set_title(f'Spektrogram {names[i]}.wav')
    ax.set_xlabel('Čas [s]\n\n')
    ax.set_ylabel('Frekvenca [Hz]')
    fig.colorbar(im, ax=ax)
    ax.set_ylim((0, 10000))
    # ax.heatmap(pxx[:pxx.shape[0]//2,])#, ax=ax)

plt.tight_layout()

plt.savefig('05/graphs/specgram.pdf', dpi=512)
plt.close(fig)


#### avtokorelacija ####


smpl_rate = 2240
n = np.linspace(0, 16*np.pi, smpl_rate)
sin = np.sin(n) + np.cos(3*n)  + np.random.pareto(5, smpl_rate) * np.random.rand(smpl_rate) * 3 * (-1)**np.random.randint(0, 2, smpl_rate)

po_def, f, t1 = auto_corr(sin)
po_def_np, f, t2 = auto_corr_np(sin)
po_fft, f, t3 = auto_corr_fft(sin)

fig, axes = plt.subplots(3, 1, figsize=(10, 8))
ax1, ax2, ax3 = axes.flatten()

ax1.plot(sin, label='signal + šum')
ax1.legend()
ax1.set_xlabel('$x$')
ax1.set_ylabel('Amplituda')
ax1.set_title('Tesni signal (zašumljeni harmoniki)')

ax2.plot(f, np.array(po_def), label='python', alpha=1)
ax2.plot(f, np.real(np.array(po_def_np)) + 1, label='numpy (+1)', alpha=1)
ax2.plot(f, np.real(np.array(po_fft)) + 2, label='FFT (+2)', alpha=1)
ax2.legend()
ax2.set_xlabel('$n$')
ax2.set_ylabel('Amplituda + premik')
ax2.set_title('Avtokorelacijska funkcija za testni signal')


ax3.plot(f, t1, label='python', alpha=1)
ax3.plot(f, t2, label='numpy', alpha=1)
ax3.plot(f, t3, label='FFT', alpha=0.5)
ax3.legend()
ax3.set_yscale('log')
ax3.set_xlabel('$ n $')
ax3.set_ylabel('Čas izračuna [s]')
ax3.set_title('Časovna zahtevnost')


plt.tight_layout()
plt.savefig('05/graphs/testni.pdf', dpi=512)
plt.close(fig)


#### sovice fft ####


fig, axes = plt.subplots(6, 1, figsize=(10, 8), sharex=True)
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

ac_fft, ac_freq = fft(data[1], sampling[1])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax1.plot(ac_freq[mask], ac_fft[mask], label=f'{names[1]}')
ax1.legend()
ax1.set_ylabel('Amplituda')

ac_fft, ac_freq = fft(data[5], sampling[5])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax2.plot(ac_freq[mask], ac_fft[mask], label=f'{names[5]}')
ax2.legend()
ax2.set_ylabel('Amplituda')

ac_fft, ac_freq = fft(data[0], sampling[0])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax3.plot(ac_freq[mask], ac_fft[mask], label=f'{names[0]}')
ax3.legend()
ax3.set_ylabel('Amplituda')

ac_fft, ac_freq = fft(data[2], sampling[2])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax4.plot(ac_freq[mask], ac_fft[mask], label=f'{names[2]}')
ax4.legend()
ax4.set_ylabel('Amplituda')

ac_fft, ac_freq = fft(data[3], sampling[3])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax5.plot(ac_freq[mask], ac_fft[mask], label=f'{names[3]}')
ax5.legend()
ax5.set_ylabel('Amplituda')

ac_fft, ac_freq = fft(data[4], sampling[4])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax6.plot(ac_freq[mask], ac_fft[mask], label=f'{names[4]}')
ax6.legend()
ax6.set_xlabel('Frekvenca [Hz]')
ax6.set_ylabel('Amplituda')


plt.savefig('05/graphs/fft.pdf', dpi=512)
plt.close(fig)


#### sovice auto ####

fig, axes = plt.subplots(6, 1, figsize=(10, 8), sharex=True)

for i, ax in enumerate(axes):
    sr = sampling[i]
    ac_fft, ns_fft = auto_corr_fft_rel(data[i])
    ac_fft = np.abs(ac_fft)
    ax.plot(ns_fft / sr, ac_fft, label=names[i])
    ax.legend()
    ax.set_ylabel('Amplituda')
    if i ==5:
        ax.set_xlabel('Čas zamika [s]')

plt.savefig('05/graphs/ac.pdf', dpi=512)
plt.close(fig)


#### sovice ac in fft ####



fig, axes = plt.subplots(6, 1, figsize=(10, 8), sharex=True)

for i, j in zip([1,5,0,2,3,4], [0,1,2,3,4,5]):
    sr = sampling[i]
    ac_fft, ns_fft = auto_corr_fft_rel(data[i])
    ac_fft, ac_freq = fft(ac_fft, sample_rate)
    ac_fft = np.abs(ac_fft)
    mask = (ac_freq > 250) & (ac_freq < 500)
    axes[j].plot(ac_freq[mask], ac_fft[mask], label=names[i])
    axes[j].set_ylabel('Amplituda')
    if j == 5:
        axes[j].set_xlabel('Frekvenca [Hz]')

plt.savefig('05/graphs/acfft.pdf', dpi=512)
plt.close(fig)



fig, axes = plt.subplots(6, 1, figsize=(10, 8), sharex=True)

for i, j in zip([1,5,0,2,3,4], [0,1,2,3,4,5]):
    sr = sampling[i]
    ac_fft, ns_fft = auto_corr_fft_rel(data[i])
    ac_fft, ns_fft = auto_corr_fft_rel(ac_fft)
    ac_fft, ac_freq = fft(ac_fft, sample_rate)
    ac_fft = np.abs(ac_fft)
    mask = (ac_freq > 250) & (ac_freq < 500)
    axes[j].plot(ac_freq[mask], ac_fft[mask], label=names[i])
    axes[j].set_ylabel('Amplituda')
    if j == 5:
        axes[j].set_xlabel('Frekvenca [Hz]')

plt.savefig('05/graphs/acacfft.pdf', dpi=512)
plt.close(fig)

#### avtokorelacija rng ####


rand = 2*np.random.random(2**25) - 1


corr = correlate(rand, rand, mode='full')
corr = corr[corr.size//2:]
corr /= corr[0]
plt.plot(corr)
plt.title('Avtokorelacija NUMPY naključnega signala')
plt.grid()
plt.savefig('05/graphs/nump.pdf', dpi=512)
plt.clf()


def randu(n, seed=1):
    a = 65539
    m = 2**31
    x = seed
    out = []
    for _ in range(n):
        x = (a * x) % m
        out.append(x / m)
    return np.array(out, float)

rand = 2*randu(2**25)-1

corr = correlate(rand, rand, mode='full')
corr = corr[corr.size//2:]
corr /= corr[0]
plt.plot(corr)
plt.title('Avtokorelacija RANDU naključnega signala')
plt.grid()
plt.savefig('05/graphs/rand.pdf', dpi=512)
plt.clf()