import rust as rs
import numpy as np
import matplotlib.pyplot as plt
import timeit
import pandas as pd
import time

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


#### speeeeed #####

# loglengs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
# lengs = 2**loglengs
# print(lengs)

# numpt = []
# numperr = []

# rustt = []
# rusterr = []

# for len in lengs:
#     h = np.random.random(len)
#     timesr = []
#     timesn = []
#     for i in range(10):
#         timesn.append(timeit.timeit(lambda: np.fft.fft(h), number=1))
#         timesr.append(timeit.timeit(lambda: rs.fft(h, False), number=1))
#     numpt.append(np.mean(timesn))
#     numperr.append(np.std(timesn))
#     rustt.append(np.mean(timesr))
#     rusterr.append(np.std(timesr))

# plt.errorbar(lengs, numpt, numperr)
# plt.errorbar(lengs, rustt, rusterr)
# # plt.xscale('log')
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



# fig, axes = plt.subplots(3, 2, figsize=(10, 15))
# axes = axes.flatten()

# norm = 1

# for i, ax in enumerate(axes):
#     d = data[i] / norm
#     sample_rate = sampling[i]
#     pxx, freqs, bins, im = ax.specgram(d, Fs=sample_rate)
#     ax.set_title(f'Spektrogram {names[i]}.wav')
#     ax.set_xlabel('Čas [s]\n\n')
#     ax.set_ylabel('Frekvenca [Hz]')
#     fig.colorbar(im, ax=ax)

# plt.tight_layout()

# plt.show()
# plt.close(fig)


#### avtokorelacija ####

def fft(x, fs, N=0, shift=False):
    if N == 0:
        N = len(x)
    X = np.fft.fft(x, N)
    f = np.fft.fftfreq(N, 1/fs)
    if shift:
        X = np.fft.fftshift(X)
        f = np.fft.fftshift(f)
    return X, f


def ifft(X, fs, N=0):
    if N == 0:
        N = len(X)
    x = np.fft.ifft(X, N)
    t = np.fft.fftfreq(N, 1/fs)
    return x, t


def pad_data(data):
    N = len(data)
    pad = np.zeros(N)
    return np.hstack((data, pad))


def auto_corr(h):
    times = []
    t0 = time.time()
    N = h.shape[0]
    ns = np.arange(N)
    fis = []
    h_pad = pad_data(h)
    for n in ns:
        h_pad = pad_data(h)
        h_pad_roll = np.roll(h_pad, n)
        A = 1 / (N - n)
        out = A * np.dot(h_pad, h_pad_roll)
        fis.append(out)
        t1 = time.time()
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
ax1.set_title('Tesni primer - kompleksno sestavljen zašumljen signal')

ax2.plot(f, np.array(po_def), label='definicija python', alpha=1)
ax2.plot(f, np.real(np.array(po_def_np)) + 1, label='definicija numpy', alpha=1)
ax2.plot(f, np.real(np.array(po_fft)) + 2, label='FFT', alpha=1)
ax2.legend()
ax2.set_xlabel('$n$')
ax2.set_ylabel('Amplituda + premik')
ax2.set_title('Avtokorelacijska funkcija za testni signal')


ax3.plot(f, t1, label='definicija python', alpha=1)
ax3.plot(f, t2, label='definicija numpy', alpha=1)
ax3.plot(f, t3, label='FFT', alpha=0.5)
ax3.legend()
ax3.set_yscale('log')
ax3.set_xlabel('$ n $')
ax3.set_ylabel('Čas izračuna [s]')
ax3.set_title('Časovna zahtevnost')


plt.tight_layout()
plt.show()
plt.close(fig)


#### sovice ####




fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot([335]*10, np.linspace(-0.2, 6, 10), zs=0, zdir='z', color='black', alpha=0.75, label='335 Hz', linestyle='dashed')
ax.plot([375]*10, np.linspace(-0.2, 6, 10), zs=0, zdir='z', color='black', alpha=0.75, label='375 Hz', linestyle='dashdot')

yticks = np.arange(6)
ord = [1, 5, 0, 2, 3, 4]
data_s = [data[i] for i in ord]
sampling_s = [sampling[i] for i in ord]
names_s = [names[i] for i in ord]
for i, d in enumerate(data_s):
    i = len(data_s) - i - 1
    sample_rate = sampling_s[i]
    ac_fft, ac_freq = fft(data_s[i], sample_rate)
    ac_fft = np.abs(ac_fft)
    mask = (ac_freq > 0) & (ac_freq < 1000)
    ax.plot(ac_freq[mask], ac_fft[mask], zs=i, zdir='y', alpha=0.8, label=names_s[i])

ax.set_xlabel('Časovni zamik [s]')
ax.set_ylabel('Različni posnetki')
ax.set_zlabel('Amplituda')
ax.legend()
ax.zaxis.labelpad=-1.2 # <- change the value here

plt.xlim(0, 1000)
plt.title('Frekvenčni spekter signalov $h$')
plt.tight_layout()
plt.show()



fig, axes = plt.subplots(6, 1, figsize=(10, 8))
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

ac_fft, ac_freq = fft(data[1], sampling[1])
ac_fft = np.abs(ac_fft)
mask = (ac_freq > 0) & (ac_freq < 1000)
ax1.plot(ac_freq[mask], ac_fft[mask], label=f'{names[1]}')

plt.show()


