import rust as rs
import numpy as np
import matplotlib.pyplot as plt
# import timeit
import time
import scienceplots
from scipy.optimize import curve_fit

plt.style.use('science')

#### testiranje ####

# a = [1,2,3,4,5,6,7,8,9,10]

# print(rs.dft(rs.gauss(1, 10, True)))
# print(rs.gauss(1, 10, True))

# plt.plot(rs.gauss(100, 1000, True))
# plt.show()
# plt.plot([z.real for z in rs.dft(rs.gauss(100, 1000, True))])
# plt.show()
# plt.plot(np.array(rs.gauss(100, 1000, True)) - np.array(rs.idft(rs.dft(rs.gauss(100, 1000, True)))))
# plt.show()

# n = 200 # Number of data points
# T = 100. # Sampling period
# dt = T/n
# tmin=0.
# tmax=dt*n
# print("sampling freq:",1./dt)
# nuc=0.5/dt
# print("critical freq:",nuc)
# t=np.linspace(tmin,tmax,n, endpoint=False)

# nu = 2.3
# signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*0.5) + np.cos(2*np.pi*0.9*t) + 0.2*(0.5 - np.random.random(len(t)))).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n, endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))

# plt.show()

# signal = np.real(rs.filter(signal, 2*nuc*dt)).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n, endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))
# plt.show()


#### dft nekaj osnovnih primerov #####

n = 100 # Number of data points
T = 100. # Sampling period
dt = T/n
tmin=0.
tmax=dt*n
print("sampling freq:",1./dt)
nuc=0.5/dt
print("critical freq:",nuc)

t=np.linspace(tmin,tmax,n,endpoint=False)
t_rep = np.linspace(tmin,tmax,100*n,endpoint=False)
# #sine
nu = nuc*0.2
signal = np.sin(2*np.pi*t*nu)

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n

nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/sin.pdf', dpi=512)
plt.close(fig)


# cosine
nu = nuc*0.2
signal = np.cos(2*np.pi*t*nu).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n

nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/cos.pdf', dpi=512)
plt.close(fig)


# combo
nu = 0.2*nuc
signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n

nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/combo.pdf', dpi=512)
plt.close(fig)


# # Non periodic
t=np.linspace(tmin,tmax,n)#,endpoint=False)

nu = 0.2*nuc
signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/non_periodic.pdf', dpi=512)
plt.close(fig)


# # previsoka frekvenca - potujitev
t=np.linspace(tmin,tmax,n, endpoint=False)

nu = 4.3*nuc
signal = np.cos(2*np.pi*t*nu).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
# ax1 = fig.add_subplot(gs[0, 0])
# ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[0, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[1, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[2, 0], sharex=ax2)

# ax1.plot(t, signal,color='g')
# ax1.set_ylabel(r'$h_k$', size = 'x-large')
# ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/potujitev_fur.pdf', dpi=512)
plt.close(fig)
# plt.show()
# plt.close(fig)


t_rep = np.linspace(tmin,tmax,100*n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
# ax1 = fig.add_subplot(gs[0, 0])
# ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[0, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[1, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[2, 0], sharex=ax2)

ax2.plot(t_rep[0:25*n], np.cos(2*np.pi*nu*t_rep[0:25*n]),color='b')
ax2.plot(t[0:n//4], signal[0:n//4],'.',color='black')
ax2.set_ylabel(r'$\sin (2\pi 4.3\nu_N t)$', size = 'x-large')
ax2.set_title('Signal $h_k$')
ax3.plot(t_rep[0:25*n], np.cos(2*np.pi*(nu-2*nuc)*t_rep[0:25*n]),color='r')
ax3.plot(t[0:n//4], signal[0:n//4],'.',color='black')
ax3.set_ylabel(r'$\sin (2\pi 2.3\nu_N t)$', size = 'x-large')
ax4.plot(t_rep[0:25*n], np.cos(2*np.pi*(nu-4*nuc)*t_rep[0:25*n]),color='y')
ax4.plot(t[0:n//4], signal[0:n//4],'.',color='black')
ax4.set_ylabel(r'$\sin (2\pi 0.3\nu_N t)$', size = 'x-large')
ax4.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.savefig('04/graphs/potujitev_sig.pdf', dpi=512)
plt.close(fig)
# plt.show()
# plt.close(fig)



# # gaussovka neperiodicna
t=np.linspace(tmin,tmax,n,endpoint=False)

signal = rs.gauss(0.1*n, n//2, False)

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n



nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/gauss_non_periodic.pdf', dpi=512)
plt.close(fig)


# # gaussovka periodicna
t=np.linspace(tmin,tmax,n,endpoint=False)

signal = rs.gauss(0.1*n, n//2, True)

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/gauss_periodic.pdf', dpi=512)
plt.close(fig)



#### casovna zahtevnost ####

times = np.zeros((10, 100))

for i, n in enumerate(np.linspace(10, 1000, 10)): # vzame veliko casa
    gauss = rs.gauss(0.1*n, int(n), True)
    for j in range(100):
        start = time.time()
        dft = rs.dft(gauss)
        stop = time.time()
        times[i][j] = stop-start

mean_times = times.mean(axis=1)
std_times = times.std(axis=1)

plt.errorbar(np.linspace(10, 1000, 10), mean_times, std_times, fmt='.', ms=2, label='Meritve')
# plt.show()

def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

popt, pcov = curve_fit(quadratic, np.linspace(10, 1000, 10), mean_times, sigma=std_times, absolute_sigma=True)

a, b, c = popt
a_err, b_err, c_err = np.sqrt(np.diag(pcov))

print(f"a = {a} ± {a_err}")
print(f"b = {b} ± {b_err}")
print(f"c = {c} ± {c_err}")

t_fit = np.linspace(10, 1000, 1000)
y_fit = quadratic(t_fit, a, b, c)

plt.plot(t_fit, y_fit, label='Kvadratični fit')
plt.grid()
plt.legend()
plt.title('Časovna zahtevnost')
plt.xlabel('Dolžina signala')
plt.ylabel('$t$ [$s$]')
plt.savefig('04/graphs/time.pdf', dpi=512)
plt.clf()


#### numericna natancnost ####

signal = np.array(rs.gauss(1000, 100, False))

plt.plot(np.abs(signal - np.array(rs.idft(rs.dft(signal))))[10:], label='Lastna implementacija')# pri fftju so ta prvi cist predobri in unicjo graf
plt.plot(np.abs(signal - np.fft.ifft(np.fft.fft(signal)))[10:], label='numpy.fft')
plt.yscale('log')
plt.grid()
plt.xlabel('št. indeksa')
plt.ylabel('Absolutna napaka')
plt.title(r'$\vert h_k - idft\left(dft\left(h_k\right)\right) \vert$')
plt.savefig('04/graphs/napaka.pdf', dpi=512)
plt.clf()
# plt.show()


#### bach ##### - uporabljamo numpy fft saj drugace traja 10 minut

t = 2.3

f, ax = plt.subplots(6,1,sharex=True)

signal = np.loadtxt('04/samples/bach44100.txt')
ax[0].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C0')
# ax[0].set_ylabel('44100 HZ')

signal = np.loadtxt('04/samples/bach11025.txt')
ax[1].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C1')
# ax[1].set_ylabel('11025 HZ')

signal = np.loadtxt('04/samples/bach5512.txt')
ax[2].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C2')
# ax[2].set_ylabel('5512 HZ')

signal = np.loadtxt('04/samples/bach2756.txt')
ax[3].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C3')
# ax[3].set_ylabel('2756 HZ')

signal = np.loadtxt('04/samples/bach1378.txt')
ax[4].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C4')
# ax[4].set_ylabel('1378 HZ')

signal = np.loadtxt('04/samples/bach882.txt')
ax[5].plot(np.linspace(-0.5*len(signal)/t, 0.5*len(signal)/t, len(signal)) ,np.abs(np.roll(np.fft.fft(signal), len(signal)//2)), color='C5')
# ax[5].set_ylabel('882 HZ')
ax[5].set_xlabel(r'$\nu$', size = 'x-large')

for a in ax:
    a.tick_params(labelleft=False)


plt.savefig('04/graphs/Bach.pdf', dpi=512)
plt.close(f)

#### filtering ####

n = 200 # Number of data points
T = 100. # Sampling period
dt = T/n
tmin=0.
tmax=dt*n
print("sampling freq:",1./dt)
nuc=0.5/dt
print("critical freq:",nuc)
t=np.linspace(tmin,tmax,n, endpoint=False)

nu = 0.3
signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*0.5) + np.cos(2*np.pi*0.9*t) + 0.5*(0.5 - np.random.random(len(t)))).tolist()
dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n

nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/unfiltered.pdf', dpi=512)
plt.close(fig)

signal = np.real(rs.filter(signal, 2*nuc*dt))

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


nus= np.linspace(-nuc,nuc,n,endpoint=False)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.4)

# Top subplot (independent x-axis)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title(r'Signal $h_k$')

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax1.set_ylabel(r'$h_k$', size = 'x-large')
ax1.set_xlabel(r'$t$ [$s$]', size = 'x-large', loc='right')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax2.set_title('Fouriereva transformiranka $H_k$')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$ [$Hz$]', size = 'x-large', loc='right')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
# plt.show()
plt.savefig('04/graphs/filtered.pdf', dpi=512)
plt.close(fig)
