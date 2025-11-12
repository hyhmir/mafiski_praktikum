import rust as rs
import numpy as np
import matplotlib.pyplot as plt
import timeit

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
ax1.set_title("Independent x-axis (top)")

# Bottom three subplots (shared x-axis)
ax2 = fig.add_subplot(gs[1, 0], sharex=None)  # first of the shared group
ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax2)

ax1.plot(t, signal,color='g')
ax2.plot(nus, np.real(Hk),color='b')
ax2.set_ylabel(r'$Re[H_k]$', size = 'x-large')
ax3.plot(nus, np.imag(Hk),color='r')
ax3.set_ylabel(r'$Im[H_k]$', size = 'x-large')
ax4.plot(nus, np.absolute(Hk)**2,color='y')
ax4.set_ylabel(r'$\vert H_k \vert ^2$', size = 'x-large')
ax4.set_xlabel(r'$\nu$', size = 'x-large')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)

fig.suptitle('Plot')
plt.show()


# # cosine
# nu = 0.5
# signal = np.cos(2*np.pi*t*nu).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n,endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))

# plt.show()


# # combo
# nu = 0.5
# signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n,endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))

# plt.show()


# # Non periodic
# t=np.linspace(tmin,tmax,n)#,endpoint=False)

# nu = 0.5
# signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n)#,endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))

# plt.show()


# # previsoka frekvenca - potujitev
# t=np.linspace(tmin,tmax,n, endpoint=False)

# nu = 2.3
# signal = np.cos(2*np.pi*t*nu).tolist()

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n, endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))

# plt.show()

# # gaussovka neperiodicna
# t=np.linspace(tmin,tmax,n,endpoint=False)

# signal = rs.gauss(0.1*n, n//2, False)

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n, endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))
# plt.show()


# # gaussovka periodicna
# t=np.linspace(tmin,tmax,n,endpoint=False)

# signal = rs.gauss(0.1*n, n//2, True)

# dft = rs.dft(signal)
# Hk=np.roll(dft,int(n/2))/n


# plt.plot(t, signal)
# plt.show()

# nus = np.linspace(-nuc, nuc, n,endpoint=False)
# plt.plot(nus, np.real(Hk))
# plt.plot(nus, np.imag(Hk))
# plt.show()



#### casovna zahtevnost ####

# times = np.zeros((10, 10))

# for i, n in enumerate(np.linspace(10, 1000, 10)): # vzame veliko casa
#     gauss = rs.gauss(0.1*n, int(n), True)
#     for j in range(10):
#         times[i][j] = timeit.timeit(lambda: rs.dft(gauss), number=1)

# mean_times = times.mean(axis=1)
# std_times = times.std(axis=1)

# plt.errorbar(np.linspace(10, 1000, 10), mean_times, std_times)
# plt.show()

# fitaj


#### numericna natancnost ####

# signal = np.array(rs.gauss(1000, 100, False))

# plt.plot(np.abs(signal - np.array(rs.idft(rs.dft(signal)))))
# plt.plot(np.abs(signal - np.fft.ifft(np.fft.fft(signal))))
# plt.yscale('log')
# plt.show()


#### bach ##### - uporabljamo numpy fft saj drugace traja 10 minut

# signal = np.loadtxt('04/samples/bach882.txt')
# plt.plot(np.abs(np.fft.fft(signal) - np.array(rs.dft(signal))))
# plt.show()

# signal = np.loadtxt('04/samples/bach44100.txt')
# plt.plot(np.roll(np.fft.fft(signal), len(signal)//2))
# plt.show() # enako za ostale

