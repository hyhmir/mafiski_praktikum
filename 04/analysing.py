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


#### dft nekaj osnovnih primerov #####

n = 200 # Number of data points
T = 100. # Sampling period
dt = T/n
tmin=0.
tmax=dt*n
print("sampling freq:",1./dt)
nuc=0.5/dt
print("critical freq:",nuc)

t=np.linspace(tmin,tmax,n,endpoint=False)
#sine
nu = 0.5
signal = np.sin(2*np.pi*t*nu).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))
plt.show()

# cosine
nu = 0.5
signal = np.cos(2*np.pi*t*nu).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))

plt.show()


# combo
nu = 0.5
signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))

plt.show()


# Non periodic
t=np.linspace(tmin,tmax,n)#,endpoint=False)

nu = 0.5
signal = (np.cos(2*np.pi*t*nu) + np.sin(2*np.pi*t*nu/2)).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n)#,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))

plt.show()


# previsoka frekvenca
t=np.linspace(tmin,tmax,n, endpoint=False)

nu = 2.3
signal = np.cos(2*np.pi*t*nu).tolist()

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n)#,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))

plt.show()

# gaussovka neperiodicna
t=np.linspace(tmin,tmax,n,endpoint=False)

signal = rs.gauss(0.1*n, n//2, False)

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n, endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))
plt.show()


# gaussovka periodicna
t=np.linspace(tmin,tmax,n,endpoint=False)

signal = rs.gauss(0.1*n, n//2, True)

dft = rs.dft(signal)
Hk=np.roll(dft,int(n/2))/n


plt.plot(t, signal)
plt.show()

nus = np.linspace(-nuc, nuc, n,endpoint=False)
plt.plot(nus, np.real(Hk))
plt.plot(nus, np.imag(Hk))
plt.show()



#### casovna zahtevnost ####

# times = np.zeros((100, 10))

# for i, n in enumerate(np.linspace(10, 1000, 100)):
#     gauss = rs.gauss(0.1*n, int(n), True)
#     for j in range(10):
#         times[i][j] = timeit.timeit(lambda: rs.dft(gauss), number=1)

# mean_times = times.mean(axis=1)
# std_times = times.std(axis=1)

# plt.errorbar(np.linspace(10, 1000, 100), mean_times, std_times)
# plt.show()


