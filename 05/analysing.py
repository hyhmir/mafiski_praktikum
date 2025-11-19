import rust as rs
import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0, 10, 1024)
h = np.sin(t)

plt.plot(rs.fft(h, True))
plt.plot(np.fft.ifft(h))
plt.show()

plt.plot(np.fft.ifft(np.fft.fft(h)) - rs.fft(rs.fft(h, False), True))
plt.show()

plt.plot(rs.fft(rs.fft(h, False), True))
plt.plot(h)
plt.show()

plt.plot(rs.fft(rs.fft(h, False), True) - h)
plt.show()