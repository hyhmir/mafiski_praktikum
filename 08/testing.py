import rust as rs
import numpy as np
import matplotlib.pyplot as plt


y = rs.shooter(0, 0, 1, 0.5, np.linspace(0, 10, 1000), 1e-3, 1000, False, False)
# print(y[0:])
plt.plot(y[0:])
plt.show()