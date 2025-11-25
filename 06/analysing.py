from rust import euler, analyt, heun, rk2a, rku4, rk45, rkf, pc4 # type: ignore
import numpy as np
import matplotlib.pyplot as plt


#### testiranje ######

# t = np.linspace(0, 10, 10)

# analytical = analyt(t, -100)

# eul = euler(t, -100)

# heu = heun(t, -100)

# ruku2a = rk2a(t, -100)

# ruku4 = rku4(t, -100)

# ruku45, err = rk45(t, -100)

# rukuf, tr = rkf(0, 10, -100, tol=1e-13)
# analitr = analyt(tr, -100)

# pc = pc4(t, -100)

# plt.plot(t, analytical)
# plt.plot(t, pc)
# # plt.plot(t, heu)
# plt.show()


# plt.plot(t, np.abs(np.array(analytical) - np.array(pc)))
# plt.yscale('log')
# plt.show()


# plt.plot(tr, analitr)
# plt.plot(tr, rukuf)
# # plt.plot(t, heu)
# plt.show()


# plt.plot(tr, np.abs(np.array(analitr) - np.array(rukuf)))
# plt.yscale('log')
# plt.show()


#### euler #####
