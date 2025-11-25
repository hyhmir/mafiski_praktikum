from rust import euler, analyt, heun, rk2a, rku4, rk45, rkf, pc4 # type: ignore
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

#### testiranje ######

t = np.linspace(0, 10, 10)

analytical = analyt(t, -100)

eul = euler(t, -100)

heu = heun(t, -100)

ruku2a = rk2a(t, -100)

ruku4 = rku4(t, -100)

ruku45, err = rk45(t, -100)

rukuf, tr = rkf(0, 80, 21, tol=1e-13)
analitr = analyt(tr, 21)

pc = pc4(t, -100)

# plt.plot(t, analytical)
# plt.plot(t, pc)
# # plt.plot(t, heu)
# plt.show()


# plt.plot(t, np.abs(np.array(analytical) - np.array(pc)))
# plt.yscale('log')
# plt.show()


plt.plot(tr, analitr)
plt.plot(tr, rukuf)
# plt.plot(t, heu)
plt.show()


plt.plot(tr, np.abs(np.array(analitr) - np.array(rukuf)))
plt.yscale('log')
plt.show()

cajta = np.linspace(0, 80, 1000)
analytical = analyt(cajta, 21)


#### euler #####

# for h in [0.1, 1, 10]:
#     t = np.linspace(0, 80, int(80/h))
#     eul = euler(t, 21)
#     plt.plot(t, eul, label=f'h = {h}')


# plt.plot(cajta, analytical, '-.',label='analytical')
# plt.grid()
# plt.legend()
# plt.title('Uspešnost Eulerjeve metode')
# plt.savefig('06/graphs/euly.pdf', dpi=512)


