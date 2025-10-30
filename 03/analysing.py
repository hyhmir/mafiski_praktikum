import rust as rs
import numpy as np
import timeit
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.special import factorial
import scienceplots

plt.style.use('science')

# t1 = []
# t2 = []
# t4 = []

# for n in range(1, 100):
#     t1.append(timeit.timeit(lambda: rs.q_4(n, '1'), number=100)) # type: ignore
#     t2.append(timeit.timeit(lambda: rs.q_4(n, '2'), number=100)) # type: ignore
#     t4.append(timeit.timeit(lambda: rs.q_4(n, '4'), number=100)) # type: ignore


# plt.plot(t1, label='$\\left[q\\right]^4$')
# plt.plot(t2, label='$\\left[q^2\\right]^2$')
# plt.plot(t4, label='$\\left[q^4\\right]$')
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
# plt.yscale('log')
# plt.grid()
# plt.xlabel('Velikost matrike')
# plt.ylabel('$t$ [$s$]')
# plt.title('Hitrost sestavljanja Hamiltonke')
# plt.savefig('03/graphs/q4_speed.pdf', dpi=512)
# plt.clf()

# mat1 = rs.q_4(10, '1')
# mat2 = rs.q_4(10, '2')
# mat4 = rs.q_4(10, '4')

# err_1 = rs.matsum(mat4, rs.mul(-1., mat1))
# err_2 = rs.matsum(mat4, rs.mul(-1., mat2))

# plt.imshow(err_1, cmap='viridis')
# plt.colorbar()
# plt.title('Napaka 10x10 matrike')# prikazana kot t.i.: \\textit{heatmap}$')
# plt.savefig('03/graphs/err1.pdf', dpi=512)
# plt.clf()


# plt.imshow(err_2, cmap='viridis')
# plt.colorbar()
# plt.title('Napaka 10x10 matrike')# prikazana kot t.i.: \\textit{heatmap}$')
# plt.savefig('03/graphs/err2.pdf', dpi=512)
# plt.clf()


def hamiltonka(n:int, lam:float):
    return rs.matsum(rs.harm(n), rs.mul(lam, rs.q_4(n, '4')))


# t_jac = []
# t_qr_hous = []
# t_lanczos = []
# t_nump = []
# t_numph = []

# for i in range(10, 100):
#     t_jac.append(timeit.timeit(lambda: rs.jacobi_eigen(hamiltonka(i, 0.5), 10000, 1e-9), number=10))
#     if i < 30:
#         t_qr_hous.append(timeit.timeit(lambda: rs.qr_eigen(hamiltonka(i, 0.5), 10000, 1e-8), number=10))
#     t_lanczos.append(timeit.timeit(lambda: rs.lanczos_smallest(hamiltonka(i, 0.5), i, None, None), number=10))
#     # t_nump.append(timeit.timeit(lambda: np.linalg.eigh(hamiltonka(i, 0.5)), number=10))
#     t_numph.append(timeit.timeit(lambda: np.linalg.eig(hamiltonka(i, 0.5)), number=10))


# plt.plot(t_jac, label='Jacobi')
# plt.plot(t_qr_hous, label='QR')
# plt.plot(t_lanczos, label='Lanzcos')
# # plt.plot(t_nump, label='nump')
# plt.plot(t_numph, label='numpy.linalg.eigh')
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
# plt.yscale('log')
# plt.grid()
# plt.xlabel('Velikost matrike')
# plt.ylabel('$t$ [$s$]')
# plt.title('Hitrost različnih algoritmov')
# plt.savefig('03/graphs/algos.pdf')
# plt.clf()
# plt.show()

# vals, vecs = rs.lanczos_smallest(hamiltonka(200, 0.5), 200, None, None)
# val, vec = np.linalg.eigh(hamiltonka(200, 0.5))

# idx = np.argsort(val)
# val = val[idx]
# vec = vec[:, idx]

# val_err = []
# vali_err = []
# # print(np.array(hamiltonka(200, 0.5)).shape)
# # print(np.array(vecs[0]).shape)
# vecs = rs.transposey(vecs)
# vecs = np.array(vecs)
# vec = vec.T #/ np.max(vec)
# # # vecs = np.array(vecs) / np.max(vecs)
# valn_err = []

# # print(vecs[0])
# for i in range(200):
#     val_err.append(np.abs(vals[i] - val[i]))
#     vali_err.append(np.linalg.norm(np.array(hamiltonka(200, 0.5)) @ np.array(vecs[i]) - vals[i] * np.array(vecs[i])))
#     valn_err.append(np.linalg.norm(np.array(hamiltonka(200, 0.5)) @ vecs[i] - val[i] * vecs[i]))

# for i in range(len(vecs[0])):
#     vec[i] = vec[i] / vec[i][np.argmax(np.abs(vec[i]))]
#     vecs[i] = vecs[i] / vecs[i][np.argmax(np.abs(vecs[i]))]

# vecs_err = rs.matsum(vec, rs.mul(-1., vecs))

# plt.plot(val_err)
# plt.show()

# plt.plot(vali_err)
# plt.grid()
# plt.xlabel('Zaporedna lastna vrednost')
# plt.ylabel('$|Hv_i - \\lambda_iv_i|_2$')
# plt.title('Napaka Lanczosovega algoritma')
# plt.savefig('03/graphs/errlanz.pdf', dpi=512)
# plt.clf()

# plt.plot(valn_err)
# plt.grid()
# plt.xlabel('Zaporedna lastna vrednost')
# plt.ylabel('$|Hv_i - \\lambda_iv_i|_2$')
# plt.title('Napaka numpy.linalg.eigh algoritma')
# plt.savefig('03/graphs/errnumpy.pdf', dpi=512)
# plt.clf()


# plt.imshow(vecs_err[:120][:120])
# plt.colorbar()
# plt.show()

# plt.imshow(vecs)
# plt.colorbar()
# plt.show()


# plt.imshow(vec)
# plt.colorbar()
# plt.show()


# vals, vecs = np.linalg.eigh(hamiltonka(4000, 1))

# plt.plot(vals)
# plt.yscale('log')
# plt.show()

# plt.imshow(vecs)
# plt.show()


####### potrebn n ######

# eigs = [[] for _ in range (1, 10)]
# inds = [i for i in range(1, 61)]
# print(inds)
# for i in range(1,61):
#     vals, vecs = np.linalg.eigh(hamiltonka(i, 1))
#     vals = np.sort(vals)
#     for j in range(9):
#         try:
#             eigs[j].append(vals[j])
#         except:
#             eigs[j].append(0)

# print(np.array(eigs).shape)
# print(np.array(inds).shape)

# for i in [0, 1, 3, 8]:
#     plt.plot(inds[i:], eigs[i][i:], label=f'$E_{i}$')

# plt.yscale('log')
# plt.legend(loc='upper right',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
# plt.grid()
# plt.title('Konvergenca lastnih energij')
# plt.ylabel('Izračunana lastna energija')
# plt.xlabel('Velikost matrike')
# plt.savefig('03/graphs/energ_konv.pdf', dpi=512)
# plt.clf()

###### e od l ######

# for l in [0, 0.01, 0.1, 0.5, 1]:
#     ham = hamiltonka(2000, l)
#     vals, vecs = np.linalg.eigh(ham)
#     vals = np.sort(vals)
#     plt.plot(vals[:100], label=f'$\\lambda = {l}$')

# plt.grid()
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
# plt.title('Lastne energije')
# plt.ylabel('Lastna energija')
# plt.xlabel('Zaporedna lasntna energija')
# plt.savefig('03/graphs/energ_.pdf', dpi=512)
# plt.clf()
#### n od i #####
# tol = 1e-6

# for l in [0, 0.01, 0.1, 0.5, 1]:
#     vals, vecs = np.linalg.eigh(hamiltonka(3000, l))
#     vals = np.sort(vals)
#     di = {0:[], 0.001:[], 0.01:[], 0.1:[], 0.5:[], 1:[]}
#     for i in tqdm(range(80)):
#         right = 2000
#         left = i
#         while right - left > 2:
#             val, vec = np.linalg.eigh(hamiltonka((right + left)//2 + 1, l))
#             val = np.sort(val)
#             if np.abs(val[i] - vals[i]) < tol:
#                 right = (right + left)//2
#             else:
#                 left = (right + left)//2 + 1
        

#             # val, vec = np.linalg.eigh(hamiltonka(j+1, l))
#             # val = np.sort(val)
#             # # print(val.shape)
#             # if np.abs(val[i] - vals[i]) < tol:
#             #     di[l].append(i+1)
#             #     break
#             # j += 1
#         di[l].append(right)
    
#     plt.plot(di[l], label=f'$\\lambda = {l}$')

# plt.grid()
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
# plt.title('Konstantna napaka')
# plt.ylabel('Potrebna velikost matrike')
# plt.xlabel('Zaporedna lasntna energija')
# plt.savefig('03/graphs/nodi.pdf', dpi=512)
# plt.clf()

##### risanje funkcij ####
def normalize(vector):
    normalized_vector = vector / np.linalg.norm(vector)
    
    dominant_index = np.argmax(np.abs(normalized_vector))

    if normalized_vector[dominant_index] < 0:
        normalized_vector = -normalized_vector
    
    return normalized_vector

# for i in range(4):
#     for l in [0, 0.1, 0.5, 1]:
#         vals, vecs = np.linalg.eigh(hamiltonka(100, l))
#         idx = np.argsort(vals)
#         vals = vals[idx]
#         vecs = vecs[:, idx].T
#         if i==3 and l in [0.5, 1]:
#             hermite = np.polynomial.hermite.Hermite(-normalize(vecs[i])*[(1 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))) for n in range(1, len(vecs[0]) + 1)])
#         else:
#             hermite = np.polynomial.hermite.Hermite(normalize(vecs[i])*[(1 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))) for n in range(1, len(vecs[0]) + 1)])
#         def wave(x):
#             return np.exp(-x**2 / 2) * hermite(x)

#         x = np.linspace(-6, 6, 1000)
#         plt.plot(x, wave(x), label=f'$\\lambda = {l}$')
#     plt.grid()
#     plt.legend(loc='best',
#                 frameon=True,
#                 framealpha=0.9,
#                 facecolor='white',
#                 edgecolor='gray'
#             )
#     plt.title(f'Lastne funkcije {i}-ega \n vzbujenega stanja')
#     plt.ylabel('$\\psi (q)$')
#     plt.xlabel('$q$')
#     plt.savefig(f'03/graphs/eig{str(i).strip('.')}.pdf', dpi=512)
#     plt.clf()


# for i in range(5):
#     vals, vecs = np.linalg.eigh(hamiltonka(100, 1))
#     idx = np.argsort(vals)
#     vals = vals[idx]
#     vecs = vecs[:, idx].T
#     hermite = np.polynomial.hermite.Hermite(normalize(vecs[i])*[(1 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))) for n in range(1, len(vecs[0]) + 1)])
#     def wave(x):
#         return np.exp(-x**2 / 2) * hermite(x)
#     x = np.linspace(-2.5, 2.5, 1000)
#     def prob(x):
#         return wave(x)**2
#     def v(x, l):
#         return 0.5*x**2 + l * x**4



#     plt.plot(x, wave(x) + vals[i], label=f'$E_{i} = {round(vals[i], 1)}$')
#     plt.text(-2.4, vals[i]+0.5, f'$E_{i} = {round(vals[i], 1)}$')
#     plt.fill_between(x, vals[i], wave(x) + vals[i], alpha=0.2)

# plt.ylim((-0.1, 12.5))
# plt.plot(x, v(x, 1), c='black')
# plt.grid()
# # plt.legend(loc='best',
# #             frameon=True,
# #             framealpha=0.9,
# #             facecolor='white',
# #             edgecolor='gray'
# #         )
# plt.title(f'Lastne funkcije')
# plt.ylabel('$E$ [$\\hbar \\omega$]')
# plt.xlabel('$q$')
# plt.savefig('03/graphs/pretty.pdf', dpi=512)
# plt.clf()


### dodatna ####

for i in [0, 2, 4, 6, 8]:
    vals, vecs = np.linalg.eigh(rs.dodatna(100))
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx].T
    hermite = np.polynomial.hermite.Hermite(normalize(vecs[i])*[(1 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))) for n in range(1, len(vecs[0]) + 1)])
    def wave(x):
        return np.exp(-x**2 / 2) * hermite(x)
    x = np.linspace(-4.6, 4.6, 1000)
    def prob(x):
        return wave(x)**2
    def v(x, l):
        return -2*x**2 + 0.1 * x**4



    plt.plot(x, wave(x) + vals[i], label=f'$E_{i} = {round(vals[i], 1)}$')
    plt.text(-3.2, vals[i]+0.5, f'$E_{i} = {round(vals[i], 1)}$')
    plt.fill_between(x, vals[i], wave(x) + vals[i], alpha=0.2)

# plt.plot(x, v(x, 1), c='black')
# plt.show()
# plt.ylim((-0.1, 12.5))
plt.plot(x, v(x, 1), c='black')
plt.grid()
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
plt.title(f'Sode lastne funkcije')
plt.ylabel('$E$ [$\\hbar \\omega$]')
plt.xlabel('$q$')
plt.savefig('03/graphs/dodatnasod.pdf', dpi=512)
plt.clf()


for i in [1, 3, 5, 7, 9]:
    vals, vecs = np.linalg.eigh(rs.dodatna(100))
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx].T
    hermite = np.polynomial.hermite.Hermite(normalize(vecs[i])*[(1 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))) for n in range(1, len(vecs[0]) + 1)])
    def wave(x):
        return np.exp(-x**2 / 2) * hermite(x)
    x = np.linspace(-4.6, 4.6, 1000)
    def prob(x):
        return wave(x)**2
    def v(x, l):
        return -2*x**2 + 0.1 * x**4



    plt.plot(x, wave(x) + vals[i], label=f'$E_{i} = {round(vals[i], 1)}$')
    plt.text(-3.2, vals[i]+0.5, f'$E_{i} = {round(vals[i], 1)}$')
    plt.fill_between(x, vals[i], wave(x) + vals[i], alpha=0.2)

# plt.plot(x, v(x, 1), c='black')
# plt.show()
# plt.ylim((-0.1, 12.5))
plt.plot(x, v(x, 1), c='black')
plt.grid()
# plt.legend(loc='best',
#             frameon=True,
#             framealpha=0.9,
#             facecolor='white',
#             edgecolor='gray'
#         )
plt.title(f'Lihe lastne funkcije')
plt.ylabel('$E$ [$\\hbar \\omega$]')
plt.xlabel('$q$')
plt.savefig('03/graphs/dodatnalih.pdf', dpi=512)
plt.clf()
