import rust as rs
import numpy as np
import timeit
import matplotlib.pyplot as plt

# t1 = []
# t2 = []
# t4 = []

# for n in range(1, 100):
#     t1.append(timeit.timeit(lambda: rs.q_4(n, '1'), number=100)) # type: ignore
#     t2.append(timeit.timeit(lambda: rs.q_4(n, '2'), number=100)) # type: ignore
#     t4.append(timeit.timeit(lambda: rs.q_4(n, '4'), number=100)) # type: ignore


# plt.plot(t1, label='t1')
# plt.plot(t2, label='t2')
# plt.plot(t4, label='t4')
# plt.legend()
# plt.yscale('log')
# plt.show()

# mat1 = rs.q_4(10, '1')
# mat2 = rs.q_4(10, '2')
# mat4 = rs.q_4(10, '4')

# err_1 = rs.matsum(mat4, rs.mul(-1., mat1))
# err_2 = rs.matsum(mat4, rs.mul(-1., mat2))

# plt.imshow(err_1, cmap='viridis')
# plt.colorbar()
# plt.show()


# plt.imshow(err_2, cmap='viridis')
# plt.colorbar()
# plt.show()


def hamiltonka(n:int, lam:float):
    return rs.matsum(rs.harm(n), rs.mul(lam, rs.q_4(n, '4')))


# t_jac = []
# t_qr_hous = []
# t_lanczos = []
# t_nump = []

# for i in range(10, 100):
#     t_jac.append(timeit.timeit(lambda: rs.jacobi_eigen(hamiltonka(i, 0.5), 10000, 1e-8), number=10))
#     if i < 30:
#         t_qr_hous.append(timeit.timeit(lambda: rs.qr_eigen(hamiltonka(i, 0.5), 10000, 1e-8), number=10))
#     t_lanczos.append(timeit.timeit(lambda: rs.lanczos_smallest(hamiltonka(i, 0.5), i, None, None), number=10))
#     t_nump.append(timeit.timeit(lambda: np.linalg.eig(hamiltonka(i, 0.5)), number=10))

# plt.plot(t_jac, label='jac')
# plt.plot(t_qr_hous, label='qr')
# plt.plot(t_lanczos, label='lanzcos')
# plt.plot(t_nump, label='nump')
# plt.legend()
# plt.yscale('log')
# plt.show()

vals, vecs = rs.lanczos_smallest(hamiltonka(200, 0.5), 200, None, None)
val, vec = np.linalg.eig(hamiltonka(200, 0.5))

idx = np.argsort(val)
val = val[idx]
vec = vec[:, idx]

val_err = []
for i in range(200):
    val_err.append(np.abs(vals[i] - val[i]))


vecs = rs.transposey(vecs)
vecs = np.array(vecs)
vec = vec.T #/ np.max(vec)
# vecs = np.array(vecs) / np.max(vecs)

# print(vecs[0])

for i in range(len(vecs[0])):
    vec[i] = vec[i] / vec[i][np.argmax(np.abs(vec[i]))]
    vecs[i] = vecs[i] / vecs[i][np.argmax(np.abs(vecs[i]))]

vecs_err = rs.matsum(vec, rs.mul(-1., vecs))

plt.plot(val_err)
plt.show()

plt.imshow(vecs_err[:120][:120])
plt.colorbar()
plt.show()

plt.imshow(vecs)
plt.colorbar()
plt.show()


plt.imshow(vec)
plt.colorbar()
plt.show()

