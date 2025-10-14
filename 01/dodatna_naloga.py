import mpmath
from mpmath import mp
import matplotlib.pyplot as plt
import scienceplots
import numpy as np


mp.dps = 100


def a_zero(x):
    z = (3/8) * mp.pi * (4 * x - 1)
    return z**(2/3) * (1 + 5/48 * z**(-2) - 5/36 * z**(-4) + 77125/82944 * z**(-6) - 108056875/6967296 * z**(-8))


def b_zero(x):
    z = (3/8) * mp.pi * (4 * x - 3)
    return z**(2/3) * (1 + 5/48 * z**(-2) - 5/36 * z**(-4) + 77125/82944 * z**(-6) - 108056875/6967296 * z**(-8))


if __name__ == '__main__':
    n = list(range(1, 101))
    ai_zeros = mp.matrix(1, len(n))
    bi_zeros = mp.matrix(1, len(n))
    ai_zeros_ref = mp.matrix(1, len(n))
    bi_zeros_ref = mp.matrix(1, len(n))

    for i in n:
        ai_zeros[i - 1] = -a_zero(i)
        bi_zeros[i - 1] = -b_zero(i)
        ai_zeros_ref[i - 1] = mp.airyaizero(i)
        bi_zeros_ref[i - 1] = mp.airybizero(i)

    
    plt.style.use(['science', 'bright'])

    plt.plot(n, ai_zeros, ls='-', label='Ai')
    plt.plot(n, bi_zeros, ls='--', label='Bi')

    plt.grid()
    plt.legend(loc='best',
        frameon=True,        # turn on legend box
        framealpha=0.9,      # 0 = transparent, 1 = opaque
        facecolor='white',   # background color
        edgecolor='gray'     # border color
    )
    plt.title('Ničle funkcij $Ai$ in $Bi$')
    plt.xlabel('n - št. ničle')
    plt.ylabel('x')

    plt.savefig('01/graphs/nicle_draw.pdf', dpi=512)
    plt.clf()


    plt.plot(n, np.abs(ai_zeros - ai_zeros_ref), label='Ai', ls='-')
    plt.plot(n, np.abs(bi_zeros - bi_zeros_ref), label='Bi', ls='--')

    plt.grid()
    plt.legend(loc='best',
        frameon=True,        # turn on legend box
        framealpha=0.9,      # 0 = transparent, 1 = opaque
        facecolor='white',   # background color
        edgecolor='gray'     # border color
    )
    plt.title('Absolutna napaka ničel')
    plt.xlabel('n - št. ničle')
    plt.ylabel('x')
    plt.yscale('log')

    plt.savefig('01/graphs/nicle_err.pdf', dpi=512)
    plt.clf()

    print(mp.nstr(ai_zeros[:10], 6))
    print(mp.nstr(ai_zeros_ref[:10], 6))
    print(mp.nstr(bi_zeros[:10], 6))
    print(mp.nstr(bi_zeros_ref[:10], 6))

