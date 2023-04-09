import Utils.utils
from RLWE.ring_element import RingElement
import random
from numpy.polynomial import Polynomial
import numpy as np


def main():
    log_n = 5
    n = 2 ** log_n
    q = 1001
    coeffs = []
    for i in range(n // 2):
        random_bit = random.randint(0, q)
        coeffs.append(random_bit)
    a = Polynomial(coeffs)

    hamming_weight = 64
    arr = [0] * (n // 2)
    for i in range(hamming_weight):
        random_index = random.randrange(n // 2)
        if random.randrange(2) == 1:
            arr[random_index] = 1
        else:
            arr[random_index] = -1
    s = Polynomial(arr)

    print('a=', a.coef[:10])
    print('s=', s.coef[:10])
    print('|a|=', Utils.utils.canonical_norm(a, n))
    print('|s|=', Utils.utils.canonical_norm(s, n))
    phim = Utils.utils.build_cyclotomic_poly(n)
    print('|a*s|=', Utils.utils.canonical_norm(a*s, n))
    print('a*s=', (a*s).coef[:10])
    poly_modulo_phi_m = Utils.utils.modulo_polynomial(a*s, phim)
    print('a*s/phim=', poly_modulo_phi_m.coef[:10])



if __name__ == '__main__':
    main()
