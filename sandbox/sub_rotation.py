from math import *
import cmath
import numpy as np
import typing
import random
from numpy.polynomial import Polynomial
from RLWE.ring_element import RingElement
from Utils import utils

random.seed(0)
n = 16
q = 10000019
k = 4
encode_mat = [[0] * (n) for i in range(n)]
decode_mat = [[0] * (n) for i in range(n)]
five_powers_order = n // 2
scale = 2**9

def test_automorphism():
    ####### test that m(\zeta^k)=m'(\zeta) for m'(X)=m(X^k).
    zeta = cmath.exp((2j * pi) / (2 * n))
    U, U_conj = utils.generate_canonical_power_of_five(2 * n)
    b = [random.randint(1, 100) for _ in range(8)]
    print(b)
    a = utils.encode(U, U_conj, 2 * n, b)
    a_r = [int(x.real * scale) for x in a]
    a_elm = RingElement(Polynomial(a_r), 2 * n, q)
    print('decode(a_elm/scale)=', utils.decode(U, U_conj, [x / scale for x in a_elm.poly.coef]))
    k = 2
    power = (5 ** k)
    eval_a_zeta_power_five = np.polyval(a_elm.poly.coef[::-1], zeta ** power)
    a_k = a_elm.automorphism(power)
    print('decode(a_k/scale)=', utils.decode(U, U_conj, [x / scale for x in a_k.poly.coef]))
    eval_a_k_zeta = np.polyval(a_k.poly.coef[::-1], zeta)
    print('eval_a_zeta_power_five: ', eval_a_zeta_power_five)
    print('eval_a_k_zeta: ', eval_a_k_zeta)
    assert abs(eval_a_zeta_power_five - eval_a_k_zeta) < 0.001


def main():
    # test_automorphism()
    k = 2
    i = 1
    print('params: k=',k)
    print('params: i=',i)
    five_powers_order = n // 2  # size of max rotation
    zeta = cmath.exp((2j * pi) / (2 * n))
    U, U_conj = utils.generate_canonical_power_of_five_cosets(2 * n, k)
    # b = [random.randint(1, 100) for _ in range(8)]
    slots = [complex(round(random.uniform(-10, 10), 2), round(random.uniform(-10, 10), 2)) for _ in
                              range(8)]

    print('slots=              ',slots)
    a = utils.encode(U, U_conj, 2 * n, slots)
    a_r = [int(x.real * scale) for x in a]
    a_elm = RingElement(Polynomial(a_r), 2 * n, q)
    decode_a_elm = utils.decode(U, U_conj, [x / scale for x in a_elm.poly.coef])
    print('decode(a_elm/scale)=', [complex(round(x.real, 2), round(x.imag, 2)) for x in decode_a_elm][:n//2])

    power = (5 ** ((five_powers_order//k)*i) % (2*n))
    eval_a_zeta_power_five = np.polyval(a_elm.poly.coef[::-1], zeta ** power)
    a_k = a_elm.automorphism(power)
    decode_a_k = utils.decode(U, U_conj, [x / scale for x in a_k.poly.coef])
    print('decode(a_k/scale)=  ', [complex(round(x.real, 2), round(x.imag, 2)) for x in decode_a_k][:n//2])
    eval_a_k_zeta = np.polyval(a_k.poly.coef[::-1], zeta)
    print('eval_a_zeta_power_five: ', eval_a_zeta_power_five)
    print('eval_a_k_zeta: ', eval_a_k_zeta)
    assert abs(eval_a_zeta_power_five - eval_a_k_zeta) < 0.001


if __name__ == '__main__':
    main()