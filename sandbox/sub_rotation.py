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

def auto(x: int, l:int ) -> int : # X -> X^{5^l mod 2n}
    return x**(5**l % (2*n))


def gen_encoding_matrix():
    global encode_mat
    global decode_mat

    zeta = cmath.exp(2j * pi / (2 * n))
    roots = []
    # generate roots with special powers
    for j in range(0, five_powers_order // k):
        for i in range(0, k):
            coset_representative = zeta**(5**j)
            zeta_coset = auto(coset_representative, (five_powers_order / k) * i)
            roots.append(zeta_coset)

    #add the conjugated roots
    for l in range(0, n // 2):
        roots.append(np.conj(roots[n//2 - 1 - l]))
    assert len(roots) == n
    # build the matrix
    row_idx = 0
    for root in roots:
        for l in range(0, n):
            decode_mat[row_idx][l] = root ** l
        row_idx += 1
    assert len(decode_mat) == n and len(decode_mat[0]) == n
    encode_mat = np.linalg.inv(np.array(decode_mat))
    assert len(decode_mat) == n and len(decode_mat[0]) == n
    # for r in encode_mat:
    #     print(r)
    # print(np.dot(encode_mat,decode_mat))


def encode(slots)-> typing.List[complex]:
    assert len(slots) == n, "len(coeffs): " + str(len(slots))
    res = []
    for row in np.array(encode_mat):
        res.append(np.inner(row, slots))
        print('encode: ', np.inner(row, slots))
    res = [round(x.real * scale) for x in res] #scale
    return res


def decode(coeffs: typing.List[float]) -> typing.List[complex]:
    assert len(coeffs) == n
    coeffs = [x / scale for x in coeffs]
    res = []
    for row in np.array(decode_mat):
        res.append(np.inner(row, coeffs))
    # slots = decode_mat.dot(np.array(coeffs))
    return res


def gen_rnd_slots():
    numbers = []
    for _ in range(8):
        real_part = random.uniform(-10, 10)
        imaginary_part = random.uniform(-10, 10)
        complex_number = complex(real_part, imaginary_part)
        numbers.append(complex_number)
    conjugates = [number.conjugate() for number in numbers[::-1]]
    numbers.extend(conjugates)
    return numbers


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
    k = 4
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

    # for i in range(n//2):
    #     print('5^',i,': ', (5**i) % (2*n))

   ########

    gen_encoding_matrix()
    rot_size = five_powers_order // k
    # x_power_rot_size_poly = [0] * rot_size
    # x_power_rot_size_poly.extend([1])
    # x_power_rot_size = Polynomial(x_power_rot_size_poly)
    # print(x_power_rot_size)
    # ctx_compose = [ctx[0].compose(x_power_k), ctx[1].compose(x_power_k)]



    # print("--test fft times inv_fft equals the identity matrix---")
    # mat = np.array(fft_mat)
    # inv_mat = np.linalg.inv(mat)
    # print(np.dot(mat,inv_mat))

    # slots = gen_rnd_slots()
    # slots = [100,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,100]
    # assert len(slots)== n
    # print('---slots---', slots)
    #
    # coeffs = encode(slots)
    # print('---*encode*---',coeffs)
    # new_slots = decode(coeffs)
    # print('---*decode*---', new_slots)
    #
    # ring_elm = RingElement(Polynomial(coeffs), 2*n, q)
    # rot_steps = five_powers_order // k
    # print('ring_elm before: ', ring_elm.poly.coef)
    # ring_elm = ring_elm.automorphism(5 ** rot_steps)
    # print('ring_elm after: ', ring_elm.poly.coef)
    # new_slots = decode(ring_elm.poly.coef)
    # print('----new_slots----', new_slots)

    # for x in new_slots:
    #     print(x)
    #
    #
    # print('---inv_fft_slots scaled---')
    # inv_fft_slots = decode(slots)
    # inv_fft_slots = [round(x.real*scale) for x in inv_fft_slots]
    # for x in inv_fft_slots:
    #     print(x)

    # print('--fft_slots round and scale ---')
    # assert len(inv_fft_slots) == n
    # for x in encode([x / scale for x in inv_fft_slots]):
    #     print(x)

    # ring_elm = RingElement(Polynomial(inv_fft_slots), 2*n, q)
    # print('---encode with rescale ring_elm.poly.coef---')
    # for x in encode(ring_elm.poly.coef):
    #     print(x)

    # print('decode before: ',inv_fft_slots)
    # print('ring_elm before: ', ring_elm.poly.coef)
    # rot_steps = five_powers_order // k
    # print('rot_steps: ',rot_steps)
    # rotated_poly = ring_elm.automorphism(5**rot_steps)
    # rotated_poly.poly.coef = np.append(rotated_poly.poly.coef, 0)  # padding with 0
    # print('rotated_poly after: ', rotated_poly.poly.coef)
    # print('rotated_poly resccaled: ', [x/scale for x in rotated_poly.poly.coef])
    # print('--- test new_slots equal rotated slots---')
    # new_slots = encode([x / scale for x in rotated_poly.poly.coef])
    # for x in new_slots:
    #     print(x)
    # print('------')

    # for i in range(n):
    #     print(encode_mat[i], " ")

    # zeta = cmath.exp(2j * pi / (2*n))
    # five_powers_order = n//2  # the powers of 5 parse half of the primitive roots. Z_{2N}=C_2xC_{N}
    # for i in range(0, five_powers_order):
    #     t = zeta ** (5 **i)
        # print('t: ',t)
    # print('----')
    # for j in range(0,five_powers_order//k):
    #     for i in range(0, k):
    #         coset_representative = auto(zeta, j)
    #         zeta_coset = auto(coset_representative, (five_powers_order / k) * i)
            # print(zeta_coset)
            # if i == k-1:
            #     print('*:', auto(zeta_coset,(five_powers_order/k)))


if __name__ == '__main__':
    main()