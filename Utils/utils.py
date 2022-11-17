from numpy.polynomial import Polynomial
import numpy as np
from cmath import exp, pi
import math
import random


def dot_prod(v1, v2):
    sum = 0
    for vi, vj in zip(v1, v2):
        sum += vi * vj
    return sum


def bit_decomp(poly_list, size: int):
    bit_decomp_poly_list = []
    for poly in poly_list:
        bit_decomp_poly_list.append(bit_decomp_poly(poly, size))
    list_of_poly_list = np.array(bit_decomp_poly_list).transpose()
    flat_list = [item for sublist in list_of_poly_list for item in sublist]
    return flat_list


def bit_decomp_poly(poly: Polynomial, size: int) -> list[Polynomial]:
    result = []
    for c in poly.coef:
        result.append(bit_decomp_int(c, size))
    result = np.array(result).transpose()
    return [Polynomial(x) for x in result]


def bit_decomp_int(z: int, size: int):
    result = []
    while size > 0:
        size = size - 1
        lsb = z % 2
        result.append(lsb)
        z = (z - lsb)//2
    return result


def powers_of_2(poly_list, size: int) -> list[Polynomial]:
    result = []
    for i in range(size):
        result.extend([2** i * x for x in poly_list])
    return result


def get_nth_primitive_roots_of_unity(n):
    c = 2j * pi / n
    return [exp(k * c) for k in range(n) if math.gcd(k, n) == 1]


def get_canonical_error(ctx, sk, plaintext_expected):
    plaintext_actual = ctx.decrypt(sk)
    diff = plaintext_actual - plaintext_expected
    return diff.canonical_norm()


def ceil(poly: Polynomial) -> Polynomial:
    return Polynomial([math.ceil(c) for c in poly.coef])


def recenter(element: int, mod: int) -> int:
    if element <= mod // 2:
        return element
    else:
        return -abs(mod - element)


def recenter_polynomial(poly : Polynomial, mod: int) -> Polynomial:
    new_coeffs = []
    for c in poly.coef:
        new_coeffs.append(recenter(c, mod))
    return Polynomial(new_coeffs)


def canonical_norm(p: Polynomial, cyclotomic_index: int):
    eval_abs = []
    for val in get_nth_primitive_roots_of_unity(cyclotomic_index):
        eval_result = np.polynomial.polynomial.polyval(val, p.coef)
        eval_abs.append(abs(eval_result))
    return max(eval_abs)


# computes f(x) modulo g(x)
def modulo_polynomial(f: Polynomial, g: Polynomial) -> Polynomial:
    return Polynomial(np.flip(np.polydiv(np.flip(f.coef), np.flip(g.coef))[1]))


# returns X^{m/2} + 1
def build_cyclotomic_poly(cyclotomic_index: int) -> Polynomial:
    assert math.log2(cyclotomic_index) == math.floor(math.log2(cyclotomic_index))
    phi_m = [0]*(cyclotomic_index//2 + 1)
    phi_m[0] = 1
    phi_m[-1] = 1
    return Polynomial(phi_m)


def generate_ternary_polynomial(degree: int, hamming_weight: int) -> Polynomial:
    arr = [0] * degree
    for i in range(hamming_weight):
        random_index = random.randrange(degree)
        if random.randrange(2) == 1:
            arr[random_index] = 1
        else:
            arr[random_index] = -1
    return Polynomial(arr)


def modulo_int(poly : Polynomial, mod: int):
    mod_coeffs = []
    for i in poly.coef:
        mod_coeffs.append(i % mod)
    return recenter_polynomial(Polynomial(mod_coeffs), mod)
