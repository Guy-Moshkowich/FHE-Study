from numpy.polynomial import Polynomial
import numpy as np
from cmath import exp, pi
import math

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


def round(poly) -> Polynomial:
    return Polynomial([math.ceil(c) for c in poly.coef])


def recenter(element: int, mod: int) -> int:
    if element <= mod // 2:
        return element
    else:
        return -abs(mod - element)
