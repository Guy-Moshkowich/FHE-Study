from numpy.polynomial import Polynomial
import numpy as np


def dot_prod(v1, v2):
    sum = 0
    for vi, vj in zip(v1, v2):
        sum += vi * vj
    return sum


def bit_decomp(poly: Polynomial, size: int):
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
