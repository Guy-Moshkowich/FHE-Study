from numpy.polynomial import Polynomial
import numpy as np


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


def bit_decomp_poly(poly: Polynomial, size: int):
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


def powers_of_2(poly_list, size: int):
    result = []
    for i in range(size):
        result.extend([2** i * x for x in poly_list])
    return result
