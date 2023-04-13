from numpy.polynomial import Polynomial
import numpy as np
from cmath import exp, pi
import math
import random
import cmath
import typing
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




def generate_canonical(n: int):
    U = [[0] * (n//2) for i in range(n // 4)]
    U_conj = [[0] * (n//2) for i in range(n // 4)]
    for k in range(n // 4):
        root_k = cmath.exp((2 * cmath.pi * 1j * (2 * k + 1)) / n)
        for t in range(n//2):
            U[k][t] = root_k ** t
            U_conj[k][t] = np.conj(U[k][t])
    return U, U_conj


def decode(U, U_conj, plaintext):
    result = []
    for row in U:
        result.append(np.inner(row, plaintext))
    for row in U_conj:
        result.append(np.inner(row, plaintext))
    return result


def encode(U, U_conj, dim, slots):
    tmp1 = []
    tmp2 = []
    for row in np.array(U_conj).transpose():
        tmp1.append((1/(dim//2))*np.inner(row, slots))
    for row in np.array(U).transpose():
        tmp2.append((1/(dim//2))*np.inner(row, [np.conj(slot) for slot in slots]))
    result = [tmp1[i]+tmp2[i] for i in range(len(tmp1))]
    return result


def get_units(n :int) -> typing.List[int]:
    return [k for k in range(n) if math.gcd(k, n) == 1]


def order(n :int, unit: int) -> int:
    order = 1
    x = int(unit)
    while x != 1:
        order += 1
        x = x**order % n
    return order


def units_cyc_gen(n: int) -> int:
    units = get_units(n)
    for unit in units:
        if order(n, unit) == len(units)//2:
            return unit


def fft_matrix(n) -> typing.List[typing.List]:
    g = units_cyc_gen(n)
    primitive = exp(2j * pi / n)
    cyc_powers = [g ** i % n for i in range(n // 4)]
    cyc_roots = []
    for k in cyc_powers:
        cyc_roots.append(primitive ** k)
    roots = list(cyc_roots)
    for r in cyc_roots:
        roots.append(np.conj(r))
    mat = []
    for root in roots:
        row = [root**k for k in range(0, n//2)]
        mat.append(row)
    return mat


def fft(n: int, p: Polynomial) -> typing.List[complex]:
    mat = np.array(fft_matrix(n))
    coef = list(p.coef)
    coef.extend([0 for i in range(n//2-len(coef))]) # padding with 0
    return list(mat.dot(np.array(coef)))

def inv_fft(n: int, a: typing.List[complex]) -> Polynomial:
    mat = np.array(fft_matrix(n))
    inv_mat = np.linalg.inv(mat)
    coef = inv_mat.dot(np.array(a))
    return Polynomial(coef)

# def NTT(p:Polynomial, q:int)->typing.List[int]:
#
#     return []