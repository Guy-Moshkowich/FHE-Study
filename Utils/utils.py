from numpy.polynomial import Polynomial
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


def get_canonical_error(ctx, context, plaintext_expected):
    plaintext_actual = context.decrypt(ctx)
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


def generate_canonical(dim: int):
    U = [[0] * (dim // 2) for i in range(dim // 4)]
    U_conj = [[0] * (dim // 2) for i in range(dim // 4)]
    for k in range(dim // 4):
        root_k = cmath.exp((2 * cmath.pi * 1j * (2 * k + 1)) / dim)
        for t in range(dim // 2):
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

def generate_canonical_power_of_five(dim: int):
    U = [[0] * (dim // 2) for _ in range(dim // 4)]
    U_conj = [[0] * (dim // 2) for _ in range(dim // 4)]
    zeta = cmath.exp((2 * cmath.pi * 1j) / dim)
    for k in range(dim // 4):
        root_k = zeta**((5**k))
        for t in range(dim // 2):
            U[k][t] = root_k ** t
            U_conj[k][t] = np.conj(U[k][t])
    return U, U_conj


def generate_canonical_power_of_five_cosets(dim: int, k: int):
    U = [[0] * (dim // 2) for _ in range(dim // 4)]
    U_conj = [[0] * (dim // 2) for _ in range(dim // 4)]
    zeta = cmath.exp((2 * cmath.pi * 1j) / dim)
    five_powers_order = dim // 4  # size of max rotation
    # k = 2 # size of sub rotation
    roots = []
    # generate roots with special powers
    for j in range(0, five_powers_order // k):
        for i in range(0, k):
            coset_representative = zeta ** (5 ** j)
            zeta_coset = coset_representative**(5**((five_powers_order / k) * i))
            roots.append(zeta_coset)
    assert len(roots) == dim // 4
    for l in range(dim // 4):
        for t in range(dim // 2):
            U[l][t] = roots[l] ** t
            U_conj[l][t] = np.conj(U[l][t])
    return U, U_conj




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


# def get_primitive_complex_root(n):
#     g = units_cyc_gen(n)
#     return = exp(2j * pi / n)

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


def generate_primitive_root_of_unity(n: int, p: int)-> int:
    assert n < p
    assert p % n == 1 # n|p-1 <=> p-1=0 mod n <=> p=1 mod n
    assert is_power_2(n)
    for x in range (2, p-1):
        i = 2
        while (x ** i) % p != 1:
            i = i + 1
        if i == n:
            return x
    assert 1 == 0 # we should not get here.

def is_power_2(n):
    return (n != 0) and ((n & (n - 1)) == 0)

#
# def special_vandermonde(n:int, roots):
#     mat = []
#     for i in range(0, n):
#         row = []
#         for k in range(0, n):
#             row.append((x**i) ** k)
#         mat.append(row)
#     return mat


def special_fft_matrix(n: int):
    x = cmath.exp((2j * pi) / (2*n))
    roots = [x**i for i in get_units(2*n)]
    mat = []
    for root in roots:
        row = [(root ** k) for k in range(0, n)]
        mat.append(row)
    return mat


def special_fft(n: int, p: Polynomial) -> typing.List[complex]:
    mat = np.array(special_fft_matrix(n))
    coef = list(p.coef)
    coef.extend([0 for i in range(n-len(coef))]) # padding with 0
    return list(mat.dot(np.array(coef)))


def inv_special_fft_matrix(n: int):
    mat = np.array(special_fft_matrix(n))
    return np.linalg.inv(mat)

def inv_special_fft(n: int, a: typing.List[complex]) -> Polynomial:
    mat = np.array(special_fft_matrix(n))
    inv_mat = np.linalg.inv(mat)
    coef = inv_mat.dot(np.array(a))
    return Polynomial(coef)


def ntt_matrix(n: int, p: int):
    w = generate_primitive_root_of_unity(2*n, p)
    roots = [w**i for i in get_units(2*n)]
    mat = []
    for root in roots:
        row = [(root ** k) % p for k in range(0, n)]
        mat.append(row)
    return mat

def ntt_mat(p:int, a: typing.List[int], ntt_mat)->typing.List[int]:
    vec = np.array(a)
    mat_mul = np.array(ntt_mat).dot(vec)
    return [(mat_mul[i] % p) for i in range(len(mat_mul))]


def ntt(n: int, p:int, a: typing.List[int])->typing.List[int]:
    mat = np.array(ntt_matrix(n, p))
    return ntt_mat(p, a, mat)


def inv_ntt_mat(p:int, a: typing.List[int], inv_ntt_mat)->typing.List[int]:
    vec = np.array(a)
    vec_res = np.array(inv_ntt_mat).dot(vec)
    res = [0] * len(vec_res)
    for i in range(len(vec_res)):
        res[i] = (vec_res[i] % p)
        if res[i] > p // 2:
            res[i] = res[i] - p
    return res


def inv_ntt(n: int, p:int, a: typing.List[int])->typing.List[int]:
    mat = np.array(inv_ntt_matrix(n, p))
    return inv_ntt_mat(p, a, mat)


def inv(x: int, p: int)->int:
    if x == 1:
        return x
    for g in range(2, p):
        if g * x % p == 1:
            break
    assert g * x % p == 1, 'g='+str(g)+', x='+str(x)
    return g


def inv_ntt_matrix(n: int, p: int):
    two_n_inv = inv(n, p)
    x = generate_primitive_root_of_unity(2 * n, p)
    inv_gen = inv(x, p)
    roots = [inv_gen**i % p for i in get_units(2 * n)]
    mat = []
    for k in range(0,n):
        row = [(two_n_inv * (root ** k)) % p for root in roots]
        mat.append(row)
    return mat


def bit_reverse(num, num_bits: int):
    reversed_num = 0

    for _ in range(num_bits):
        reversed_num <<= 1
        reversed_num |= num & 1
        num >>= 1

    return reversed_num


def read_uint64_from_file(filename):
    data = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Convert each line to an integer and append to the list
                data.append(int(line.strip()))
    except FileNotFoundError:
        print(f"Error: The file {filename} does not exist.")
    except ValueError:
        print("Error: All lines in the file must contain only numbers.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return data



# def inplaceNegacyclicNTT( N:int, vals: typing.List[complex]):
#     t = N #N
#     n = int(math.log(N))
#     logt = n+1
#     m=1
#     root2N = 2j * pi / 2*N
#     while m < N:
#         t = int(t >> 1)
#         logt = logt-1
#         for i in range(0, m):
#             j1 = i * (2**logt)
#             j2 = j1+t-1
#             W = root2N**bit_reverse(m+i, n)
#             for j in range(j1, j2+1):
#                 tmp = W*vals[j+t]
#                 vals[j+t] = vals[j]-tmp
#                 vals[j]=vals[j]-tmp
#         m = m << 1
