import cmath
import numpy as np
from typing import List


def decode(U: List[List[complex]],
           U_conj: List[List[complex]],
           plaintext: List[int]):
    result = []
    for row in U:
        result.append(np.inner(row, plaintext))
    for row in U_conj:
        result.append(np.inner(row, plaintext))
    return result


def encode(U: List[List[complex]],
           U_conj: List[List[complex]],
           dim: int,
           slots: List[complex]):
    tmp1 = []
    tmp2 = []
    for row in np.array(U_conj).transpose():
        tmp1.append((1/(dim//2))*np.inner(row, slots))
    for row in np.array(U).transpose():
        tmp2.append((1/(dim//2))*np.inner(row, [np.conj(slot) for slot in slots]))
    result = [tmp1[i]+tmp2[i] for i in range(len(tmp1))]
    return result


def generate_canonical_power_of_five_cosets(dim: int, k: int):
    U = [[0] * (dim // 2) for _ in range(dim // 4)]
    U_conj = [[0] * (dim // 2) for _ in range(dim // 4)]
    zeta = cmath.exp((2 * cmath.pi * 1j) / dim)
    five_powers_order = dim // 4  # size of max rotation

    # generate roots with special powers
    roots = []
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