from math import *


def rot(vec, k):
    n = len(vec)
    k = -k % n
    return vec[-k:] + vec[:-k]

def mul(vec1, vec2):
    return [a * b for a, b in zip(vec1, vec2)]


def add(vec1,vec2):
    return [a + b for a, b in zip(vec1, vec2)]


def linear_trans(diags, vec):
    sum = [0]*len(vec)
    for step, diag in diags.items():
        rot_vec = rot(vec, step)
        prod = mul(diag, rot_vec)
        sum = add(sum, prod)
    return sum


def linear_trans_bsgs_tau(diags, vec):
    d = 4
    cnt = 0
    delta = ceil(sqrt(d-1))
    outer_sum = [0]*len(vec)
    for i in range(delta):
        inner_sum = [0]*len(vec)
        for j in range(delta):
            cnt += 1
            if cnt > d:
                break
            a = rot(diags[d*i*delta + j*d], -d*i*delta)
            b = rot(vec, j*d)
            c = mul(a, b)
            inner_sum = add(inner_sum, c)
        outer_sum = add(outer_sum, rot(inner_sum, d*i*delta))
    return outer_sum


def main():

    diags = {0: [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
            4: [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
            8: [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
            12: [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]
    }
    vec = []
    for i in range(16):
        vec.append(i)
    expected = [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11]
    assert(expected == linear_trans(diags, vec))
    assert(expected == linear_trans_bsgs_tau(diags, vec))


if __name__ == '__main__':
    main()