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


def linear_trans_bsgs_transpose(diags, vec):
    d = 4
    d1 = d-1
    cnt = 0
    delta = ceil(sqrt(2*d-1))
    I = ceil((2*d-1)/delta)
    outer_sum = [0]*len(vec)
    for i in range(I):
        inner_sum = [0]*len(vec)
        for j in range(delta):
            cnt += 1
            if cnt > 2*d-1:
                break
            a = rot(diags[(d1*i*delta) + (d1*(j-d+1))], -d1*i*delta)
            b = rot(vec, d1*(j - d + 1))
            c = mul(a, b)
            inner_sum = add(inner_sum, c)
        outer_sum = add(outer_sum, rot(inner_sum, d1*i*delta))
    return outer_sum


def main():

    diags = {-9: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            -6: [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
            -3: [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
            0: [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
            3: [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            6: [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            9: [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    }
    vec = []
    for i in range(16):
        vec.append(i)
    expected = [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15 ]
    assert(linear_trans(diags, vec) == expected)
    assert(linear_trans_bsgs_transpose(diags, vec) == expected)


if __name__ == '__main__':
    main()