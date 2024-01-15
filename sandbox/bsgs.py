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


def linear_trans_bsgs(diags, vec):
    d = 4
    delta = ceil(sqrt(2*d-1))
    I = ceil((2*d-1)/delta)
    print('d=', d)
    print('I=', I)
    print('delta=', delta)
    outer_sum = [0]*len(vec)
    for i in range(I):
        inner_sum = [0]*len(vec)
        for j in range(delta):
            if i*delta+j-d+1 >= d:
                break
            a = rot(diags[i*delta+j-d+1], -i*delta)
            b = rot(vec, j - d + 1)
            c = mul(a, b)
            inner_sum = add(inner_sum, c)
        outer_sum = add(outer_sum, rot(inner_sum, i*delta))
    return outer_sum

def main():
    diags = {-3: [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
             -2: [0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0],
             -1: [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1],
             0: [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
             1: [0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0],
             2: [0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0],
             3: [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
             }
    vec = []
    for i in range(16):
        vec.append(i)
    expected = [0, 1, 2, 3, 5, 6, 7, 4, 10, 11, 8, 9, 15, 12, 13, 14]
    assert(linear_trans(diags, vec) == expected)
    assert(linear_trans_bsgs(diags, vec) == expected)



if __name__ == '__main__':
    main()