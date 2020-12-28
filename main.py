import random
import numpy as np
import sympy
import math


def main():
    n = 10
    p = build_p(n)
    epsilon = 1
    m = round((1 + epsilon)*(n + 1)*math.log(p))
    alpha = 1/(math.sqrt(n)*(math.log(n)**2))
    print('alpha: ', alpha)
    print('sample: ', sample_chi(p, alpha))
    s = sample_uniform(n, p)
    pk = build_public_key(m, s, p, n)
    for i in range(1):
        plaintext = random.randint(0,1)
        print("plaintext: ", plaintext)
        ciphertext = encrypt(plaintext, pk, m, n, p, s)
        print("ciphertext: ", ciphertext)
        dec = decrypt(ciphertext, s, p)
        print("decryption: ", dec)
        print()


def sample_chi(p, alpha):
    normal_sample = np.random.normal(0, alpha, 1)
    return(normal_sample)



def build_p(n):
    p = list(sympy.primerange(n**2, 2*(n**2)))[0]
    return p


def decrypt(ciphertext_bit, s, p):
    a = ciphertext_bit[0]
    b = ciphertext_bit[1]
    d = (b - (np.dot(a,s))) % p
    # print("a: ", a)
    # print("b: ", b)
    #
    # print("minus: ", b - (np.dot(a,s)))
    # print("d: ", d)
    # print("s: ", s)
    # print("dot: ", np.dot(a,s))
    if d - 0 < abs(d - np.floor(p/2)):
        return 0
    else:
        return 1


def encrypt(plaintext_bit, pk, m, n, p,s):
    subset_size = random.randint(0, len(pk))
    subset_m = np.random.choice(m, subset_size)
    sum_a_i = [0]*n
    for i in subset_m:
        sum_a_i = np.add(pk[i][0], sum_a_i)
    sum_b_i = 0
    for i in subset_m:
        sum_b_i += pk[i][1]
    if plaintext_bit == 1:
        sum_b_i += np.floor(p/2)

    d = (sum_b_i - (np.dot(sum_a_i, s))) % p
    # print("ebc_d:", d)
    return [sum_a_i, sum_b_i]


def build_public_key(m, s, p, n):
    pk = []
    for i in range(m):
        a_i = sample_uniform(n, p)
        e_i = sample_uniform_modulo(p)
        # print("e_i: ", e_i)
        b_i = np.dot(a_i, s) + e_i
        # print("b_i-a_i: ", b_i-np.dot(a_i, s))
        # print("pk_bi: ", b_i)
        pk.append([a_i, b_i])
    print("pk:", pk)
    return pk


def sample_uniform_modulo(p):
    return random.choice(range(0, p))


def sample_uniform(n, p):
    return [sample_uniform_modulo(p) for i in range(n)]




if __name__ == '__main__':
    main()
