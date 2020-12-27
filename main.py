import random
import numpy as np


def main():
    n = 10
    p = 17
    m = 15
    s = sample_uniform(n, p)
    pk = build_public_key(m, s, p, n)
    for i in range(10):
        plaintext = random.randint(0,1)
        print("plaintext: ", plaintext)
        ciphertext = encrypt(plaintext, pk, m, n, p)
        print("ciphertext: ", ciphertext)
        dec = decrypt(ciphertext, s, p)
        print("decryption: ", dec)
        print()


def decrypt(ciphertext_bit, s, p):
    a = ciphertext_bit[0]
    b = ciphertext_bit[1]
    d = (b - np.dot(a,s)) % p
    print("d: ", d)
    if d - 0 < abs(d - np.floor(p/2)):
        return 0
    else:
        return 1

def encrypt(plaintext_bit, pk, m, n, p):
    subset_size = random.randint(0, len(pk))
    subset_m = np.random.choice(m, subset_size)
    sum_a_i = [0]*n
    for i in subset_m:
        sum_a_i = np.add(pk[i][0], sum_a_i)
    sum_a_i = np.mod(sum_a_i, [p]*n)
    sum_b_i = 0
    for i in subset_m:
        sum_b_i += pk[i][1]
    sum_b_i = sum_b_i % p
    if plaintext_bit == 1:
        sum_b_i += np.floor(p/2)
    return [sum_a_i, sum_b_i]


def build_public_key(m, s, p, n):
    pk = []
    for i in range(m):
        a_i = sample_uniform(n, p)
        e_i = sample_uniform_modulo(p)
        b_i = (np.dot(a_i, s) + e_i) % p
        pk.append([a_i, b_i])
    print("pk:", pk)
    return pk


def sample_uniform_modulo(p):
    return random.choice(range(0, p))


def sample_uniform(n, p):
    return [sample_uniform_modulo(p) for i in range(n)]




if __name__ == '__main__':
    main()
