import random
from numpy.polynomial import Polynomial
import numpy as np


class BGV:
    minor_range = 5

    def __init__(self, m, q, p):
        self.m = m
        self.q = q
        self.p = p
        self.phi_m = self.get_cyclotomic()

    def generate_key(self):
        self.secret_key = Polynomial([random.randrange(0, self.q) for i in range(0, self.m//2)])

    def get_cyclotomic(self):
        phi_m = [0]*(self.m//2 + 1)
        phi_m[0] = 1
        phi_m[-1]=1
        return Polynomial(phi_m)

    def encrypt(self, plaintext):
        c1 = Polynomial([random.randrange(0, self.q) for i in range(0, self.m // 2)])
        major_noise = self.get_noise(c1)
        minor_noise = Polynomial([random.randrange(0, self.minor_range) for i in range(0, self.m // 2)])
        c0 = self.modulo(plaintext - major_noise + self.p*minor_noise, self.q)
        return Ciphertext(c0, c1)

    def get_noise(self, val):
        multi_modulo_phim = np.polydiv((self.secret_key * val).coef, self.phi_m.coef)[1]
        multi_modulo_phim_modulo_q = self.modulo(Polynomial(multi_modulo_phim), self.q)
        return multi_modulo_phim_modulo_q

    def decrypt(self, ciphertext):
        noise = self.get_noise(ciphertext.c1)
        return self.modulo(self.modulo(ciphertext.c0 + noise, self.q), self.p)

    def modulo(self, poly, mod):
        return Polynomial([x % mod for x in poly.coef])

    def add(self, ctx1, ctx2):
        return Ciphertext(ctx1.c0 + ctx2.c0, ctx1.c1 + ctx2.c1)


class Ciphertext:

    def __init__(self, c0, c1):
        self.c0 = c0
        self.c1 = c1

    def __add__(self, other):
        return Ciphertext(self.c0 + other.c0, self.c1 + other.c1)


    def __str__(self):
        return 'c0: ' + str(self.c0) + ', c1: ' + str(self.c1)