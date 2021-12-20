import random
from numpy.polynomial import Polynomial
import numpy as np
random.seed(0)


class BGV:
    minor_range = 5

    def __init__(self, m_power, q, p, N):
        self.m = 2**m_power # Phi_m(X)=x^m+1
        self.q = q
        self.p = p
        self.N = N
        self.phi_m = self.get_cyclotomic()
        t = self.get_random_poly(self.q)
        self.secret_key = (1, t)
        self.pk_A = self.generate_public_key()

    def generate_public_key(self):
        A = []
        for i in range(self.N):
            b = self.get_random_poly(self.q)
            e = self.get_random_poly(self.minor_range)
            A.append([self.modulo(b*self.secret_key[1]+2*e, self.q), -b])
        return A

    def get_cyclotomic(self):
        phi_m = [0]*(self.m//2 + 1)
        phi_m[0] = 1
        phi_m[-1]=1
        return Polynomial(phi_m)

    def encrypt(self, plaintext):
        c1 = self.get_random_poly(self.q)
        major_noise = self.get_major_noise(c1)
        minor_noise = Polynomial([random.randrange(0, self.minor_range) for i in range(0, self.m // 2)])
        c0 = self.modulo(plaintext - major_noise + self.p*minor_noise, self.q)
        return Ciphertext(c0, c1)

    def get_random_poly(self, coeff_range):
        return Polynomial([random.randrange(0, coeff_range) for i in range(0, self.m // 2)])

    def dot_prod(self, v1, v2):
        sum = 0
        for vi,vj in zip(v1, v2):
            sum += vi*vj
        return sum

    def get_major_noise(self, val):
        return self.modulo(self.secret_key[1] * val, self.q)

    def decrypt(self, ciphertext):
        noise = self.get_major_noise(ciphertext.c1)
        return self.modulo(self.modulo(ciphertext.c0 + noise, self.q), self.p)

    def modulo(self, poly, mod):
        poly_modulo_phim = np.flip(np.polydiv(np.flip(poly.coef), np.flip(self.phi_m.coef))[1])
        return Polynomial([x % mod for x in poly_modulo_phim])

    def add(self, ctx1, ctx2):
        ctx3_0 = self.modulo(ctx1.c0 + ctx2.c0, self.q)
        ctx3_1 = self.modulo(ctx1.c1 + ctx2.c1, self.q)
        return Ciphertext(ctx3_0, ctx3_1)


class Ciphertext:

    def __init__(self, c0, c1):
        self.c0 = c0
        self.c1 = c1

    def __add__(self, other):
        return Ciphertext(self.c0 + other.c0, self.c1 + other.c1)

    def __str__(self):
        return 'c0: ' + str(self.c0) + ', c1: ' + str(self.c1)