import random
from numpy.polynomial import Polynomial
from ring_element import RingElement
import numpy as np


class BGV:
    epsilon = 5

    def __init__(self, m_power, q, p, N):
        self.m = 2**m_power # Phi_m(X)=x^m+1
        self.q = q
        self.p = p
        self.N = N
        self.int_2 = RingElement(Polynomial(2), self.m, self.q)
        self.int_0 = RingElement(Polynomial(0), self.m, self.q)
        t = RingElement.random(m=self.m, mod=self.q)
        self.secret_key = (1, t)
        self.pk = self.generate_public_key()

    def generate_public_key(self):
        A = []
        for i in range(self.N):
            b = RingElement.random(self.m, self.q)
            e = RingElement.random(self.m, self.q, max_range=self.epsilon)
            A.append(np.array([b*self.secret_key[1]+self.int_2*e, self.int_0 - b]))
        return np.array(A)

    def encrypt(self, plaintext_2: RingElement):
        plaintext_q = plaintext_2.change_modulo(self.q)
        r = np.array([RingElement.random(self.m, self.q) for i in range(self.N)])
        ctx = np.matmul(self.pk.transpose(), r) + [plaintext_q, self.int_0]
        return Ciphertext(ctx[0], ctx[1])

    def decrypt(self, ciphertext):
        noise = self.secret_key[1] * ciphertext.c1
        return (ciphertext.c0 + noise).change_modulo(2)


class Ciphertext:

    def __init__(self, c0, c1):
        self.c0 = c0
        self.c1 = c1

    def __add__(self, other):
        return Ciphertext(self.c0 + other.c0, self.c1 + other.c1)

    def __sub__(self, other):
        return Ciphertext(self.c0 - other.c0, self.c1 - other.c1)

    def __mul__(self, other):
        # (c0+s*c1)*(d0+s*d1) = c0d0+s*(c1*d0+c0d1)+s^2*c1d1 =
        # = (c0d0,c1*d0+c0d1,c1d1)*(1,s,s^2)
        mult = []
        return

    def __str__(self):
        return 'c0: ' + str(self.c0) + ', c1: ' + str(self.c1)


