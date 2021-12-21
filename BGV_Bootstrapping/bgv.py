import random
from numpy.polynomial import Polynomial
random.seed(0)
from ring_element import RingElement


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
        self.pk_A = self.generate_public_key()

    def generate_public_key(self):
        A = []
        for i in range(self.N):
            b = RingElement.random(self.m, self.q)
            e = RingElement.random(self.m, self.q, max_range=self.epsilon)
            A.append([b*self.secret_key[1]+self.int_2*e, self.int_0-b])
        return A

    def encrypt(self, plaintext_2: RingElement):
        plaintext_q = plaintext_2.change_modulo(self.q)
        c1 = RingElement.random(self.m, self.q)
        major_noise = self.get_major_noise(c1)
        minor_noise = RingElement.random(self.m, self.q, self.epsilon)
        c0 = plaintext_q - major_noise + self.int_2 * minor_noise
        return Ciphertext(c0, c1)

    def get_major_noise(self, val):
        return self.secret_key[1] * val

    def decrypt(self, ciphertext):
        noise = self.get_major_noise(ciphertext.c1)
        return (ciphertext.c0 + noise).change_modulo(2)

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


