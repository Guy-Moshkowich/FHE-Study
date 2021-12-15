import random
from numpy.polynomial import Polynomial
import numpy as np


p = 2
q = 1000
m = 16 # Phi_16(X)=x^8+1
secret_key = Polynomial([random.randrange(0,q) for i in range(0, m//2)])

def get_cyclotomic(m):
    phi_m = [0]*(m//2 + 1)
    phi_m[0] = 1
    phi_m[-1]=1
    return Polynomial(phi_m)

phi_m = get_cyclotomic(m)


class Ciphertext:
    plaintext = []
    c0 = []
    c1 = []

    def __init__(self, plaintext):
        self.plaintext = plaintext
        self.c1 = Polynomial([random.randrange(0, q) for i in range(0, m//2)])
        major_noise = self.get_noise(self.c1)
        minor_noise = Polynomial([random.randrange(0, 5) for i in range(0, m//2)])
        self.c0 = self.mod_q(plaintext - major_noise + minor_noise)

    def get_noise(self, val):
        multi_modulo_phim = np.polydiv((secret_key * val).coef, phi_m.coef)[1]
        multi_modulo_phim_modulo_q = self.mod_q(Polynomial(multi_modulo_phim))
        return multi_modulo_phim_modulo_q

    def mod_q(self, poly):
        return Polynomial([x % q for x in poly.coef])

    def noise_level(self):
        diff = [c_i - p_i for p_i, c_i in zip(self.plaintext, self.ciphertext)]
        return max(diff)

    def decrypt(self):
        noise = self.get_noise(self.c1)
        return self.c0 + noise

def encrypt(plaintext):
    return Ciphertext(plaintext)


def decrypt(ciphertext):
    return ciphertext.decrypt()


