from numpy.polynomial import Polynomial
import random
import numpy as np


class RingElement:

    def __init__(self, poly, m, mod):
        self.mod = mod
        self.m = m
        self.phi_m = self.get_cyclotomic()
        self.poly = self.modulo(poly, mod)

    def get_cyclotomic(self):
        phi_m = [0]*(self.m//2 + 1)
        phi_m[0] = 1
        phi_m[-1]=1
        return Polynomial(phi_m)

    def __add__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        result = self.modulo(self.poly + other.poly, self.mod)
        return RingElement(result, self.m, self.mod)

    def __mul__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        result = self.modulo(self.poly * other.poly, self.mod)
        return RingElement(result, self.m, self.mod)

    def random(self):
        poly =  Polynomial([random.randrange(0, self.mod) for i in range(0, self.m // 2)])
        return RingElement(poly, self.m, self.mod)

    def change_modulo(self, other, new_modulo):
        return

    def modulo(self, poly, mod):
        poly_modulo_phim = np.flip(np.polydiv(np.flip(poly.coef), np.flip(self.phi_m.coef))[1])
        poly_modulo_phim_q =  Polynomial([x % mod for x in poly_modulo_phim])
        return poly_modulo_phim_q

    def __str__(self):
        return str(self.poly)

    def __eq__(self, other):
        return all(v == 0 for v in (self.poly-other.poly).coef)
