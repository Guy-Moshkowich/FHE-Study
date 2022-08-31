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

    def __sub__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        result = self.modulo(self.poly - other.poly, self.mod)
        return RingElement(result, self.m, self.mod)

    def __mul__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        result = self.modulo(self.poly * other.poly, self.mod)
        return RingElement(result, self.m, self.mod)

    def __str__(self):
        return '[' + str(self.poly) + ', ' + str(self.m) + ', ' + str(self.mod) + ']'

    def __eq__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        return all(v == 0 for v in (self.poly-other.poly).coef)

    @classmethod
    def random(cls, m, mod, max_range=0):
        if max_range == 0:
            max_range = m//2
        poly = Polynomial([random.randint(0, max_range) for i in range(0, m)])
        return cls(poly, m, mod)

    def change_modulo(self, new_modulo):
        return RingElement(self.poly, self.m, new_modulo)

    def modulo(self, poly, mod):
        poly_modulo_phim = np.flip(np.polydiv(np.flip(poly.coef), np.flip(self.phi_m.coef))[1])
        poly_modulo_phim_q =  Polynomial([x % mod for x in poly_modulo_phim])
        return poly_modulo_phim_q