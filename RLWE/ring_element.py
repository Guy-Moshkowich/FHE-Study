import numpy.polynomial.polynomial
from numpy.polynomial import Polynomial
import random
import numpy as np
import math
import Utils.utils


class RingElement:
    primitive_roots =[]

    def __init__(self, poly: Polynomial, m, mod):
        self.mod = mod
        self.m = m
        self.phi_m = self.get_cyclotomic()
        self.poly = self.modulo(poly, mod)
        if len(RingElement.primitive_roots) == 0:
            RingElement.primitive_roots = Utils.utils.get_nth_primitive_roots_of_unity(self.m)


    def get_cyclotomic(self): #X^{m/2}+1
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
            max_range = m #TODO: replace with m
        poly = Polynomial([random.randint(0, max_range) for i in range(0, m)])
        return cls(poly, m, mod)

    @classmethod
    def random_binary(cls, m, mod):
        return RingElement.random(m, mod, max_range=1)


    @classmethod
    def small_gauss(cls, m, mod, mu=1, sigma=0.5):
        small_poly = Polynomial([int(random.gauss(mu=1, sigma=0.5)-1) for i in range(0, m//2)])
        return cls(small_poly, m, mod)

    def change_modulo(self, new_modulo):
        return RingElement(self.poly, self.m, new_modulo)

    def modulo(self, poly, mod):
        poly_modulo_phim = np.flip(np.polydiv(np.flip(poly.coef), np.flip(self.phi_m.coef))[1])
        poly_modulo_phim_q =  Polynomial([x % mod for x in poly_modulo_phim])
        return poly_modulo_phim_q

    def compose(self, g: Polynomial):  # compute composition of self.poly and g mod (X^{m//2}+1, mod)
        g_poly1d = np.poly1d(g.coef[::-1])
        compose_poly1d = np.poly1d(self.poly.coef[::-1])(g_poly1d)
        return RingElement(Polynomial(compose_poly1d.coefficients[::-1]), self.m, self.mod)

    def canonical_norm(self):
        evals_abs = []
        for x in RingElement.primitive_roots:
            coef_transformed = [self.transform_mod_div_2_center(x) for x in self.poly.coef]
            eval = numpy.polynomial.polynomial.polyval(x, coef_transformed)
            evals_abs.append(abs(eval))
        return max(evals_abs)

    def transform_mod_div_2_center(self, element):
        if element <= self.mod // 2:
            return element
        else:
            return abs(self.mod - element)

    def __str__(self):
        return str(self.poly.coef[:10])


