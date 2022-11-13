import numpy.polynomial.polynomial
from numpy.polynomial import Polynomial
import random
import numpy as np
from Utils import utils


class RingElement:
    primitive_roots = []
    poly: Polynomial

    def __init__(self, poly: Polynomial, m: int, mod: int):
        self.mod = mod
        self.m = m
        self.phi_m = self.get_cyclotomic()
        self.poly = self.modulo(poly, self.mod, self.phi_m)
        if len(RingElement.primitive_roots) == 0:
            RingElement.primitive_roots = utils.get_nth_primitive_roots_of_unity(self.m)

    def get_cyclotomic(self) -> Polynomial: #X^{m/2}+1
        phi_m = [0]*(self.m//2 + 1)
        phi_m[0] = 1
        phi_m[-1] = 1
        return Polynomial(phi_m)

    def __add__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        return RingElement(self.poly + other.poly, self.m, self.mod)

    def __sub__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        return RingElement(self.poly - other.poly, self.m, self.mod)

    def __mul__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        return RingElement(self.poly * other.poly, self.m, self.mod)


    def __eq__(self, other):
        assert self.mod == other.mod
        assert self.m == other.m
        return all(v == 0 for v in (self.poly-other.poly).coef)

    @classmethod
    def random(cls, m, mod, max_range=0):
        if max_range == 0:
            max_range = mod//2  # modulo range -mod/2 .. +mod/2
        coeffs = []
        for i in range(m//2):
            random_bit = random.randint(0, max_range)
            coeffs.append(random_bit)
        poly = Polynomial(coeffs)
        return cls(poly, m=m, mod=mod)

    @classmethod
    def random_ternary(cls, m, mod):
        poly = Polynomial([random.randrange(2) - 1 for i in range(0, m//2)])
        return cls(poly, m=m, mod=mod)

    @classmethod
    def random_binary(cls, m, mod):
        return RingElement.random(m, mod, max_range=1)


    @classmethod
    def small_gauss(cls, m, mod):
        small_poly = Polynomial([int(random.gauss(mu=1, sigma=0.5)-1) for i in range(0, m//2)])
        return cls(small_poly, m, mod)

    def change_modulo(self, new_modulo):
        return RingElement(self.poly, self.m, new_modulo)

    def modulo(self, poly: Polynomial, mod: int, phi_m: Polynomial):
        poly_modulo_phi_m = np.flip(np.polydiv(np.flip(poly.coef), np.flip(phi_m.coef))[1])
        poly_modulo_phi_m_q = [x % mod for x in poly_modulo_phi_m]
        poly_modulo_phi_m_q_recentered = [utils.recenter(x, mod) for x in poly_modulo_phi_m_q]
        return Polynomial(poly_modulo_phi_m_q_recentered)

    def compose(self, g: Polynomial):  # compute composition of self.poly and g mod (X^{m//2}+1, mod)
        g_poly1d = np.poly1d(g.coef[::-1])
        compose_poly1d = np.poly1d(self.poly.coef[::-1])(g_poly1d)
        return RingElement(Polynomial(compose_poly1d.coefficients[::-1]), self.m, self.mod)

    def canonical_norm(self):
        evals_abs = []
        for x in RingElement.primitive_roots:
            eval = numpy.polynomial.polynomial.polyval(x, self.poly.coef)
            evals_abs.append(abs(eval))
        return max(evals_abs)


    def __str__(self):
        return 'first 10 coefficients:' + str(self.poly.coef[:10]) + ", degree: " + str(self.m) + ', modulo: ' + str(self.mod)

    @classmethod
    def const(cls, val, n, q):
        return RingElement(Polynomial([val]), n, q)
