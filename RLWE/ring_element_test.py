import unittest
import sys
import os
sys.path.append(os.path.abspath("../RLWE/"))
sys.path.append(os.path.abspath("../Utils/"))
from ring_element import RingElement
from numpy.polynomial import Polynomial
from Utils import utils
import numpy as np

class TestRingElement(unittest.TestCase):

    def setUp(self):
        self.m = 2**2
        self.q = 10

    def test_modulo_number(self):
        r1 = RingElement(Polynomial([11, 12]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2,0]), self.m, self.q)
        self.assertEqual(r1, r2)

    def test_modulo_round_around_zero(self):
        r1 = RingElement(Polynomial([9]), self.m, self.q)
        r2 = RingElement(Polynomial([-1]), self.m, self.q)
        self.assertEqual(r1, r2)

    def test_modulo_phi(self):
        poly = Polynomial([1, 2, 3])
        r1 = RingElement(poly, self.m, self.q)
        r2 = RingElement(poly + 5 * utils.build_cyclotomic_poly(self.m), self.m, self.q)
        self.assertTrue(r1 == r2)


    def test_add(self):
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2]), self.m, self.q)
        expected = RingElement(Polynomial([3, 5]), self.m, self.q)
        self.assertTrue(expected == r1 + r2)

    def test_add_with_rounding(self):
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 12]), self.m, self.q)
        expected = RingElement(Polynomial([3, 5]), self.m, self.q)
        self.assertTrue(expected == r1 + r2)

    def test_sub(self):
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2]), self.m, self.q)
        expected = RingElement(Polynomial([1, 1]), self.m, self.q)
        self.assertTrue(expected == r1 - r2)

    def test_multi(self):
        # (2+3x)(1+2x)=2 + 7x + 6x^2
        # (6 + 7x) + (x^2 + 1)6
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2]), self.m, self.q)
        expected = RingElement(Polynomial([6, 7]), self.m, self.q)
        self.assertTrue(expected == r1 * r2)

    def test_change_modulo(self):
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        result = r1.change_modulo(2)
        expected = RingElement(Polynomial([0, 1]), self.m, 2)
        self.assertEqual(expected, result)

    def test_compose(self):
        r = RingElement(Polynomial([15, 123, 7, 99]), 8, 256)
        x_power_3 = Polynomial([0, 0, 0, 1])
        result = r.compose(x_power_3)
        expected = RingElement(Polynomial([15, 99, 249, 123]), 8, 256)
        self.assertEqual(expected, result)

    def test_norm_canonical(self):
        m = 2**10
        mod = 10000
        r = RingElement.small_gauss(m, mod)
        self.assertTrue(r.canonical_norm() <= 20, "error is to big: error=" + str(r.canonical_norm()))

    def test_automorphism(self):
        dim = 16
        h = 4
        q = 1000
        t = utils.generate_ternary_polynomial(dim // 2, h)
        sk = RingElement(t, dim, q)
        U, U_conj = utils.generate_canonical(dim)
        tmp1 = list(sk.poly.coef)
        tmp1.extend([0 for i in range(len(tmp1), dim//2)])
        slots = utils.decode(U, U_conj, tmp1)
        sk_new = sk.automorphism(5)
        tmp2 = list(sk_new.poly.coef)
        tmp2.extend([0 for i in range(len(tmp2), dim//2)])
        slots_new = utils.decode(U, U_conj, tmp2)
        for slot in slots:
            found = False
            for slot_new in slots_new:
                if np.abs(slot - slot_new) > 0.0001:
                    found = True
                    break
            self.assertTrue(found)






    # def test_ae_size_for_binary_a(self):
    #     n = 1024
    #     q = 1000
    #     a = RingElement.random_binary(n, q)
    #     a_div_2 = Polynomial([c/4 for c in a.poly.coef])
    #     e = RingElement.small_gauss(n, q)
    #     a_div_2_times_e = RingElement(utils.round(a_div_2 * e.poly), n, q)
    #     self.assertTrue((a_div_2_times_e).canonical_norm() < 50)

    # def test_ae_size_for_non_binary_a(self):
    #     n = 1024
    #     q = 100
    #     a = RingElement.random(n, q)
    #     e = RingElement.small_gauss(n, q)
    #     a_bit_decomp_list = utils.bit_decomp_poly(a.poly, q)
    #     ae = RingElement.const(0, n, q)
    #     for a_binary in a_bit_decomp_list:
    #         print('a_binary: ', a_binary)
    #         a_binary_div_2 = Polynomial([c/2 for c in a_binary.coef])
    #         a_binary_div_2_times_e = RingElement(utils.round(a_binary_div_2 * e.poly), n, q)
    #         print('a_binary_div_2_times_e_norm():', a_binary_div_2_times_e.canonical_norm())
    #
    #         ae += a_binary_div_2_times_e
    #     print('(ae).canonical_norm():', (ae).canonical_norm())
    #     self.assertTrue((ae).canonical_norm() < 50)


    # def test_norm_canonical_binary_a(self):
    #     a = RingElement.random_binary(m=2**10, mod=100)
    #     e = RingElement.small_gauss(m=2**10, mod=100)
    #     print("e.canonical_norm(): ", e.canonical_norm());
    #     print("a.canonical_norm(): ", a.canonical_norm());
    #     print("ae.canonical_norm(): ", (a*e).canonical_norm());
