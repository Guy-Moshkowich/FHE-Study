import unittest
from bgv import RingElement
from numpy.polynomial import Polynomial


class TestRingElement(unittest.TestCase):

    def setUp(self):
        self.m = 2**2
        self.q = 10

    def test_modulo_number(self):
        r1 = RingElement(Polynomial([11, 12]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2,0]), self.m, self.q)
        self.assertEqual(r1, r2)

    def test_modulo_phi(self):
        poly = Polynomial([1, 2, 3])
        r1 = RingElement(poly, self.m, self.q)
        r2 = RingElement(poly + 5* r1.get_cyclotomic(self.m), self.m, self.q)
        self.assertTrue(r1 == r2)

    def test_add(self):
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2]), self.m, self.q)
        expected = RingElement(Polynomial([3, 5]), self.m, self.q)
        self.assertTrue(expected == r1 + r2)

    def test_multi(self):
        # (2+3x)(1+2x)=2 + 7x + 6x^2
        # (6 + 7x) + (x^2 + 1)6
        r1 = RingElement(Polynomial([2, 3]), self.m, self.q)
        r2 = RingElement(Polynomial([1, 2]), self.m, self.q)
        expected = RingElement(Polynomial([6, 7]), self.m, self.q)
        self.assertTrue(expected == r1 * r2)

