import unittest
from utils import *


class TestUtils(unittest.TestCase):

    def test_dot_prod(self):
        result = dot_prod([1, 2, 3], [4, 5, 6])
        expected = 1 * 4 + 2 * 5 + 3 * 6
        self.assertEqual(expected, result)

    def test_bit_decomp_int(self):
        self.assertEqual([0, 1, 1], bit_decomp_int(z=6, size=3))

    def test_bit_decomp_int_size(self):
        self.assertEqual([0, 1, 1, 0, 0], bit_decomp_int(z=6, size=5))

    def test_bit_decomp_deg_0(self):
        result = bit_decomp(Polynomial([6]), 3)
        expected = [Polynomial([0]), Polynomial([1]), Polynomial([1])]
        self.assertEqual(expected, result)

    def test_bit_decomp(self):
        orig_poly = Polynomial([3,4,5])
        result = bit_decomp(orig_poly, 3) #3+4x+5x^2
        reconstruct = Polynomial([0])
        for i in range(len(result)):
            for a in result[i].coef:
                self.assertTrue(a == 1 or a == 0)
            reconstruct += 2**i * result[i]
        self.assertEqual(reconstruct, orig_poly)
