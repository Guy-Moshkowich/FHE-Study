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
        result = bit_decomp_poly(Polynomial([6]), 3)
        expected = [Polynomial([0]), Polynomial([1]), Polynomial([1])]
        self.assertEqual(expected, result)

    def test_bit_decomp_poly(self):
        orig_poly = Polynomial([3,4,5])
        result = bit_decomp_poly(orig_poly, 3) #3+4x+5x^2
        # print('result: ', result)
        reconstruct = Polynomial([0])
        for i in range(len(result)):
            for a in result[i].coef:
                self.assertTrue(a == 1 or a == 0)
            reconstruct += 2**i * result[i]
        self.assertEqual(reconstruct, orig_poly)


    def test_bit_decomp_poly_large_size(self):
        orig_poly = Polynomial([3,4,5])
        result = bit_decomp_poly(orig_poly, 10) #3+4x+5x^2
        reconstruct = Polynomial([0])
        for i in range(len(result)):
            reconstruct += 2**i * result[i]
        self.assertEqual(reconstruct, orig_poly)

    def test_bit_decomp(self):
        size = 8
        orig_poly_list = [Polynomial([3,4,5]), Polynomial([1,2,3])]
        result = bit_decomp(orig_poly_list, size)
        for j in range(len(orig_poly_list)):
            reconstruct = Polynomial([0])
            for l in range(size):
                reconstruct += 2**l * result[j + len(orig_poly_list)*l]
            self.assertEqual(reconstruct, orig_poly_list[j])

    def test_powers_of_2(self):
        size = 2
        orig_poly_list = [Polynomial([1, 2, 3]), Polynomial([4, 5, 6])]
        result = powers_of_2(orig_poly_list, size)
        expected = [Polynomial([1, 2, 3]), Polynomial([4, 5, 6]), Polynomial([2, 4, 6]), Polynomial([8, 10, 12])]
        self.assertEqual(expected, result)