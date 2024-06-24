import unittest
from crt import *


class TestCrt(unittest.TestCase):

    def test_crt_to_coef(self):
        poly = [0, 96, 0, 1, 1, 0, 1, 96]
        poly_crt = coef_to_crt(poly, qi)
        exp = [0, 96, 0, 1, 1, 0, 1, 96, 0, 96, 0, 1, 1, 0, 1, 96]
        self.assertEqual(poly_crt, exp)

    def test_coef_to_crt(self):
        qi = [97, 193]
        poly_crt = [0, 96, 0, 1, 1, 0, 1, 96, 0, 96, 0, 1, 1, 0, 1, 96]
        poly_coef = crt_to_coef(poly_crt, qi)
        exp = [0, 96, 0, 1, 1, 0, 1, 96]
        self.assertEqual(poly_coef, exp)

    def test_add(self):
        qi = [5,7]
        tmp1 = coef_to_crt([1,2,3,4,5,6,7,8], qi)
        tmp2 = coef_to_crt([8,7,6,5,4,3,2,1], qi)
        res = add(tmp1, tmp2, qi)
        exp = [9,9,9,9,9,9,9,9]
        self.assertEqual(crt_to_coef(res, qi), exp)


    def test_mul(self):
        qi = [97, 193]
        tmp1 = coef_to_crt([2, 0, 0, 0, 0, 0, 0, 0], qi)
        tmp2 = coef_to_crt([3, 0, 0, 0, 0, 0, 0, 0], qi)
        res = mul(tmp1, tmp2, qi)
        exp = [6,0,0,0,0,0,0,0]
        self.assertEqual(crt_to_coef(res, qi), exp)

    def test_mul_scalar(self):
        qi = [97, 193]
        poly = coef_to_crt([2, 0, 0, 0, 0, 0, 0, 0], qi)
        scalar = 2
        res = mul_scalar(scalar, poly, qi)
        exp = [4, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(crt_to_coef(res, qi), exp)

