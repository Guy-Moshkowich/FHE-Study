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

    def test_sub(self):
        qi = [5, 7]
        tmp1 = coef_to_crt([1, 2, 3, 4, 5, 6, 7, 8], qi)
        tmp2 = coef_to_crt([8, 7, 6, 5, 4, 3, 2, 1], qi)
        res = sub(tmp1, tmp2, qi)
        exp = [28, 30, 32, 34, 1, 3, 5, 7]
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

    def test_mod_up(self):
        qi = [97, 193]
        pi = [101,103]
        poly = [201,202,203,204,205,206,207,208]
        poly_crt = coef_to_crt(poly, qi)
        poly_qipi = mod_up(poly_crt,qi,pi)
        exp = [7, 8, 9, 10, 11, 12, 13, 14, 8, 9, 10, 11, 12, 13, 14, 15, 100, 0, 1, 2, 3, 4, 5, 6, 98, 99, 100, 101,
               102, 0, 1, 2]
        self.assertEqual(poly_qipi, exp)

    def test_gen_relin_key(self):
        qi = [97, 193]
        pi = [101,103]
        qipi = [97,193,101,103]
        sk_qi = gen_sk(qi)
        sk_qipi = mod_up(sk_qi,qi,pi)
        relin_key_ax, relin_key_bx  = gen_relin_key(sk_qi, qi, pi)
        sk_sqr_qi = mul(sk_qi,sk_qi,qi)
        sk_sqr_piqi = mod_up(sk_sqr_qi, qi, pi)
        sk_sqr_times_p_piqi = mul_scalar(P, sk_sqr_piqi, qipi)
        ax_times_sk_qipi = mul(relin_key_ax, sk_qipi, qipi)
        res = sub(relin_key_bx, ax_times_sk_qipi , qipi)
        self.assertEqual(sk_sqr_times_p_piqi, res) # relin encrypts P*sk^2


