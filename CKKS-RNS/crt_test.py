import unittest
from scheme import *

class TestCrt(unittest.TestCase):

    def test_crt_to_coef(self):
        qi = [97,193]
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

    def test_coef_to_crt_elm(self):
        qi = [97, 193]
        elm = 200
        out = coef_to_crt_elm(elm, qi)
        exp = [6, 7]
        self.assertEqual(out, exp)

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
        pi = [101,103]
        qipi = [97,193,101,103]
        P = 101*103
        sk_qipi = gen_sk(qipi)
        relin_key_ax, relin_key_bx  = gen_relin_key(sk_qipi, qipi, pi)
        sk_sqr_qipi = mul(sk_qipi,sk_qipi,qipi)
        sk_sqr_times_p_qipi = mul_scalar(P, sk_sqr_qipi, qipi)
        ax_times_sk_qipi = mul(relin_key_ax, sk_qipi, qipi)
        res = sub(relin_key_bx, ax_times_sk_qipi , qipi)
        self.assertEqual(sk_sqr_times_p_qipi, res) # relin encrypts P*sk^2

    def test_enc_dec(self):
        qi = [97, 193]
        sk_qi = gen_sk(qi)
        pt_qi = coef_to_crt([1, 2, 3, 4, 5, 6, 7, 8], qi)
        ax_qi, bx_qi = encrypt(pt_qi, sk_qi, qi)
        pt_dec_qi = decrypt(ax_qi, bx_qi, sk_qi, qi)
        self.assertEqual(pt_dec_qi, pt_qi)

    def test_he_add(self):
        qi = [97, 193]
        sk_qi = gen_sk(qi)
        pt1_qi = coef_to_crt([1, 2, 3, 4, 5, 6, 7, 8], qi)
        pt2_qi = coef_to_crt([9, 10, 11, 12, 13, 14, 15, 16], qi)
        exp_qi = [10, 12, 14, 16, 18, 20, 22, 24]
        ax1_qi, bx1_qi = encrypt(pt1_qi, sk_qi, qi)
        ax2_qi, bx2_qi = encrypt(pt2_qi, sk_qi, qi)
        res_ax, res_bx = he_add(ax1_qi, bx1_qi,ax2_qi, bx2_qi, qi)
        res_dec_qi = decrypt(res_ax, res_bx, sk_qi, qi)
        self.assertEqual(crt_to_coef(res_dec_qi,qi), exp_qi)

    def test_mod_down_elm(self):
        qi = [97, 193]
        qipi = [97, 193, 101, 103]
        P = 101*103
        elm_coef = 2*P+2
        elm_qipi = coef_to_crt_elm(elm_coef, qipi)
        elm_qi = mod_down_elm(elm_qipi, qipi, qi)
        exp = round(elm_coef/P)
        self.assertEqual(crt_to_coef_elm(elm_qi, qi), exp)

    def test_mod_down(self):
        qi = [97, 193]
        qipi = [97, 193, 101, 103]
        P=101*103
        poly_coef = [i*P+10*i for i in range(n)]
        poly_qipi = coef_to_crt(poly_coef, qipi)
        poly_qi = mod_down(poly_qipi, qipi, qi)
        exp = [0, 1, 2, 3, 4, 5, 6, 7]
        self.assertEqual(exp, crt_to_coef(poly_qi, qi))




    def test_he_mul(self):

        '''
        import numpy as np

        a = [Q*P-1]
        b = [Q*P-1]
        print(a)

        d =np.polymul(a, b)
        print('d: ', d.tolist())
        cyc2 = np.poly1d([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0])
        quotient, remainder = np.polynomial.polynomial.polydiv(d, cyc2)
        print('res: ', quotient, remainder.tolist())

        c = utils.modulo_polynomial(Polynomial(a)*Polynomial(b), cyc)
        print('c: ',c)
        print('c%97: ',c.coef[0] %97)
        print('np.poly1d(a): ',np.poly1d(a))
        cyc2 = np.poly1d([1,0,0,0,0,0,0,0,1])
        print('cyc2: ')
        print(cyc2)
        print('np.poly1d(a)*np.poly1d(b): ')
        tmp = np.poly1d(a)*np.poly1d(b)
        print(np.poly1d(a)*np.poly1d(b))
        quotient, remainder = np.polydiv(d,cyc2)
        print('res: ', quotient, remainder[0]%97)

        for i, coeff in enumerate(tmp):
            if abs(coeff) >= 1e10 or abs(coeff) < 1e-10:
                print(f"Coefficient of x^{remainder.order - i}: {coeff:.22f}")  # Print with 18 decimal places
            else:
                print(f"Coefficient of x^{remainder.order - i}: {coeff}")
        '''

        # qi = [97, 193]
        # scale = 6 #100
        # qipi = [97, 193, 101, 103]
        qi = [5,7]
        pi = [11,13]
        qipi = [5,7,11,13]
        sk_qipi = gen_sk(qipi, debug=True)
        sk_qi = mod_down(sk_qipi, qipi, qi)
        print("sk_qipi: ", sk_qipi)
        print("sk_qi: ", sk_qi)
        # pt1_qi = coef_to_crt([0,scale,0,0,0,0,0,0], qi)
        # pt2_qi = coef_to_crt([0,0,scale,0,0,0,0,0], qi)
        pt1_qi = coef_to_crt([0, 2, 0, 0, 0, 0, 0, 0], qi)
        pt2_qi = coef_to_crt([0, 0, 3, 0, 0, 0, 0, 0], qi)
        print("pt1_qi: ", pt1_qi)
        print("pt2_qi: ", pt2_qi)

        # exp_qi = [0,0,0,scale*scale,0,0,0,0]
        ax1_qi, bx1_qi = encrypt(pt1_qi, sk_qi, qi, debug=True)
        # print(crt_to_coef(decrypt(ax1_qi, bx1_qi, sk_qi, qi),qi))
        ax2_qi, bx2_qi = encrypt(pt2_qi, sk_qi, qi, debug=True)
        # print(crt_to_coef(decrypt(ax2_qi, bx2_qi, sk_qi, qi),qi))

        relin_key_ax_qipi, relin_key_bx_qipi = gen_relin_key(sk_qipi, qi, qipi, debug=True)

        res_ax_qi, res_bx_qi = he_mul(ax1_qi, bx1_qi, ax2_qi, bx2_qi, relin_key_ax_qipi, relin_key_bx_qipi, qi,pi)
        res_dec_qi = decrypt(res_ax_qi, res_bx_qi, sk_qi, qi)
        print('res=', crt_to_coef(res_dec_qi,qi))
        # self.assertEqual(crt_to_coef(res_dec_qi, qi), exp_qi)
