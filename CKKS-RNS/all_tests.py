import unittest
from scheme import *

qi = [1009, 1153]
pi = [1201, 1217]
scheme_1000 = Scheme(n, qi, pi)


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


    def test_mul2(self):
        qi = [97, 113]
        tmp1 = coef_to_crt([0, 2, 0, 0, 0, 0, 0, 0], qi)
        tmp2 = coef_to_crt([0, 0, 3, 0, 0, 0, 0, 0], qi)
        scheme = Scheme(n,qi)
        res = scheme.polyEval.mul(tmp1, tmp2, qi)
        exp = [0,0,0,6,0,0,0,0]
        self.assertEqual(crt_to_coef(res, qi), exp)

    def test_mul_minus(self):
        qi = [97, 113]
        scheme = Scheme(n,qi)

        Q = prod(qi)
        sk = [9796, 0, 0, 0, 0, 0, 0, 0]
        A = [7202, 0, 0, 0, 0, 0, 0, 0]
        sk_qi = coef_to_crt(sk, qi)
        print('sk_qi: ', sk)
        A_p=Polynomial(A)
        print('A_p: ',A_p)
        sk_p = Polynomial(sk)
        print('sk_p: ',sk_p)
        res = scheme.modulo(A*sk_p, Q)
        print('res coef: ', res)
        print('res qi: ', coef_to_crt(res,qi))

        print(9796*7202%9797)

        A_qi = coef_to_crt(A, qi)
        print(A_qi)
        eval = Evaluator(n, qi)
        res2 = eval.mul(A_qi, sk_qi, qi)
        print('res2 qi: ', res2)
        print('res2 coef: ',crt_to_coef(res2,qi))


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
        # qi = [1009, 1153]
        # pi = [1201, 1217]
        # scheme = Scheme(n, qi, pi)

        qipi = qi.copy()
        qipi.extend(pi)
        P = prod(pi)
        sk_qipi = scheme_1000.gen_sk(qipi, debug=Debug.POSITIVE_SK)
        relin_key_ax, relin_key_bx  = scheme_1000.gen_relin_key(sk_qipi, qipi, pi)
        sk_sqr_qipi = scheme_1000.polyEval.mul(sk_qipi, sk_qipi, qipi)
        sk_sqr_times_p_qipi = mul_scalar(P, sk_sqr_qipi, qipi)
        # ax_times_sk_qipi = mul2(relin_key_ax, sk_qipi, qipi)
        ax_times_sk_qipi = scheme_1000.polyEval.mul(relin_key_ax, sk_qipi, qipi)
        res = sub(relin_key_bx, ax_times_sk_qipi , qipi)
        self.assertEqual(sk_sqr_times_p_qipi, res) # relin encrypts P*sk^2

    def test_enc_dec(self):
        qi = [97, 193]
        scheme = Scheme(n, qi)
        sk_qi = scheme.gen_sk(qi,Debug.POSITIVE_SK)
        pt_qi = coef_to_crt([1, 2, 3, 4, 5, 6, 7, 8], qi)
        ax_qi, bx_qi = scheme.encrypt(pt_qi, sk_qi, qi, debug=Debug.DISABLE_NOISE)
        pt_dec_qi = scheme.decrypt(ax_qi, bx_qi, sk_qi, qi)
        self.assertEqual(pt_dec_qi, pt_qi)

    def test_he_add(self):
        qi = [97, 193]
        scheme = Scheme(n, qi)
        sk_qi = scheme.gen_sk(qi,Debug.POSITIVE_SK)
        pt1_qi = coef_to_crt([1, 2, 3, 4, 5, 6, 7, 8], qi)
        pt2_qi = coef_to_crt([9, 10, 11, 12, 13, 14, 15, 16], qi)
        exp_qi = [10, 12, 14, 16, 18, 20, 22, 24]
        ax1_qi, bx1_qi = scheme.encrypt(pt1_qi, sk_qi, qi, debug=Debug.DISABLE_NOISE)
        ax2_qi, bx2_qi = scheme.encrypt(pt2_qi, sk_qi, qi, debug=Debug.DISABLE_NOISE)
        res_ax, res_bx = scheme.he_add(ax1_qi, bx1_qi,ax2_qi, bx2_qi, qi)
        res_dec_qi = scheme.decrypt(res_ax, res_bx, sk_qi, qi)
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
        # qi = [999983, 999979]
        # pi = [1000003, 1000033]

        qi = [1009, 1153]
        pi = [1201, 1217]
        # scheme = Scheme(n, qi, pi)

        # qi = [97, 101]
        # pi = [103, 107]

        # primes near 100 s.t. P = 1 mod 16
        # pi = [97, 113]
        # qi = [183, 241]

        qipi = qi.copy()
        qipi.extend(pi)
        P = prod(pi)
        Q = prod(qi)

        print('qipi: ', qipi)
        print('Q=', Q)
        print('P=', P)
        print('QP=', Q*P)

        sk_qi = scheme_1000.gen_sk(qi, Debug.POSITIVE_SK)
        sk_qipi = mod_up(sk_qi, qi, pi)
        print("sk_qipi: ", sk_qipi)
        print("sk_qi: ", sk_qi)
        print("sk_coef: ", crt_to_coef(sk_qi, qi))

        # pt1_coef = [0, 200, 0, 0, 0, 0, 0, 0]
        # pt2_coef = [0, 0, 30, 0, 0, 0, 0, 0]
        # pt1_qi = coef_to_crt(pt1_coef, qi)
        # pt2_qi = coef_to_crt(pt2_coef, qi)

        pt1_qi = gen_rand_poly_crt(qi)
        pt2_qi = gen_rand_poly_crt(qi)
        pt1_coef = crt_to_coef(pt1_qi,qi)
        pt2_coef = crt_to_coef(pt2_qi,qi)

        print("pt1_qi: ", pt1_qi)
        print("pt2_qi: ", pt2_qi)
        print("pt1_coef: ", crt_to_coef(pt1_qi,qi))
        print("pt2_coef: ", crt_to_coef(pt2_qi,qi))
        ax1_qi, bx1_qi = scheme_1000.encrypt(pt1_qi, sk_qi, qi, debug=Debug.DISABLE_NOISE)
        ax2_qi, bx2_qi = scheme_1000.encrypt(pt2_qi, sk_qi, qi, debug=Debug.DISABLE_NOISE)
        print('bx1: ', crt_to_coef(bx1_qi, qi))
        print('bx2: ', crt_to_coef(bx2_qi, qi))
        print(crt_to_coef(scheme_1000.decrypt(ax1_qi, bx1_qi, sk_qi, qi), qi))
        print(crt_to_coef(scheme_1000.decrypt(ax2_qi, bx2_qi, sk_qi, qi), qi))

        relin_key_ax_qipi, relin_key_bx_qipi = scheme_1000.gen_relin_key(sk_qipi, qipi, pi)

        res_ax_qi, res_bx_qi = scheme_1000.he_mul(ax1_qi, bx1_qi, ax2_qi, bx2_qi, relin_key_ax_qipi, relin_key_bx_qipi, qi, pi, sk_qi, pt1_qi, pt2_qi)
        print(res_ax_qi,res_ax_qi)
        res_dec_qi = scheme_1000.decrypt(res_ax_qi, res_bx_qi, sk_qi, qi)
        res_dec_coef = crt_to_coef(res_dec_qi,qi)
        print('res_dec_qi: ',res_dec_qi)
        print('res=', res_dec_coef)
        exp_coef = scheme_1000.modulo(Polynomial(pt1_coef) * Polynomial(pt2_coef), Q)
        exp_coef.extend([0]*(len(res_dec_coef)-len(exp_coef)))
        print('exp_coef=', exp_coef)
        error = 50
        for i in range(len(res_dec_coef)):
            c = res_dec_coef[i]
            d = exp_coef[i]
            if (c>Q/2):
                c=Q-c
            if (d>Q/2):
                d=Q-d
            assert_cond = abs(c- d) < error
            if not assert_cond:
                print(i, c, d)
            self.assertTrue(assert_cond)

    def test_ntt_multiplication(self):
        # precondition: p mod 2n =1, 17 mod 8 =1
        n = 8
        p = 1009
        a = [0,1,0,0,0,0,0,0] # x
        b = [0,0,1,0,0,0,0,0] # x^2
        ntt_a = utils.ntt(n, p, a)
        ntt_b = utils.ntt(n, p, b)
        actual = [ntt_a[i]*ntt_b[i] % p for i in range(len(a))]
        self.assertEqual([0,0,0,1,0,0,0,0], utils.inv_ntt(n,p,actual))