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

    def test_mul_minus(self):
        qi = [97, 101]
        Q=97*101
        sk = [9796, 0, 0, 0, 0, 0, 0, 0]
        A = [7202, 0, 0, 0, 0, 0, 0, 0]
        sk_qi = coef_to_crt(sk, qi)
        print('sk_qi: ', sk)
        A_p=Polynomial(A)
        print('A_p: ',A_p)
        sk_p = Polynomial(sk)
        print('sk_p: ',sk_p)
        res = modulo(A*sk_p, Q)
        print('res coef: ', res)
        print('res qi: ', coef_to_crt(res,qi))

        print(9796*7202%9797)

        A_qi = coef_to_crt(A, qi)
        print(A_qi)
        res2 = mul(A_qi, sk_qi, qi)
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

    def test_mod_down_396810(self):
        qi = [97, 101]
        scale = 100
        pi = [103, 107]
        qipi = [97, 101, 103, 107]
        poly_coef = [396810]
        poly_qipi = coef_to_crt(poly_coef, qipi)
        print(poly_qipi)
        poly_qi = mod_down(poly_qipi, qipi, qi)
        print(poly_qi)

    # def test_he_mul_sk_positive_debug_true(self):
    #     qi = [5, 7]
    #     pi = [11, 13]
    #     qipi = [5, 7, 11, 13]
    #     sk_qi = gen_sk(qi,debug=True)
    #     sk_qipi = mod_up(sk_qi, qi, pi)
    #     pt1_qi = coef_to_crt([random.randint(0, 10) for _ in range(n)], qi)
    #     pt2_qi = coef_to_crt([random.randint(0, 10) for _ in range(n)], qi)
    #     ax1_qi, bx1_qi = encrypt(pt1_qi, sk_qi, qi)
    #     ax2_qi, bx2_qi = encrypt(pt2_qi, sk_qi, qi)
    #     relin_key_ax_qipi, relin_key_bx_qipi = gen_relin_key(sk_qipi, qipi, pi)
    #     res_ax_qi, res_bx_qi = he_mul(ax1_qi, bx1_qi, ax2_qi, bx2_qi, relin_key_ax_qipi, relin_key_bx_qipi, qi, pi)
    #     res_dec_qi = decrypt(res_ax_qi, res_bx_qi, sk_qi, qi)
    #     exp = mul(pt1_qi, pt2_qi,qi)
    #     self.assertEqual(exp, res_dec_qi)


    # def test_he_mul_debug(self):
    #     # qi = [97, 193]
    #     # scale = 6 #100
    #     # qipi = [97, 193, 101, 103]
    #     qi = [5,7]
    #     pi = [11,13]
    #     qipi = [5,7,11,13]
    #     sk_qi = gen_sk(qi, debug=True)
    #     sk_qipi = mod_up(sk_qi, qi, pi)
    #     print("sk_qipi: ", sk_qipi)
    #     print("sk_qi: ", sk_qi)
    #     print("sk_coef: ", crt_to_coef(sk_qi,qi))
    #
    #     pt1_qi = coef_to_crt([10, 0, 0, 0, 0, 0, 0, 0], qi)
    #     pt2_qi = coef_to_crt([12, 0, 0, 0, 0, 0, 0, 0], qi)
    #     print("pt1_qi: ", pt1_qi)
    #     print("pt2_qi: ", pt2_qi)
    #
    #     # exp_qi = [0,0,0,scale*scale,0,0,0,0]
    #     ax1_qi, bx1_qi = encrypt(pt1_qi, sk_qi, qi)
    #     print(crt_to_coef(decrypt(ax1_qi, bx1_qi, sk_qi, qi),qi))
    #     ax2_qi, bx2_qi = encrypt(pt2_qi, sk_qi, qi)
    #     print(crt_to_coef(decrypt(ax2_qi, bx2_qi, sk_qi, qi),qi))
    #
    #     relin_key_ax_qipi, relin_key_bx_qipi = gen_relin_key(sk_qipi, qipi, pi)
    #
    #     res_ax_qi, res_bx_qi = he_mul(ax1_qi, bx1_qi, ax2_qi, bx2_qi, relin_key_ax_qipi, relin_key_bx_qipi, qi,pi)
    #     res_dec_qi = decrypt(res_ax_qi, res_bx_qi, sk_qi, qi)
    #     print('res=', crt_to_coef(res_dec_qi,qi))
    #     # self.assertEqual(crt_to_coef(res_dec_qi, qi), exp_qi)

    def test_he_mul_debug2(self):
        # qi = [999983, 999979]
        # pi = [1000003, 1000033]

        # qi = [991, 997]
        # pi = [1009,1013]

        # primes near 100 s.t. P = 1 mod 8.
        qi = [73, 89]
        pi = [97, 113]

        # qi = [97, 101]
        # pi = [103, 107]


        qipi = qi.copy()
        qipi.extend(pi)
        P = prod(pi)
        Q = prod(qi)

        print('qipi: ', qipi)
        print('Q=', Q)
        print('P=', P)
        print('QP=', Q*P)

        sk_qi = gen_sk(qi, debug=True)
        sk_qipi = mod_up(sk_qi, qi, pi)
        print("sk_qipi: ", sk_qipi)
        print("sk_qi: ", sk_qi)
        print("sk_coef: ", crt_to_coef(sk_qi,qi))
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
        ax1_qi, bx1_qi = encrypt(pt1_qi, sk_qi, qi, debug=True)
        ax2_qi, bx2_qi = encrypt(pt2_qi, sk_qi, qi, debug=True)
        print('bx1: ', crt_to_coef(bx1_qi, qi))
        print('bx2: ', crt_to_coef(bx2_qi, qi))
        print(crt_to_coef(decrypt(ax1_qi, bx1_qi, sk_qi, qi), qi))
        print(crt_to_coef(decrypt(ax2_qi, bx2_qi, sk_qi, qi),qi))

        relin_key_ax_qipi, relin_key_bx_qipi = gen_relin_key(sk_qipi, qipi, pi)

        res_ax_qi, res_bx_qi = he_mul(ax1_qi, bx1_qi, ax2_qi, bx2_qi, relin_key_ax_qipi, relin_key_bx_qipi, qi, pi,sk_qi,pt1_qi, pt2_qi)
        print(res_ax_qi,res_ax_qi)
        res_dec_qi = decrypt(res_ax_qi, res_bx_qi, sk_qi, qi)
        res_dec_coef = crt_to_coef(res_dec_qi,qi)
        print('res_dec_qi: ',res_dec_qi)
        print('res=', res_dec_coef)
        exp_coef = modulo(Polynomial(pt1_coef)*Polynomial(pt2_coef), Q)
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
