import unittest
from ckks import CKKS
from RLWE.ring_element import RingElement
import math
from numpy.polynomial import Polynomial
from Utils import utils


class TestCkks(unittest.TestCase):

    def setUp(self):
        self.max_error = 10
        self.ckks = CKKS(log_n=3, q=256, max_added_noise=self.max_error // 2)

    def test_generate_secret_key(self):
        ckks = CKKS(log_n=10, q=1000)
        s = ckks.generate_secret_key()
        count_non_zero_coefs = 0
        for coef in s.poly.coef:
            if coef != 0:
                count_non_zero_coefs += 1
        self.assertLessEqual(count_non_zero_coefs, ckks.h)

    def test_encrypt_decrypt(self):
        ckks = CKKS(log_n=10, q=1000)
        plaintext = RingElement.random_ternary(ckks.n, ckks.q)
        ct = ckks.encrypt(plaintext)
        self.assert_equal(ct, ckks.secret_key, plaintext, 20)

    def test_encrypt_decrypt_bad_secret_key(self):
        self.ckks = CKKS(log_n=10, q=1000)
        plaintext_expected = RingElement.random(self.ckks.n, self.ckks.q)
        ct = self.ckks.encrypt(plaintext_expected)
        bad_secret_key = RingElement.random(self.ckks.n, self.ckks.q)
        plaintext_actual = ct.decrypt(bad_secret_key)
        diff = plaintext_actual - plaintext_expected
        self.assertTrue(diff.canonical_norm() > 1000, diff.canonical_norm())

    def test_encrypt_core_binary_a(self):
        self.ckks = CKKS(log_n=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        a = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        sk1 = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        ct = self.ckks.encrypt_core(plaintext, a, sk1, e)
        self.assert_equal(ct, sk1, plaintext, 25)

    def test_generate_swk_core_with_binary_a(self):
        self.ckks = CKKS(log_n=10, q=1000)
        a = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        s = RingElement.random(self.ckks.n, self.ckks.q)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q)
        swk = self.ckks.generate_swk_core_bit_decomp(s_prime, s, a, e)
        self.assert_equal(swk, s, s_prime, 20)

    # def test_generate_swk(self):
    #     ckks = CKKS(logn=10, q=1000)
    #     s = RingElement.random2(ckks.n, ckks.q)
    #     s_prime = RingElement.random2(ckks.n, ckks.q)
    #     swk = ckks.generate_swk(s_prime, s)
    #     for i in range(math.ceil(math.log2(ckks.q))):
    #         two_power_i = RingElement.const(2**i, ckks.n, ckks.q)
    #         self.assert_equal(swk[i], s_prime, two_power_i*s, 1000)

    # def test_switch_key_basic_for_binary_a(self):
    #     self.ckks = CKKS(log_n=10, q=1000)
    #     plaintext = RingElement.random(self.ckks.n, self.ckks.q)
    #     s = self.ckks.generate_secret_key()
    #     s_prime = self.ckks.generate_secret_key()
    #     a_prime = RingElement.random(self.ckks.n, self.ckks.q)
    #     e_prime = RingElement.small_gauss(self.ckks.n, self.ckks.q)
    #     swk_from_s_prime_to_s = self.ckks.generate_swk_core_bit_decomp(s_prime, s, a_prime, e_prime)
    #     self.assert_equal(swk_from_s_prime_to_s, s, s_prime, 20)
    #
    #     a = RingElement.random_binary(self.ckks.n, self.ckks.q)
    #     e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
    #
    #     ct_wrt_s_prime = self.ckks.encrypt_core(plaintext, a, s_prime, e)
    #     self.assert_equal(ct_wrt_s_prime, s_prime, plaintext, 20)
    #
    #     ct_wrt_s = ct_wrt_s_prime.switch_key_bit_decomp_basic(swk_from_s_prime_to_s)
    #     self.assert_equal(ct_wrt_s, s, plaintext, 150)

    def test_add(self):
        ckks = CKKS(log_n=10, q=1000)
        plaintext1 = RingElement.random(ckks.n, ckks.q)
        plaintext2 = RingElement.random(ckks.n, ckks.q)
        ctx1 = ckks.encrypt(plaintext1)
        ctx2 = ckks.encrypt(plaintext2)
        self.assert_equal(ctx1 + ctx2, ckks.secret_key, plaintext1 + plaintext2, 20)

    def test_add_many(self):
        q = 1000
        ckks = CKKS(log_n=10, q=q)
        zero = RingElement(Polynomial([0]), ckks.n, ckks.q)
        ctx_acc = ckks.encrypt(zero)
        expected_plaintext = zero
        for i in range(10):
            plaintext = RingElement.random(ckks.n, ckks.q)
            ctx = ckks.encrypt(plaintext)
            ctx_acc += ctx
            expected_plaintext += plaintext
        self.assert_equal(ctx_acc, ckks.secret_key, expected_plaintext, 50)

    # def test_binary_switch_key(self):
    #     ckks = CKKS(log_n=10, q=100)
    #     plaintext = RingElement.random(ckks.n, ckks.q)
    #     s1 = RingElement.random(ckks.n, ckks.q)
    #     s = ckks.secret_key
    #     swk_from_s_to_s1 = ckks.generate_swk_bit_decomp(s, s1)
    #     self.assert_equal(swk_from_s_to_s1[0], s1, s, 20)
    #     ct_wrt_s = ckks.encrypt(plaintext)
    #     ct_wrt_s1 = ct_wrt_s.switch_key_bit_decomp_basic(swk_from_s_to_s1)
    #     self.assert_equal(ct_wrt_s1, s1, plaintext, 200)

    def test_swk_gen(self):
        ckks = CKKS(log_n=10, q=1009, p=1013)
        s1 = RingElement.random(ckks.n, ckks.q)
        s = ckks.secret_key
        swk_from_s_to_s1 = ckks.generate_swk(s, s1)
        plaintext_poly = swk_from_s_to_s1.c0.poly-swk_from_s_to_s1.c1.poly*s1.poly
        error = (RingElement(plaintext_poly-ckks.p*s.poly, ckks.n, ckks.q*ckks.p).canonical_norm())
        self.assertTrue(error < 20)

    def test_switch_key(self):
        # ckks = CKKS(log_n=10, q=0x3fff4001, p=0x3ffe8001)
        ckks = CKKS(log_n=10, q=1009, p=7789)
        plaintext = RingElement.random(ckks.n, ckks.q)
        s1 = ckks.generate_secret_key()
        s = ckks.secret_key
        self.assertNotEqual(s1, s)
        swk_from_s_to_s1 = ckks.generate_swk(s, s1)
        ct_wrt_s = ckks.encrypt(plaintext)
        self.assert_almost_equal(ct_wrt_s.decrypt(s), plaintext, eps=20)
        ct_wrt_s1 = ct_wrt_s.switch_key(swk_from_s_to_s1, ckks.p)
        expected_c1_after_key_switch = RingElement(utils.ceil((-1/ckks.p)*swk_from_s_to_s1.c1.poly * ct_wrt_s.c1.poly), ckks.n, ckks.q)
        self.assert_almost_equal(ct_wrt_s1.c1, expected_c1_after_key_switch, eps=20)
        expected_c0_after_key_switch = RingElement(utils.ceil((-1/ckks.p)*(swk_from_s_to_s1.c1.poly * ct_wrt_s.c1.poly * s1.poly) + plaintext.poly)
                                                   , ckks.n, ckks.q)
        self.assert_almost_equal(ct_wrt_s1.c0, expected_c0_after_key_switch, eps=250)
        plaintext_result_poly = ct_wrt_s1.c0.poly - (ct_wrt_s1.c1.poly * s1.poly)
        plaintext_result = RingElement(plaintext_result_poly, ckks.n, ckks.q)
        diff = (plaintext_result - plaintext).canonical_norm()
        self.assertLess(diff, 400)

    def assert_equal(self, ctx, sk, plaintext_expected, max_error):
        error = utils.get_canonical_error(ctx, sk, plaintext_expected)
        self.assertTrue(error <= max_error, "actual diff " + str(error))

    def assert_almost_equal(self, elm1: RingElement, elm2: RingElement, eps):
        diff = RingElement(elm1.poly - elm2.poly, elm1.m, elm1.mod).canonical_norm()
        self.assertLess(diff, eps, "error: "+str(diff))
