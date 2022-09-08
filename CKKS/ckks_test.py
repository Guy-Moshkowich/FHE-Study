import unittest
from ckks import CKKS
from numpy.polynomial import Polynomial
from RLWE.ring_element import RingElement
import random

class TestCkks(unittest.TestCase):

    def setUp(self):
        self.max_error = 10
        self.ckks = CKKS(logn=3, q=256, max_added_noise=self.max_error//2)

    def test_encrypt_decrypt(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        ct = self.ckks.encrypt(plaintext)
        self.assert_equal(ct, self.ckks.secret_key, plaintext, 20)

        # plaintext_decrypted = self.ckks.decrypt(ct, self.ckks.secret_key)
        # diff = plaintext_decrypted - plaintext
        # self.assertTrue(diff.norm() <= self.max_error)

    def test_encrypt_core_binary_a(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        a = RingElement.random_binary(self.ckks.n, self.ckks.q, max_range=3)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        sk1 = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        ct = self.ckks.encrypt_core(plaintext, a, sk1, e)
        self.assert_equal(ct, sk1, plaintext, 20)

    def test_generate_swk_core_with_small_a(self):
        self.ckks = CKKS(logn=10, q=1000)
        a = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        s = RingElement.random(self.ckks.n, self.ckks.q)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q)
        swk = self.ckks.generate_swk_core(s_prime, s, a, e)
        self.assert_equal(swk, s_prime, s, 10)

    def test_generate_swk(self):
        s = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        swk = self.ckks.generate_swk(s_prime, s)
        self.assert_equal(swk, s_prime, s, 20)


        # plaintext_decrypted = self.ckks.decrypt(swk, s)
        # diff = plaintext_decrypted - s_prime
        # result = diff.norm() <= self.max_error
        # self.assertTrue(result, diff.norm())

    def test_switch_key_for_small_a(self):
        self.ckks = CKKS(logn=10, q=10000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        s = RingElement.random(self.ckks.n, self.ckks.q)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q)
        a_prime = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        e_prime = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        swk_from_s_prime_to_s = self.ckks.generate_swk_core(s_prime, s, a_prime, e_prime)
        # swk_from_s_prime_to_s = self.ckks.generate_swk(s_prime, s)
        self.assert_equal(swk_from_s_prime_to_s, s_prime, s, 20)

        a = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        e = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)

        print("e_prime=", e_prime)
        print("a=", a)
        print("e_prime*a=",e_prime*a)
        # assert (e_prime*a).norm() < self.max_error, "max|e'*a|="+str((e_prime*a).norm()) + ", avg|e'*a|="+str((e_prime*a).norm_avg())

        ct_wrt_s_prime = self.ckks.encrypt_core(plaintext, a, s_prime, e)
        self.assert_equal(ct_wrt_s_prime, s_prime, plaintext, 10)

        ct_wrt_s = self.ckks.switch_key(ct_wrt_s_prime, swk_from_s_prime_to_s)
        # self.assert_equal(ct_wrt_s, s, plaintext, 10)


    # def test_switch_key(self):
    #     for i in range(100):
    #         print('i: ', i)
    #         self.ckks = CKKS(logn=10, q=2**10, max_added_noise=6)
    #         plaintext = RingElement.random(self.ckks.n, self.ckks.q, max_range=self.ckks.q//10)
    #         s = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
    #         s_prime = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
    #         swk = self.ckks.generate_swk(s_prime, s)
    #         print('swk: ', swk[1].poly.coef[:10])
    #
    #         a = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
    #         # print('a: ',a)
    #         e = RingElement.random(self.ckks.n, self.ckks.q, max_range=5)
    #         ct = self.ckks.encrypt_core(plaintext, a, s_prime, e)
    #         print('ct: ', ct[1].poly.coef[:10])
    #         ct_wrt_s = self.ckks.switch_key(ct, swk)
    #         print('ct_wrt_s: ', ct_wrt_s[1].poly.coef[:10])
    #
    #         plaintext_decrypted = self.ckks.decrypt(ct_wrt_s, s)
    #         self.assert_poly_diff(plaintext_decrypted, plaintext, 10)

    def assert_equal(self, ctx, sk, plaintext_expected, max_error):
        plaintext_actual = self.ckks.decrypt(ctx, sk)
        diff = plaintext_actual - plaintext_expected
        result = diff.canonical_norm() <= max_error
        self.assertTrue(result, "actual diff " + str(diff.canonical_norm()))