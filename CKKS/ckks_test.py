import unittest
from ckks import CKKS
from RLWE.ring_element import RingElement


class TestCkks(unittest.TestCase):

    def setUp(self):
        self.max_error = 10
        self.ckks = CKKS(logn=3, q=256, max_added_noise=self.max_error//2)

    def test_encrypt_decrypt(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        ct = self.ckks.encrypt(plaintext)
        self.assert_equal(ct, self.ckks.secret_key, plaintext, 20)

    def test_encrypt_decrypt_bad_secret_key(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext_expected = RingElement.random(self.ckks.n, self.ckks.q)
        ct = self.ckks.encrypt(plaintext_expected)
        bad_secret_key = RingElement.random(self.ckks.n, self.ckks.q)
        plaintext_actual = ct.decrypt(bad_secret_key)
        diff = plaintext_actual - plaintext_expected
        self.assertTrue(diff.canonical_norm() > 1000, diff.canonical_norm())


    def test_encrypt_core_binary_a(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        a = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        sk1 = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        ct = self.ckks.encrypt_core(plaintext, a, sk1, e)
        self.assert_equal(ct, sk1, plaintext, 25)

    def test_generate_swk_core_with_binary_a(self):
        self.ckks = CKKS(logn=10, q=1000)
        a = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        s = RingElement.random(self.ckks.n, self.ckks.q)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q)
        swk = self.ckks.generate_swk_core(s_prime, s, a, e)
        self.assert_equal(swk, s, s_prime, 20)

    def test_generate_swk(self):
        s = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q, max_range=1)
        swk = self.ckks.generate_swk(s_prime, s)
        self.assert_equal(swk, s_prime, s, 100)


    def test_switch_key_for_binary_a(self):
        self.ckks = CKKS(logn=10, q=1000)
        plaintext = RingElement.random(self.ckks.n, self.ckks.q)
        s = RingElement.random(self.ckks.n, self.ckks.q)
        s_prime = RingElement.random(self.ckks.n, self.ckks.q)
        a_prime = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e_prime = RingElement.small_gauss(self.ckks.n, self.ckks.q)
        swk_from_s_prime_to_s = self.ckks.generate_swk_core(s_prime, s, a_prime, e_prime)
        self.assert_equal(swk_from_s_prime_to_s, s, s_prime, 20)

        a = RingElement.random_binary(self.ckks.n, self.ckks.q)
        e = RingElement.small_gauss(self.ckks.n, self.ckks.q)

        ct_wrt_s_prime = self.ckks.encrypt_core(plaintext, a, s_prime, e)
        self.assert_equal(ct_wrt_s_prime, s_prime, plaintext, 30)

        ct_wrt_s = ct_wrt_s_prime.switch_key(swk_from_s_prime_to_s)
        self.assert_equal(ct_wrt_s, s, plaintext, 150)

    def assert_equal(self, ctx, sk, plaintext_expected, max_error):
        plaintext_actual = ctx.decrypt(sk)
        diff = plaintext_actual - plaintext_expected
        result = diff.canonical_norm() <= max_error
        self.assertTrue(result, "actual diff " + str(diff.canonical_norm()))