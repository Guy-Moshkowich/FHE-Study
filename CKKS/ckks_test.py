import unittest
from ckks import CKKS
from numpy.polynomial import Polynomial
from RLWE.ring_element import RingElement


class TestBgv(unittest.TestCase):

    def setUp(self):
        self.max_error = 10
        self.ckks = CKKS(logn=4, q=256, max_added_noise=self.max_error//2)

    def test_encrypt_decrypt(self):
        plaintext = RingElement(Polynomial([15, 123, 7, 99+256]), self.ckks.n, self.ckks.q)
        ct = self.ckks.encrypt(plaintext)
        plaintext_decrypted = self.ckks.decrypt(ct, self.ckks.secret_key)
        diff = plaintext_decrypted - plaintext
        self.assertTrue([coef < self.max_error for coef in diff.poly.coef])
