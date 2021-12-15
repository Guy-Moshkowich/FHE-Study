import unittest
import bgv
from numpy.polynomial import Polynomial


class TestBgv(unittest.TestCase):

    def test_bgv_encrypt_decrypt(self):
        plaintext = Polynomial([1, 0, 1, 0, 1, 0, 1, 0])
        self.assertTrue(bgv.decrypt(bgv.encrypt(plaintext)) - plaintext) == Polynomial([0])
