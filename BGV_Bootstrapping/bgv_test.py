import unittest
from bgv import BGV

from numpy.polynomial import Polynomial


class TestBgv(unittest.TestCase):

    def test_bgv_encrypt_decrypt(self):
        p = 2
        q = 1000
        m = 16  # Phi_16(X)=x^8+1
        bgv = BGV(m, q, p)
        bgv.generate_key()
        plaintext = Polynomial([1, 0, 1, 0, 1, 0, 1, 0])
        plaintext_decrypted = bgv.decrypt(bgv.encrypt(plaintext))
        print(plaintext)
        print(plaintext_decrypted)
        self.assertTrue(plaintext_decrypted == plaintext)
