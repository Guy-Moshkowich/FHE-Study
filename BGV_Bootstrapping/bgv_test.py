import unittest
from bgv import BGV
from numpy.polynomial import Polynomial


class TestBgv(unittest.TestCase):

    def setUp(self):
        p = 2
        q = 1000
        m = 16  # Phi_16(X)=x^8+1
        self.bgv = BGV(m, q, p)
        self.bgv.generate_key()

    def test_encrypt_decrypt(self):
        plaintext = Polynomial([1, 0, 1, 0, 1, 0, 1, 0])
        plaintext_decrypted = self.bgv.decrypt(self.bgv.encrypt(plaintext))
        print(plaintext)
        print(plaintext_decrypted)
        self.assertTrue(plaintext_decrypted == plaintext)

    def test_add(self):
        plaintext1 = Polynomial([1, 0, 1, 0, 1, 0, 1, 0])
        plaintext2 = Polynomial([0, 1, 0, 1, 0, 1, 0, 1])
        expected = Polynomial([1, 1, 1, 1, 1, 1, 1, 1])
        ctx1 = self.bgv.encrypt(plaintext1)
        ctx2 = self.bgv.encrypt(plaintext2)
        result = self.bgv.decrypt(ctx1 + ctx2)
        print(result)
        self.assertTrue(result, expected)

