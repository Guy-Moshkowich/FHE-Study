import unittest
from bgv import BGV
from numpy.polynomial import Polynomial
from ring_element import RingElement

class TestBgv(unittest.TestCase):


    def setUp(self):
        self.bgv = BGV(m_power=4, q=100, p=2, N=10)

    def test_encrypt_decrypt(self):
        plaintext = RingElement(Polynomial([1, 0, 1,0]), self.bgv.m, 2)
        plaintext_decrypted = self.bgv.decrypt(self.bgv.encrypt(plaintext))
        self.assertEqual(plaintext_decrypted,  plaintext)

    def test_add(self):
        plaintext1 = RingElement(Polynomial([1, 0, 1, 0, 1, 0, 1, 0]), self.bgv.m, 2)
        plaintext2 = RingElement(Polynomial([0, 1, 0, 1, 0, 1, 0, 1]), self.bgv.m, 2)
        expected = RingElement(Polynomial([1, 1, 1, 1, 1, 1, 1, 1]), self.bgv.m, 2)
        ctx1 = self.bgv.encrypt(plaintext1)
        ctx2 = self.bgv.encrypt(plaintext2)
        result = self.bgv.decrypt(ctx1 + ctx2)
        self.assertEqual(expected, result)

    def test_sub(self):
        plaintext1 = RingElement(Polynomial([1, 1, 1, 1, 0, 0, 0, 0]), self.bgv.m, 2)
        plaintext2 = RingElement(Polynomial([1, 1, 1, 1, 0, 0, 0, 0]), self.bgv.m, 2)
        expected = RingElement(Polynomial([0]), self.bgv.m, 2)
        ctx1 = self.bgv.encrypt(plaintext1)
        ctx2 = self.bgv.encrypt(plaintext2)
        result = self.bgv.decrypt(ctx1 - ctx2)
        self.assertEqual(expected, result)

    def test_pk(self):
        result = []
        for i in range(self.bgv.N):
            ring_element_q = self.bgv.pk_A[i][0] + self.bgv.secret_key[1]*self.bgv.pk_A[i][1]
            ring_element_2 = ring_element_q.change_modulo(2)
            result.append(ring_element_2)
            self.assertEqual(RingElement(Polynomial(0),self.bgv.m, 2), ring_element_2)
