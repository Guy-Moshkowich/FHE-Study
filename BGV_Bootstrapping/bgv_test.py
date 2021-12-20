import unittest
from bgv import BGV
from numpy.polynomial import Polynomial


class TestBgv(unittest.TestCase):

    def setUp(self):
        self.bgv = BGV(m_power=4, q=100, p=2, N=10)

    def test_encrypt_decrypt(self):
        plaintext = Polynomial([1, 0, 1,0])
        plaintext_decrypted = self.bgv.decrypt(self.bgv.encrypt(plaintext))
        print('plaintext:', plaintext)
        print('plaintext_decrypted: ', plaintext_decrypted)
        self.assertTrue(all(v == 0 for v in plaintext_decrypted - plaintext))

    def test_add(self):
        plaintext1 = Polynomial([1, 0, 1, 0, 1, 0, 1, 0])
        plaintext2 = Polynomial([0, 1, 0, 1, 0, 1, 0, 1])
        expected = Polynomial([1, 1, 1, 1, 1, 1, 1, 1])
        ctx1 = self.bgv.encrypt(plaintext1)
        ctx2 = self.bgv.encrypt(plaintext2)
        result = self.bgv.decrypt(ctx1 + ctx2)
        print(result)
        self.assertTrue(expected, result)

    def test_dot_prod(self):
        result = self.bgv.dot_prod([1,2,3],[4,5,6])
        expected = 1*4+2*5+3*6
        self.assertEqual(expected, result)

    def test_pk(self):
        result = []
        for i in range(self.bgv.N):
            poly = self.bgv.pk_A[i][0] + self.bgv.secret_key[1]*self.bgv.pk_A[i][1]
            poly_mod_q = self.bgv.modulo(poly, self.bgv.q)
            poly_mod_q_mod_p = self.bgv.modulo(poly_mod_q, self.bgv.p)
            result.append(poly_mod_q_mod_p)
            print(result)
            self.assertTrue(all(v == 0 for v in poly_mod_q_mod_p.coef))
