import unittest
from bgv import BGV
from numpy.polynomial import Polynomial


class TestBgv(unittest.TestCase):

    def setUp(self):
        # Phi_16(X)=x^8+1
        self.bgv = BGV(m_power=2, q=10, p=2, N=1)

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
        self.assertTrue(expected, result)

    def test_dot_prod(self):
        result = self.bgv.dot_prod([1,2,3],[4,5,6])
        expected = 1*4+2*5+3*6
        self.assertEqual(expected, result)


    def test_modulo(self):
        import numpy as np
        p1 = Polynomial([1,2,1,1])
        p2 = Polynomial([1,0,1])

        print(p1)
        print(np.array(p1.coef))
        poly_modulo_phim = np.polydiv(np.flip(p1.coef), np.flip(p2.coef))
        p3 = Polynomial(np.flip(poly_modulo_phim[0]))
        print(p3)
        p4 = Polynomial(np.flip(poly_modulo_phim[1]))
        print('p4: ', p4)

        # poly = Polynomial([1,2,11,12])
        # poly1 = Polynomial([1,2])
        # poly2 = Polynomial([11,12])
        #
        #
        # poly_mod_q = self.bgv.modulo(poly1*poly2, self.bgv.q)
        #
        # print(poly_mod_q)

    def test_pk(self):
        result = []
        for i in range(self.bgv.N):
            poly = self.bgv.pk_A[i][0] + self.bgv.secret_key[1]*self.bgv.pk_A[i][1]
            poly_mod_q = self.bgv.modulo(poly, self.bgv.q)
            poly_mod_q_mod_p = self.bgv.modulo(poly_mod_q, self.bgv.p)
            result.append(poly_mod_q_mod_p)

        # A_trans = [[self.bgv.pk_A[j][i] for j in range(len(self.bgv.pk_A))] for i in range(len(self.bgv.pk_A[0]))]
        # result = [self.bgv.modulo(x+y*self.bgv.secret_key[1], self.bgv.q) for x,y in zip(A_trans[0], A_trans[1])]
        expected = [Polynomial(0) for i in range(self.bgv.N)]
        print(result)
        # print(A_trans[0])
        # print(A_trans[1])
        self.assertEqual(expected, result)