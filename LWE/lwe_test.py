import unittest
from lwe import *


class TestLwe(unittest.TestCase):
    eps = 10

    def test_gen_secret_key(self):
        params = LweParams(modulo=100, dim=5)
        sk = params.gen_secret_key()
        m = 50
        ct = params.encrypt(m, sk)
        self.assertTrue(params.decrypt(ct, sk)-m < self.eps)
