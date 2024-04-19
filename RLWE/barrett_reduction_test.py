import unittest
from barrett_reduction import *
import random


class TestBarrett(unittest.TestCase):

    def test_modulo_number(self):
        for i in range(1000):
            n = random.randint(2, 2**25) # we do not use 2**32 in a 2**64 word size due to approximation error
            x = random.randint(1, 2**25 )
            barret = Barrett(n)
            self.assertEqual(barret.mod(x), x % n)
