import unittest
import CRT


class MyTestCase(unittest.TestCase):
    def test_crtToBigInt_small(self):
        primes = [3, 5, 7]
        crt = CRT.CRT(primes)
        some_int = 10
        crt2 = [some_int % primes[i] for i in range(crt.L)]
        self.assertEqual(some_int, crt.crtToBigInt(crt2, crt.L))

    def test_rescale_small(self):
        primes = [3, 5, 7]
        crt = CRT.CRT(primes)
        some_int = 15
        crt2 = [some_int % primes[i] for i in range(crt.L)]
        crt2_rescaled = crt.rescale(crt2)
        big_int = crt.crtToBigInt(crt2_rescaled, crt.L - 1)
        self.assertEqual(big_int, some_int//primes[crt.L - 1])

    def test_rescale_huge(self):
        primes = [9007199255019521, 9007199255347201, 18014398510645249, 18014398510661633]
        crt = CRT.CRT(primes)
        big_int = 10741460466601571728957078983049652939238945085500696422459756533
        expected = big_int//primes[crt.L-1]
        big_int_crt = [big_int % prime for prime in primes]
        print("big_int_crt=",big_int_crt)
        crt_rescaled = crt.rescale(big_int_crt)
        print("crt_rescaled=", crt_rescaled)

        big_int_rescaled = crt.crtToBigInt(crt_rescaled, crt.L-1)
        # poly_rescaled_elm0 = _coeffInCrt::[0] = [7241850199221906, 5367263381624476, 9464641570783496, 0, ]
        self.assertEqual(big_int_rescaled, expected)


if __name__ == '__main__':
    unittest.main()
