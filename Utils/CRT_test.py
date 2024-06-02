import unittest
import CRT
import Utils.utils


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
        # primes in C++ 9007199255019521, 9007199255347201,18014398510645249, 18014398510661633
        primes = [9007199255019521, 9007199255347201, 18014398510645249, 18014398510661633]
        crt = CRT.CRT(primes)
        # bigInt C++: 10741460466601571728957078983049652939238945085500696422459756533
        big_int = 10741460466601571728957078983049652939238945085500696422459756533
        expected = big_int//primes[crt.L-1]
        print('expected=',expected)
        #c++  596270836367051106536170319143162274945881454247
        #here 596270836367051106536170319143162274945881454247

        big_int_crt = [big_int % prime for prime in primes]
        print("big_int_crt=", big_int_crt)
        # in C++: [3947265819819504, 2932785485109328, 368783843801005, 17638066271951182
        # here:   [3947265819819504, 2932785485109328, 368783843801005, 17638066271951182]
        print(3947265819819504-17638066271951182 + 2*9007199255019521)
        print(3947265819819504-17638066271951182 + 2*9007199255019521)
        crt_rescaled = crt.rescale(big_int_crt)
        print("big_int_rescaled_crt=", crt_rescaled)

        big_int_rescaled = crt.crtToBigInt(crt_rescaled, crt.L-1)
        # diff in C++: [4323597487479972, 3309116482336516, 745114891311712, 0, ]
        #              [4323598057907364, 3309117723852548, 745116082495072, 0, ]
        # diff here:   [4323598057907364, 3309117723852548, 745116082495072, 0]
        self.assertEqual(big_int_rescaled, expected)


    def test_decrypt(self):
        Utils.utils.read_uint64_from_file()

if __name__ == '__main__':
    unittest.main()


[18433053273257419938, 18432038792922709762, 745116082495072, 0, ]
[4323597487479972, 3309116482336516, 745114891311712, 0, ]
[4323598057907364, 3309117723852548, 745116082495072, 0]