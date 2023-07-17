import unittest
from utils import *
import numpy
import cmath

class TestUtils(unittest.TestCase):

    def test_dot_prod(self):
        result = dot_prod([1, 2, 3], [4, 5, 6])
        expected = 1 * 4 + 2 * 5 + 3 * 6
        self.assertEqual(expected, result)

    def test_bit_decomp_int(self):
        self.assertEqual([0, 1, 1], bit_decomp_int(z=6, size=3))

    def test_bit_decomp_int_size(self):
        self.assertEqual([0, 1, 1, 0, 0], bit_decomp_int(z=6, size=5))

    def test_bit_decomp_deg_0(self):
        result = bit_decomp_poly(Polynomial([6]), 3)
        expected = [Polynomial([0]), Polynomial([1]), Polynomial([1])]
        self.assertEqual(expected, result)

    def test_bit_decomp_poly(self):
        orig_poly = Polynomial([3,4,5])
        result = bit_decomp_poly(orig_poly, 3) #3+4x+5x^2
        reconstruct = Polynomial([0])
        for i in range(len(result)):
            for a in result[i].coef:
                self.assertTrue(a == 1 or a == 0)
            reconstruct += 2**i * result[i]
        self.assertEqual(reconstruct, orig_poly)


    def test_bit_decomp_poly_large_size(self):
        orig_poly = Polynomial([3,4,5])
        result = bit_decomp_poly(orig_poly, 10) #3+4x+5x^2
        reconstruct = Polynomial([0])
        for i in range(len(result)):
            reconstruct += 2**i * result[i]
        self.assertEqual(reconstruct, orig_poly)

    def test_bit_decomp(self):
        size = 8
        orig_poly_list = [Polynomial([3,4,5]), Polynomial([1,2,3])]
        result = bit_decomp(orig_poly_list, size)
        for j in range(len(orig_poly_list)):
            reconstruct = Polynomial([0])
            for l in range(size):
                reconstruct += 2**l * result[j + len(orig_poly_list)*l]
            self.assertEqual(reconstruct, orig_poly_list[j])

    def test_powers_of_2(self):
        size = 2
        orig_poly_list = [Polynomial([1, 2, 3]), Polynomial([4, 5, 6])]
        result = powers_of_2(orig_poly_list, size)
        expected = [Polynomial([1, 2, 3]), Polynomial([4, 5, 6]), Polynomial([2, 4, 6]), Polynomial([8, 10, 12])]
        self.assertEqual(expected, result)

    def test_bit_comp_with_powers_of_2(self):
        size = 8
        a = [Polynomial([1, 2, 3]), Polynomial([7,7,7])]
        b = [Polynomial([4, 5, 6]), Polynomial([8,8,8])]
        result_bit_decomp = bit_decomp(a, size)
        result_powers_of_2 = powers_of_2(b, size)
        result = np.dot(result_bit_decomp, result_powers_of_2)
        expected = np.dot(a, b)
        self.assertEqual(expected, result)

    def test_nth_primitive_roots_of_unity(self):
        roots = get_nth_primitive_roots_of_unity(4)
        numpy.testing.assert_almost_equal(roots[0].real, 0, 0.0001)
        numpy.testing.assert_almost_equal(roots[0].imag, 1, 0.0001)
        numpy.testing.assert_almost_equal(roots[1].real , 0, 0.0001)
        numpy.testing.assert_almost_equal(roots[1].imag, -1, 0.0001)

    def test_conjugation_of_nth_primitive_roots_of_unity(self):
        m = 8
        roots = get_nth_primitive_roots_of_unity(m)
        for i in range(len(roots)):
            found_conju_pair = False;
            for j in range(len(roots)):
                if cmath.isclose(roots[i].conjugate(), roots[j]):
                    found_conju_pair = True;
            self.assertTrue(found_conju_pair)

    def test_recenter(self):
        self.assertTrue(recenter(5, 1000) == 5)
        self.assertTrue(recenter(999, 1000) == -1, 'recenter: ' + str(recenter(999, 1000)))

    def test_modulo_polynomial(self):
        f = Polynomial([1, 2, 0, 3, 4])
        g = Polynomial([2, 1, 1])
        res = modulo_polynomial(f, g)
        self.assertTrue((res.coef == [15, 11]).all())

    def test_modulo_polynomial(self):
        dim = 16
        poly = Polynomial([0,0,0,0,0,0,0,0,1])  # poly=X**8
        phi_m = build_cyclotomic_poly(dim)  # phi(16)=X**8+1
        poly = modulo_polynomial(poly, phi_m)  # X**8 mod X**8+1 = -1
        self.assertTrue((poly.coef == [-1]).all())

    def test_modulo_coeff(self):
        f = Polynomial([1, 2, 0, 3, 4])
        res = modulo_int(f, 2)
        self.assertTrue((res.coef == [1, 0, 0, 1, 0]).all())


    def test_decode(self):
        U,U_conj = generate_canonical(8)
        coeffs_rep = [10/4, math.sqrt(2), 10/4, math.sqrt(2)/2]
        expected = [3+4j, 2-1j, 3-4j, 2+1j]
        numpy.testing.assert_almost_equal(decode(U,U_conj, coeffs_rep), expected, 0.001)

    def test_encode(self):
        dim = 8
        U, U_conj = generate_canonical(dim)
        slots = [3+4j, 2-1j]
        expected = [10/4, math.sqrt(2), 10/4, math.sqrt(2)/2]
        actual = encode(U,U_conj, dim, slots)
        numpy.testing.assert_almost_equal(expected, actual, 0.001)

    def test_get_units(self):
        actual = get_units(8)
        expected = [1, 3, 5, 7]
        self.assertEqual(expected, actual)

    def test_order(self):
        self.assertEqual(1, order(16, 1))
        self.assertEqual(4, order(16, 3))  # 3^2=9, 3^3=27=11, 3^4 =33=1
        self.assertEqual(4, order(16, 5))  # 5^2=9, 5^3=27=11, 5^4 =1
        self.assertEqual(2, order(16, 7))  # 7^2=49=1
        self.assertEqual(2, order(16, 9))  # 9^2=81
        self.assertEqual(4, order(16, 11))  # 11^2=121=9, 11^3=99=3, 11^4=33=1
        self.assertEqual(4, order(16, 13))  # 13^2=169=9,
        self.assertEqual(2, order(16, 15))

    def test_generator(self):
        self.assertTrue(units_cyc_gen(8) in [3, 5, 11, 13])

    def test_fft(self):
        n = 8
        p = Polynomial([0, 1])  # p(x)=x
        res = fft(n, p)
        self.assertEqual(4, len(res))
        self.assertEqual(res[0], np.conj(res[2]))
        self.assertEqual(res[1], np.conj(res[3]))
        self.assertEqual(res, [(0.7071067811865476+0.7071067811865475j), (-0.7071067811865474+0.7071067811865477j), (0.7071067811865476-0.7071067811865475j), (-0.7071067811865474-0.7071067811865477j)])

    def test_fft_multiplication(self):
        n = 8
        a = Polynomial([0, 1])  # a(x)=x
        b = Polynomial([0, 0, 1]) # b(x)=x^2
        c = a*b  # c(x)=x^3 ==> [0,0,0,1]
        a_fft = fft(n, a)
        b_fft = fft(n, b)
        expected = fft(n, c)
        actual = [a_fft[i]*b_fft[i]  for i in range(n//2)]
        self.assertEqual(actual, expected)

    def test_fft_multiplication_negacyclic(self):
        n = 8
        a = Polynomial([1, 1])  # a(x)=x+1
        b = Polynomial([1, 0, 0, 1]) # b(x)=x^3+1
        c = Polynomial([0,1,0,1])  # c(x)=x^4+x^3+x+1=x^3+x ==> [0,1,0,1]
        a_fft = fft(n, a)
        b_fft = fft(n, b)
        expected = fft(n, c)
        actual = [a_fft[i]*b_fft[i]  for i in range(n//2)]
        np.testing.assert_almost_equal(expected, actual, 0.0001)

    def test_fft_matrix(self):
        res = fft_matrix(8)
        self.assertEqual(4, len(res))
        for k in range(len(res)):
            self.assertEqual(4, len(res[k]))
        for k in range(len(res)):
            self.assertEqual(1, res[k][0])
        self.assertEqual((0.7071067811865476+0.7071067811865475j), res[0][1])
        self.assertEqual( (-0.7071067811865474+0.7071067811865477j), res[1][1])

    def test_fft_and_inv_ftt(self):
        n = 8
        p = Polynomial([0, 1])  # p(x)=x
        res = inv_fft(n, fft(n, p))
        numpy.testing.assert_almost_equal((p-res).coef, [0]*4, 0.001)

    def test_is_power_2(self):
        self.assertTrue(is_power_2(8))
        self.assertFalse(is_power_2(7))

    def test_generate_primitive_root_of_unity(self):
        n = 4
        p = 13
        x = generate_primitive_root_of_unity(n, p)
        for k in range (1, n):
            self.assertNotEqual(1, (x ** k) % p, "x="+str(x)+", k="+str(k))
        self.assertEqual(1, (x**n) % p)

    def test_ntt_multiplication(self):
        # precondition: p mod 2n =1, 17 mod 8 =1
        n = 4
        p = 17
        a = [1,1,0,0] # x+1
        b = [0,1,0,0] # x
        expected = ntt(n, p, [0, 1, 1, 0]) # x*(x+1)=x^2+x => [0,1,1,0]
        ntt_a = ntt(n, p, a)
        ntt_b = ntt(n, p, b)
        actual = [ntt_a[i]*ntt_b[i] % p for i in range(len(a))]
        self.assertEqual(expected, actual)

    def test_ntt_multiplication_negacyclic(self):
        n = 4 # phi(x)=x^4+1
        p = 17
        a = [1,1,0,0] # x+1
        b = [1,0,0,1] # x^3+1
        expected = ntt(n, p, [0, 1, 0, 1]) # (x+1)(x^3+1)=x^4+x^3+x+1=x^3+x
        ntt_a = ntt(n, p, a)
        ntt_b = ntt(n, p, b)
        actual = [ntt_a[i]*ntt_b[i] % p for i in range(len(a))]
        self.assertEqual(expected, actual)

    def test_ntt_mat_mul_inv_ntt_mat_is_identity_mat(self):
        n = 4
        p = 17
        inv_ntt = inv_ntt_matrix(n, p)
        ntt = ntt_matrix(n, p)
        actual = (np.array(inv_ntt).dot(np.array(ntt)))
        for i in range(len(actual)):
            r = actual[i]
            for j in range(len(r)):
                if j==i:
                    self.assertEqual(1, r[j]%p)
                else:
                    self.assertEqual(0, r[j]%p)


    def test_ntt_and_inv_ntt(self):
        n = 4
        p = 17
        a = [1, 0, 0, 0]
        ntt_a = ntt(n, p, a)
        print('ntt_a: ', ntt_a)
        actual = inv_ntt(n, p, ntt_a)
        print('inv_ntt: ', actual)
        expected = a
        np.testing.assert_almost_equal(expected, actual, 0.0001)

    def test_ntt_and_inv_ntt_centering_modulo(self):
        n = 4
        p = 73
        a = [14.0, 3.0, 0.0, -3.0]
        ntt_a = ntt(n, p, a)
        print('ntt_a: ', ntt_a)
        actual = inv_ntt(n, p, ntt_a)
        print('inv_ntt: ', actual)
        expected = a
        np.testing.assert_almost_equal(expected, actual, 0.0001)

    def test_special_fft_matrix(self):
        actual = special_fft_matrix(n=4)
        expected = [[(1+0j), (0.7+0.7j), (0+1j), (-0.7+0.7j)], [(1+0j), (-0.7+0.7j), (0-1j), (0.7+0.7j)], [(1+0j), (-0.7-0.7j), (0+1j), (0.7-0.7j)], [(1+0j), (0.7-0.7j), (0-1j), (-0.7-0.7j)]]
        np.testing.assert_almost_equal(expected, actual, 0.0001)


    def test_special_fft_multiplication(self):
        n = 4
        a = Polynomial([0, 1])  # a(x)=x
        b = Polynomial([0, 0, 1]) # b(x)=x^2
        c = a*b  # c(x)=x^3 ==> [0,0,0,1]
        a_sfft = special_fft(n, a)
        b_sfft = special_fft(n, b)
        expected = special_fft(n, c)
        actual = [a_sfft[i]*b_sfft[i]  for i in range(n)]
        np.testing.assert_almost_equal(expected, actual, 0.0001)

    def test_special_fft_multiplication_negacyclic(self):
        n = 4
        a = Polynomial([1, 1])  # a(x)=x+1
        b = Polynomial([1, 0, 0, 1]) # b(x)=x^3+1
        c = Polynomial([0,1,0,1])  # c(x)=x^4+x^3+x+1=x^3+x ==> [0,1,0,1]
        a_sfft = special_fft(n, a)
        b_sfft = special_fft(n, b)
        expected = special_fft(n, c)
        actual = [a_sfft[i]*b_sfft[i]  for i in range(n)]
        np.testing.assert_almost_equal(expected, actual, 0.0001)

    def test_special_fft_and_inv_special_ftt(self):
        n = 4
        poly = Polynomial([0, 1, 0, 0])  # p(x)=x
        actual = inv_special_fft(n, special_fft(n, poly))
        numpy.testing.assert_almost_equal(actual.coef, poly.coef, 0.001)

    def test_special_inv_fft_real_poly(self):
        n = 4
        actual = inv_special_fft(n, [1,2,2,1])
        for c in actual.coef:
            if math.fabs(c.imag) > 0.000001:
                self.fail('c.imag=', c.imag)

    def test_inv(self):
        actual = inv(x=4, p=17)
        expected = 13
        self.assertEqual(expected, actual)

    def test_bit_reverse(self):
        self.assertEqual(bit_reverse(0, 1), 0)  # 0 should remain unchanged (1 bit)
        self.assertEqual(bit_reverse(1, 1), 1)  # 1 should remain unchanged (1 bit)
        self.assertEqual(bit_reverse(10, 4), 5)  # Binary: 0000 1010 -> 0000 0101 (4 bits)
        self.assertEqual(bit_reverse(255, 8), 255)  # Binary: 1111 1111 -> 1111 1111 (8 bits)
        self.assertEqual(bit_reverse(65535, 16), 65535)  # Binary: 1111 1111 1111 1111 -> 1111 1111 1111 1111 (16 bits)

    # def test_inplace(self):
    #     n = 4
    #     a = Polynomial([1, 1, 0, 0])  # a(x)=x+1
    #     b = Polynomial([1, 0, 0, 1])  # b(x)=x^3+1
    #     c = Polynomial([0, 1, 0, 1])  # c(x)=x^4+x^3+x+1=x^3+x ==> [0,1,0,1]
    #     expected = special_fft(n, a)
    #     inplaceNegacyclicNTT(n, [0,0,1,1])
    #     actual = a.coef #inplace
    #     np.testing.assert_almost_equal(actual, expected, 0.0001)
