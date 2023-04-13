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

    # def test_generate_fft_matrix(self):
    #     fft_rep, inv_fft_rep = generate_fft(8)
    #     numpy.testing.assert_almost_equal(fft_rep[0][0], 1, 0.0001)
    #     numpy.testing.assert_almost_equal(fft_rep[1][0], 1, 0.0001)
    #     numpy.testing.assert_almost_equal(fft_rep[2][0], 1, 0.0001)
    #     numpy.testing.assert_almost_equal(fft_rep[3][0], 1, 0.0001)
        #
        # numpy.testing.assert_almost_equal(fft_rep[0][1], 0.707 - 0.707j, 0.0001)
        # numpy.testing.assert_almost_equal(fft_rep[1][1], -0.707 - 0.707j, 0.0001)
        # numpy.testing.assert_almost_equal(fft_rep[2][1], -0.707 + 0.707j, 0.0001)
        # numpy.testing.assert_almost_equal(fft_rep[3][1], 0.707 + 0.707j, 0.0001)
        #
        # numpy.testing.assert_almost_equal(inv_fft_rep[0][0], 1, 0.0001)
        # numpy.testing.assert_almost_equal(inv_fft_rep[1][0], 1, 0.0001)
        # numpy.testing.assert_almost_equal(inv_fft_rep[2][0], 1, 0.0001)
        # numpy.testing.assert_almost_equal(inv_fft_rep[3][0], 1, 0.0001)

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
        numpy.testing.assert_almost_equal(encode(U,U_conj, dim, slots), expected, 0.001)

    def test_get_units(self):
        units = get_units(8)
        expected = [1, 3, 5, 7]
        self.assertTrue((units == expected))

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

    def test_fft_matrix(self):
        res = fft_matrix(8)
        self.assertEqual(4, len(res))
        for k in range(len(res)):
            self.assertEqual(4, len(res[k]))
        for k in range(len(res)):
            self.assertEqual(1, res[k][0])
        self.assertEqual((0.7071067811865476+0.7071067811865475j), res[0][1])
        self.assertEqual( (-0.7071067811865474+0.7071067811865477j), res[1][1])

    def test_inv_ntt(self):
        n = 8
        p = Polynomial([0, 1])  # p(x)=x
        res = inv_fft(n, fft(n, p))
        numpy.testing.assert_almost_equal((p-res).coef, [0]*4, 0.001)
