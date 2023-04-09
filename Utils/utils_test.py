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






    # def test_bit_comp_with_powers_of_2_different_dimensions(self):
    #     bgv = BGV(m_power=4, q=256, p=2, N=10)
    #     size = 8
    #     t = Polynomial([1, 2, 3])
    #     sk_before = [1, t*t, t**2]
    #     sk_after = [1, t]
    #
    #     plaintext1 = RingElement(Polynomial([1, 1]), bgv.m, 2)
    #     ctx1 = bgv.encrypt(plaintext1)
    #     ctx1_as_poly = [ctx1.c0.poly, ctx1.c1.poly]
    #     plaintext2 = RingElement(Polynomial([1, 0]), bgv.m, 2)
    #     ctx2 = bgv.encrypt(plaintext1)
    #     ctx2_as_poly = [ctx2.c0.poly, ctx2.c1.poly]
    #     ctx_multi = [ctx1_as_poly[0] * ctx2_as_poly[0],
    #                  ctx1_as_poly[0] * ctx2_as_poly[1] + ctx1_as_poly[1]*ctx2_as_poly[0],
    #                  ctx1_as_poly[1] * ctx2_as_poly[1]]
    #     ctx_multi_linearized = bit_decomp(ctx_multi, size=bgv.linearization_bit_size)
        # ctx_after =
        # b = [Polynomial([4, 5, 6]), Polynomial([8,8,8])]
        # result_bit_decomp = bit_decomp(a, size)
        # result_powers_of_2 = powers_of_2(b, size)
        # result = np.dot(result_bit_decomp, result_powers_of_2)
        # expected = np.dot(a, b)
        # self.assertEqual(expected, result)