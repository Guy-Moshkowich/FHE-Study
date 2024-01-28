import unittest
from context import Context
from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial
from Utils import utils


class TestContext(unittest.TestCase):

    def test_generate_secret_key(self):
        context = Context(log_n=10, q=1000)
        s = context.generate_secret_key()
        count_non_zero_coefs = 0
        for coef in s.poly.coef:
            if coef != 0:
                count_non_zero_coefs += 1
        self.assertLessEqual(count_non_zero_coefs, context.h)

    def test_generate_eval_key(self):
        context = Context(log_n=10, q=11, p=7)
        ek = context.eval_key
        sk = context.secret_key
        n = context.n
        q = context.q
        p = context.p
        p_ring = RingElement.const(p, n, q)
        dec = RingElement(ek.c0.poly-ek.c1.poly*sk.poly, n, q)
        m = p_ring*sk*sk
        diff = (dec - m).canonical_norm()
        self.assertLess(diff, 20)

    def test_generate_eval_key_ceil(self):
        context = Context(log_n=10, q=11, p=7)
        ek = context.eval_key
        sk = context.secret_key
        n = context.n
        q = context.q
        p = context.p
        inv_p = 1/p
        dec = RingElement(utils.ceil(inv_p*(ek.c0.poly - ek.c1.poly * sk.poly)), n, q)
        m = sk * sk
        diff = (dec - m).canonical_norm()
        self.assertLess(diff, 20)

    def test_encrypt_decrypt(self):
        context = Context(log_n=10, q=1000)
        plaintext = RingElement.random_ternary(context.n, context.q)
        ct = context.encrypt(plaintext)
        self.assert_equal(ct, context, plaintext, 20)

    def test_encrypt_decrypt_bad_secret_key(self):
        context = Context(log_n=10, q=1000)
        plaintext_expected = RingElement.random(context.n, context.q)
        ct = context.encrypt(plaintext_expected)
        bad_secret_key = RingElement.random(context.n, context.q)
        plaintext_actual = context.decrypt_core(ct, bad_secret_key)
        diff = plaintext_actual - plaintext_expected
        self.assertTrue(diff.canonical_norm() > 1000, diff.canonical_norm())

    def test_add(self):
        context = Context(log_n=10, q=1000)
        plaintext1 = RingElement.random(context.n, context.q)
        plaintext2 = RingElement.random(context.n, context.q)
        ctx1 = context.encrypt(plaintext1)
        ctx2 = context.encrypt(plaintext2)
        print(type(ctx1+ctx2))
        self.assert_equal(ctx1 + ctx2, context, plaintext1 + plaintext2, 30)

    def test_add_many(self):
        q = 1000
        context = Context(log_n=10, q=q)
        zero = RingElement(Polynomial([0]), context.n, context.q)
        ctx_acc = context.encrypt(zero)
        expected_plaintext = zero
        for i in range(10):
            plaintext = RingElement.random(context.n, context.q)
            ctx = context.encrypt(plaintext)
            ctx_acc += ctx
            expected_plaintext += plaintext
        self.assert_equal(ctx_acc, context, expected_plaintext, 50)

    def test_mul(self):
        context = Context(log_n=3, q=17, p=2, is_debug=True)
        plaintext1 = RingElement(Polynomial([2]), context.n, context.q)
        plaintext2 = RingElement(Polynomial([3]), context.n, context.q)
        ctx1 = context.encrypt(plaintext1)
        ctx2 = context.encrypt(plaintext2)
        self.assert_equal(ctx1 * ctx2, context, plaintext1 * plaintext2, 30)

    def test_swk_gen(self):
        context = Context(log_n=10, q=1009, p=1013)
        s1 = RingElement.random(context.n, context.q)
        s = context.secret_key
        swk_from_s_to_s1 = context.generate_swk(s, s1)
        plaintext_poly = swk_from_s_to_s1.c0.poly-swk_from_s_to_s1.c1.poly*s1.poly
        error = (RingElement(plaintext_poly-context.p*s.poly, context.n, context.q*context.p).canonical_norm())
        self.assertTrue(error < 20)

    def test_switch_key(self):
        context = Context(log_n=10, q=1009, p=1013)
        plaintext = RingElement.random(context.n, context.q)
        s1 = context.generate_secret_key()
        s = context.secret_key
        self.assertNotEqual(s1, s)
        swk_from_s_to_s1 = context.generate_swk(s, s1)
        # swk_from_s_to_s1_err = RingElement(swk_from_s_to_s1.c0.poly-swk_from_s_to_s1.c1.poly*s1.poly -context.p*s.poly,context.n, context.p*context.q)
        # print('swk_from_s_to_s1_err: ', swk_from_s_to_s1_err.canonical_norm())
        ct_wrt_s = context.encrypt(plaintext)
        # ct_wrt_s_err = ct_wrt_s.decrypt(s) - plaintext
        # print('ct_wrt_s_err: ', ct_wrt_s_err.canonical_norm())

        self.assert_almost_equal(context.decrypt_core(ct_wrt_s, s), plaintext, eps=20)
        ct_wrt_s1 = ct_wrt_s.switch_key(swk_from_s_to_s1, context.p)
        expected_c1_after_key_switch = RingElement(utils.ceil((-1/context.p)*swk_from_s_to_s1.c1.poly * ct_wrt_s.c1.poly), context.n, context.q)
        self.assert_almost_equal(ct_wrt_s1.c1, expected_c1_after_key_switch, eps=20)
        expected_c0_after_key_switch = RingElement(-utils.ceil((1/context.p)*(swk_from_s_to_s1.c1.poly * ct_wrt_s.c1.poly * s1.poly)) + plaintext.poly
                                                   , context.n, context.q)
        self.assert_almost_equal(ct_wrt_s1.c0, expected_c0_after_key_switch, eps=300)
        plaintext_result_poly = ct_wrt_s1.c0.poly - (ct_wrt_s1.c1.poly * s1.poly)
        plaintext_result = RingElement(plaintext_result_poly, context.n, context.q)
        diff = (plaintext_result - plaintext).canonical_norm()
        self.assertLess(diff, 400)




    def assert_equal(self, ctx, context, plaintext_expected, max_error):
        error = utils.get_canonical_error(ctx, context, plaintext_expected)
        self.assertTrue(error <= max_error, "actual diff " + str(error))

    def assert_almost_equal(self, elm1: RingElement, elm2: RingElement, eps):
        diff = RingElement(elm1.poly - elm2.poly, elm1.dim, elm1.mod).canonical_norm()
        self.assertLess(diff, eps, "error: "+str(diff))





    # def test_generate_swk(self):
    #     Context = Context(logn=10, q=1000)
    #     s = RingElement.random2(Context.n, Context.q)
    #     s_prime = RingElement.random2(Context.n, Context.q)
    #     swk = Context.generate_swk(s_prime, s)
    #     for i in range(math.ceil(math.log2(Context.q))):
    #         two_power_i = RingElement.const(2**i, Context.n, Context.q)
    #         self.assert_equal(swk[i], s_prime, two_power_i*s, 1000)

    # def test_switch_key_basic_for_binary_a(self):
    #     self.Context = Context(log_n=10, q=1000)
    #     plaintext = RingElement.random(self.Context.n, self.Context.q)
    #     s = self.Context.generate_secret_key()
    #     s_prime = self.Context.generate_secret_key()
    #     a_prime = RingElement.random(self.Context.n, self.Context.q)
    #     e_prime = RingElement.small_gauss(self.Context.n, self.Context.q)
    #     swk_from_s_prime_to_s = self.Context.generate_swk_core_bit_decomp(s_prime, s, a_prime, e_prime)
    #     self.assert_equal(swk_from_s_prime_to_s, s, s_prime, 20)
    #
    #     a = RingElement.random_binary(self.Context.n, self.Context.q)
    #     e = RingElement.small_gauss(self.Context.n, self.Context.q)
    #
    #     ct_wrt_s_prime = self.Context.encrypt_core(plaintext, a, s_prime, e)
    #     self.assert_equal(ct_wrt_s_prime, s_prime, plaintext, 20)
    #
    #     ct_wrt_s = ct_wrt_s_prime.switch_key_bit_decomp_basic(swk_from_s_prime_to_s)
    #     self.assert_equal(ct_wrt_s, s, plaintext, 150)

 # def test_binary_switch_key(self):
    #     Context = Context(log_n=10, q=100)
    #     plaintext = RingElement.random(Context.n, Context.q)
    #     s1 = RingElement.random(Context.n, Context.q)
    #     s = Context.secret_key
    #     swk_from_s_to_s1 = Context.generate_swk_bit_decomp(s, s1)
    #     self.assert_equal(swk_from_s_to_s1[0], s1, s, 20)
    #     ct_wrt_s = Context.encrypt(plaintext)
    #     ct_wrt_s1 = ct_wrt_s.switch_key_bit_decomp_basic(swk_from_s_to_s1)
    #     self.assert_equal(ct_wrt_s1, s1, plaintext, 200)
