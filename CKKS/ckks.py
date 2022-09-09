from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial


class CKKS:

    def __init__(self, logn, q, max_added_noise=5):
        self.n = 2**logn
        self.q = q
        self.max_added_noise = max_added_noise
        self.secret_key = RingElement.random(self.n, self.q)
        self.k = 2
        self.max_error = 10#**(-5)

    def encrypt(self, plaintext: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.small_gauss(self.n, self.q)
        return self.encrypt_core(plaintext, a, self.secret_key, e)

    # ciphertext:= [a * secret_key + plaintext + e, a]
    def encrypt_core(self, plaintext:RingElement, a:RingElement, secret_key:RingElement, e:RingElement):
        ct0 = (a * secret_key + plaintext + e).change_modulo(self.q)
        ct1 = a
        ctx = [ct0, ct1]
        # assert (self.decrypt(ctx,secret_key)-plaintext).norm_canonical() <= 50, "fail encrypting.error is too big: " + str((self.decrypt(ctx,secret_key)-plaintext).norm_canonical())
        return ctx

    def decrypt(self, ciphertext, secret_key):
        return ciphertext[0] - (ciphertext[1]*secret_key).change_modulo(self.q)

    def rotate(self, ctx, k):  # compute f(X^k) mod (X^{m//2}+1, q)
        x_power_k_poly = [0]*(k+1)
        x_power_k_poly.extend(1)
        x_power_k = Polynomial(x_power_k_poly)
        print(x_power_k)
        ctx_compose = [ctx[0].compose(x_power_k), ctx[1].compose(x_power_k)]
        # TODO: key switch.

    def generate_swk(self, source_key: RingElement, target_key: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.small_gauss(self.n, self.q)
        return self.generate_swk_core(source_key, target_key, a, e)

    def generate_swk_core(self, source_key:RingElement, target_key:RingElement, a:RingElement, e:RingElement):
        return self.encrypt_core(plaintext=source_key, a=a, secret_key=target_key, e=e)

    def switch_key(self, ct, swk,s_prime,s, m):
        debug=  (ct[0] - ct[1]*s_prime-m).canonical_norm()
        assert  debug<20, debug # ct=[a_prime*s +s^\prime+e,a_prime]
        debug =  (swk[0] - swk[1]*s-s_prime).canonical_norm()
        assert  debug<20, debug # swk=[a*s^\prime+s+e,a]

        q = ct[0].mod
        n = ct[0].m
        minus_one = RingElement(Polynomial([-1]), n, q)
        self.assert_equal(ct, s_prime, m, 10)
        # a*s_prime +m+e - a*(a_prime*s_prime+s+e)=m+e-a*a_prime-a_prime*e
        a = ct[1]
        a_prime = swk[1]
        debug = ((a*s_prime + m -a*(a_prime*s + s_prime))-(m-a*a_prime*s)).canonical_norm()
        assert debug<10, debug

        debug = ((ct[0] - a * (a_prime * s + s_prime)) - (m - a * a_prime * s)).canonical_norm()
        assert debug < 10, debug

        ct_wrt_target = [ct[0]-(ct[1]*swk[0]),  minus_one* ct[1]*swk[1]]
        debug = (ct_wrt_target[0] +a*a_prime*s-m).canonical_norm()
        assert debug<100, debug


        debug = (ct_wrt_target[0] - ct_wrt_target[1]*s - m).canonical_norm()
        assert debug<100, debug
        return ct_wrt_target

    def assert_equal(self, ctx, sk, plaintext_expected, max_error):
        plaintext_actual = self.decrypt(ctx, sk)
        diff = plaintext_actual - plaintext_expected
        result = diff.canonical_norm() <= max_error
        assert result, "actual diff " + str(diff.canonical_norm())

if __name__ == '__main__':
    main()