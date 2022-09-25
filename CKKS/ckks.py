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
        c0 = (a * secret_key + plaintext + e).change_modulo(self.q)
        c1 = a
        return Ciphertext(c0, c1)

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


class Ciphertext:
    def __init__(self, c0, c1):
        self.c0 = c0
        self.c1 = c1
        self.q = self.c0.mod
        self.n = self.c0.m

    def switch_key(self, swk):
        minus_one = RingElement(Polynomial([-1]), self.n, self.q)
        new_ctx = Ciphertext(self.c0-(self.c1*swk.c0),  minus_one*self.c1*swk.c1)
        return new_ctx

    def decrypt(self, secret_key):
        return self.c0 - (self.c1*secret_key).change_modulo(self.q)

if __name__ == '__main__':
    main()