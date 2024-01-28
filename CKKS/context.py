import Utils.utils
from RLWE.ring_element import RingElement
from ciphertext import *


class Context:
    h: int  # secret key Hamming weight

    def __init__(self, log_n: int, q: int, p: int = 1, max_added_noise=5, is_debug=False):
        self.n = 2 ** log_n
        self.q = q
        self.p = p
        self.max_added_noise = max_added_noise
        self.is_debug = is_debug
        self.k = 2
        self.max_error = 10#**(-5)
        self.h = 64
        self.secret_key = self.generate_secret_key()
        self.eval_key = self.generate_eval_key(self.secret_key)

    def generate_secret_key(self):
        t = Utils.utils.generate_ternary_polynomial(self.n//2, self.h)
        return RingElement(t, self.n, self.q)

    # eval_key = (a*sk + P*sk^2 + e , a)
    def generate_eval_key(self, sk):
        a = RingElement.random(self.n, self.q * self.p)
        if not self.is_debug:
            e = RingElement.small_gauss(self.n, self.q * self.p)
        else:
            e = RingElement.const(0, self.n, self.q * self.p)
        plaintext = self.p * sk.poly * sk.poly # Ps^2
        c0 = RingElement(a.poly * sk.poly + plaintext + e.poly, self.n, self.q * self.p)
        c1 = a
        return Ciphertext(self, c0, c1)

    def get_eval_key(self):
        return self.eval_key

    def encrypt(self, plaintext: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.const(0, self.n, self.q)
        if not self.is_debug:
            while e == RingElement.const(0,self.n, self.q):
                e = RingElement.small_gauss(self.n, self.q)
        return self.encrypt_core(plaintext, a, self.secret_key, e)

    # ciphertext:= [a * secret_key + plaintext + e, a]
    def encrypt_core(self, plaintext: RingElement, a: RingElement, secret_key: RingElement, e: RingElement):
        c0 = a * secret_key + plaintext + e
        c1 = a
        return Ciphertext(self, c0, c1)

    def decrypt(self, ciphertext) :
        return self.decrypt_core(ciphertext, self.secret_key)

    def decrypt_core(self, ciphertext, sk):
        return ciphertext.c0 - (ciphertext.c1 * sk)



    """compute f(X^k) mod (X^{m//2}+1, q) """
    def rotate(self, ctx, k):
        x_power_k_poly = [0]*(k+1)
        x_power_k_poly.extend(1)
        x_power_k = Polynomial(x_power_k_poly)
        print(x_power_k)
        ctx_compose = [ctx[0].compose(x_power_k), ctx[1].compose(x_power_k)]
        # TODO: key switch.

    #
    # def generate_swk_core_bit_decomp(self, source_key:RingElement, target_key:RingElement, a:RingElement, e:RingElement):
    #     return self.encrypt_core(plaintext=source_key, a=a, secret_key=target_key, e=e)

    def generate_swk(self, source_key: RingElement, target_key: RingElement):
        a = RingElement.random(self.n, self.q * self.p)
        e = RingElement.small_gauss(self.n, self.q * self.p)
        plaintext = self.p * source_key.poly
        c0 = RingElement(a.poly * target_key.poly + plaintext + e.poly, self.n, self.q * self.p)
        c1 = a
        return Ciphertext(self, c0, c1)

