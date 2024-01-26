import Utils.utils
from RLWE.ring_element import RingElement
from Utils.utils import *
from Utils.utils import ceil

class Context:
    h: int  # secret key Hamming weight

    def __init__(self, log_n: int, q: int, p: int = 1, max_added_noise=5):
        self.n = 2 ** log_n
        self.q = q
        self.p = p
        self.max_added_noise = max_added_noise
        self.k = 2
        self.max_error = 10#**(-5)
        self.h = 64
        self.secret_key = self.generate_secret_key()
        self.eval_key = self.generate_eval_key(self.secret_key)

    def generate_secret_key(self):
        t = Utils.utils.generate_ternary_polynomial(self.n//2, self.h)
        return RingElement(t, self.n, self.q)

    def generate_eval_key(self, sk):
        a = RingElement.random(self.n, self.q * self.p)
        e = RingElement.small_gauss(self.n, self.q * self.p)
        plaintext = self.p * sk.poly * sk.poly # Ps^2
        c0 = RingElement(a.poly * sk.poly + plaintext + e.poly, self.n, self.q * self.p)
        c1 = a
        return Ciphertext(self, c0, c1)

    def get_eval_key(self):
        return self.eval_key

    def encrypt(self, plaintext: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.const(0, self.n, self.q)
        while e == RingElement.const(0,self.n, self.q):
            e = RingElement.small_gauss(self.n, self.q)
        return self.encrypt_core(plaintext, a, self.secret_key, e)

    # ciphertext:= [a * secret_key + plaintext + e, a]
    def encrypt_core(self, plaintext: RingElement, a: RingElement, secret_key: RingElement, e: RingElement):
        c0 = a * secret_key + plaintext + e
        c1 = a
        return Ciphertext(self, c0, c1)

    """compute f(X^k) mod (X^{m//2}+1, q) """
    def rotate(self, ctx, k):
        x_power_k_poly = [0]*(k+1)
        x_power_k_poly.extend(1)
        x_power_k = Polynomial(x_power_k_poly)
        print(x_power_k)
        ctx_compose = [ctx[0].compose(x_power_k), ctx[1].compose(x_power_k)]
        # TODO: key switch.


    def generate_swk_core_bit_decomp(self, source_key:RingElement, target_key:RingElement, a:RingElement, e:RingElement):
        return self.encrypt_core(plaintext=source_key, a=a, secret_key=target_key, e=e)

    def generate_swk(self, source_key: RingElement, target_key: RingElement):
        a = RingElement.random(self.n, self.q * self.p)
        e = RingElement.small_gauss(self.n, self.q * self.p)
        plaintext = self.p * source_key.poly
        c0 = RingElement(a.poly * target_key.poly + plaintext + e.poly, self.n, self.q * self.p)
        c1 = a
        return Ciphertext(self, c0, c1)


class Ciphertext:
    c0: RingElement
    c1: RingElement
    context: Context


    def __init__(self, context, c0: RingElement, c1: RingElement):
        self.context = context
        self.c0 = c0
        self.c1 = c1

    def __add__(self, other):
        return Ciphertext(self.context, self.c0 + other.c0, self.c1 + other.c1)

    def __mul__(self, other):
        b1 = self.c0
        b2 = other.c0
        a1 = self.c1
        a2 = other.c1
        eval_key = self.context.eval_key
        n = self.context.n
        q = self.context.q
        p_inv = 1/self.context.p
        new_c0 = RingElement(-ceil(p_inv*a1.poly*a2.poly*eval_key.c0.poly), n, q)
        new_c1 = RingElement(-ceil(p_inv*a1.poly*a2.poly*eval_key.c1.poly), n, q)
        ctx1 = Ciphertext(self, b1*b2, a1*b2+a2*b1)
        ctx2 = Ciphertext(self, new_c0, new_c1)
        # ctx2 = Ciphertext(round(p_inv*a1*a2*eval_key.c0), round(p_inv*a1*a2*eval_key.c1))
        return ctx1 + ctx2

    def switch_key(self, swk, p):
        p_inv = 1/p
        new_c0 = RingElement(self.c0.poly - ceil(p_inv * self.c1.poly * swk.c0.poly), self.context.n, self.context.q)
        new_c1 = RingElement(-ceil(p_inv * self.c1.poly * swk.c1.poly), self.context.n, self.context.q)
        return Ciphertext(self, new_c0, new_c1)

    def switch_key_bit_decomp(self, swk: list["Ciphertext"]):
        c1_decomp = bit_decomp([self.c1.poly], math.ceil(math.log2(self.q)))
        zero = RingElement(Polynomial([0]), self.n, self.q)
        for i in range(len(c1_decomp)):
            tmp_ctx = Ciphertext(zero, RingElement(c1_decomp[i], self.n, self.q))
            tmp_ctx.switch_key_bit_decomp_basic(swk[i])
            if i == 0:
                new_ctx = tmp_ctx
            else:
                new_ctx += tmp_ctx
        new_ctx += Ciphertext(self.c0, zero)
        return new_ctx

    def switch_key_bit_decomp_basic(self, swk):
        minus_one = RingElement(Polynomial([-1]), self.n, self.q)
        new_ctx = Ciphertext(self.c0-(self.c1*swk.c0),  minus_one*self.c1*swk.c1)
        return new_ctx

    def decrypt(self, secret_key):
        return self.c0 - (self.c1*secret_key)




 # def generate_swk_bit_decomp(self, source_key: RingElement, target_key: RingElement):
    #     swk = []
    #     for i in range(math.ceil(math.log2(self.q))):
    #         a = RingElement.random(self.n, self.q)
    #         e = RingElement.small_gauss(self.n, self.q)
    #         two_power_i = RingElement.const(2**i, self.n, self.q)
    #         a_swk = self.generate_swk_core_bit_decomp(two_power_i * source_key, target_key, a, e)
    #         error = get_canonical_error(a_swk, target_key, two_power_i * source_key)
    #         assert error < 20, "error: " + str(error)
    #         swk.append(a_swk)
    #     return swk
