from CKKS.context import *
from Utils.utils import *
from RLWE.ring_element import RingElement

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

    # (d0, d1) + round(P^{-1}*d2*eval_key) mod Q
    # Where (d0,d1,d2) := (b1*b2, a1*b2+a2*b1, a1*a2)
    # eval_key = (a*sk + P*sk^2 + e, a)
    def __mul__(self, other):
        b1 = self.c0
        b2 = other.c0
        a1 = self.c1
        a2 = other.c1
        eval_key = self.context.eval_key
        n = self.context.n
        q = self.context.q
        p_inv = 1/self.context.p
        new_c0 = RingElement(ceil(p_inv*a1.poly*a2.poly*eval_key.c0.poly), n, q)
        new_c1 = RingElement(ceil(p_inv*a1.poly*a2.poly*eval_key.c1.poly), n, q)
        ctx1 = Ciphertext(self, b1*b2, a1*b2+a2*b1)
        ctx2 = Ciphertext(self, new_c0, new_c1)
        mul = ctx1 + ctx2
        return mul

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

    def rescale_by(self, p):
        p_inv = 1/p
        c0_new = Polynomial([round(p_inv*x) for x in self.c0.poly.coef])
        c1_new = Polynomial([round(p_inv*x) for x in self.c1.poly.coef])
        self.c0 = RingElement(c0_new, self.context.n,round(self.context.q*p_inv))
        self.c1 = RingElement(c1_new, self.context.n,round(self.context.q*p_inv))



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
