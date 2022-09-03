from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial
import numpy as np
class CKKS:

    def __init__(self, logn, q, max_added_noise=5):
        self.n = 2**logn
        self.q = q
        self.max_added_noise = max_added_noise
        self.secret_key = RingElement.random(self.n, self.q)
        self.k = 2

    def encrypt(self, plaintext: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.random(self.n, self.q, max_range=self.max_added_noise)
        ct0 = (a * self.secret_key + plaintext + e).change_modulo(self.q)
        ct1 = a
        return [ct0, ct1]

    def decrypt(self, ciphertext, secret_key):
        return (ciphertext[0] - ciphertext[1]*secret_key).change_modulo(self.q)

    def compose(self, ctx, k):  # compute f(X^k) mod (X^{m//2}+1, q)
        x_power_k_poly = [0]*(k+1)
        x_power_k_poly.extend(1)
        x_power_k = Polynomial(x_power_k_poly)
        return [ctx[0].compose(x_power_k),ctx[1].compose(x_power_k)]

    def switch_key(ct, new_secret_key):
        ct_new = ct
        return ct_new


if __name__ == '__main__':
    main()