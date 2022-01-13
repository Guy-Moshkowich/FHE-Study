from ring_element import RingElement
from utils import *
import math


class BGV:
    epsilon = 5

    def __init__(self, m_power, q, p, N):
        self.m = 2**m_power # Phi_m(X)=x^m+1
        self.q = q
        self.p = p
        self.N = N
        self.int_2 = RingElement(Polynomial(2), self.m, self.q)
        self.int_1 = RingElement(Polynomial(1), self.m, self.q)
        self.int_0 = RingElement(Polynomial(0), self.m, self.q)
        t = RingElement.random(m=self.m, mod=self.q)
        self.secret_key = (self.int_1, t)
        self.public_key = self.generate_public_key(N)
        self.linearization_bit_size = int(math.log2(self.q)) + 1
        self.linearization_matrix = self.linearization_matrix()
        self.modulo_chain = [17,13,11]

    def linearization_matrix(self):
        t = self.secret_key[1].poly
        s1 = [Polynomial(1), t, t*t] #secret key tensor
        B = self.generate_public_key(len(s1) * self.linearization_bit_size)
        s1_powers_of_2_poly = powers_of_2(s1, self.linearization_bit_size)
        s1_powers_of_2 = [RingElement(p, self.m, self.q) for p in s1_powers_of_2_poly]
        A = np.array([B[0] + s1_powers_of_2, B[1]]).transpose()
        return A

    def generate_public_key(self, size):
        A = []
        for i in range(size):
            b = RingElement.random(self.m, self.q)
            e = RingElement.random(self.m, self.q, max_range=self.epsilon)
            A.append(np.array([b*self.secret_key[1]+self.int_2*e, self.int_0 - b]))
        return np.array(A).transpose()

    def encrypt(self, plaintext: RingElement):
        plaintext_q = plaintext.change_modulo(self.q)
        r = np.array([RingElement.random(self.m, self.q) for i in range(self.N)])
        ctx = np.matmul(self.public_key, r) + [plaintext_q, self.int_0]
        return Ciphertext(ctx[0], ctx[1], self)

    def decrypt(self, ciphertext, secret_key):
        noise = secret_key[1] * ciphertext.c1
        return (ciphertext.c0 + noise).change_modulo(2)


class Ciphertext:

    def __init__(self, c0: RingElement, c1: RingElement, bgv: BGV):
        self.c0 = c0
        self.c1 = c1
        self.bgv = bgv
        self.modulo_chain_index = 0

    def __add__(self, other):
        return Ciphertext(self.c0 + other.c0, self.c1 + other.c1, self.bgv)

    def __sub__(self, other):
        return Ciphertext(self.c0 - other.c0, self.c1 - other.c1, self.bgv)

    def __mul__(self, other):
        multi_ctx_3_params_poly = [self.c0.poly * other.c0.poly,
                                   self.c0.poly * other.c1.poly + self.c1.poly * other.c0.poly,
                                   self.c1.poly * other.c1.poly]
        ctx_bit_decomp_poly = bit_decomp(multi_ctx_3_params_poly,
                                         size=self.bgv.linearization_bit_size)
        ctx_bit_decomp = np.array([[RingElement(p, self.bgv.m, self.bgv.q)] for p in ctx_bit_decomp_poly])
        new_ciphertext = np.matmul(ctx_bit_decomp.transpose(), self.bgv.linearization_matrix)
        return Ciphertext(new_ciphertext[0][0], new_ciphertext[0][1], self.bgv)

    def __str__(self):
        return 'c0: ' + str(self.c0) + ', c1: ' + str(self.c1)

    def scale(self):
        next_modulo_chain_index = self.modulo_chain_index + 1
        assert next_modulo_chain_index < len(self.bgv.modulo_chain)
        q_next = self.bgv.modulo_chain[next_modulo_chain_index]
        q_current = self.bgv.modulo_chain[self.modulo_chain_index]
        self.bgv.scale([self.c0, self.c1], q_next,q_current)
        self.modulo_chain_index = next_modulo_chain_index




