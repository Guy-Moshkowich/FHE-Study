from RLWE.ring_element import RingElement

class CKKS:

    def __init__(self, logn, q, max_added_noise=5):
        self.n = 2**logn
        self.q = q
        self.max_added_noise = max_added_noise
        self.secret_key = RingElement.random(self.n, self.q)

    def encrypt(self, plaintext: RingElement):
        a = RingElement.random(self.n, self.q)
        e = RingElement.random(self.n, self.q, max_range=self.max_added_noise)
        ct0 = (a * self.secret_key + plaintext + e).change_modulo(self.q)
        ct1 = a
        return [ct0, ct1]

    def decrypt(self, ciphertext, secret_key):
        return (ciphertext[0] - ciphertext[1]*secret_key).change_modulo(self.q)


    def switch_key(ct, new_secret_key):
        ct_new = ct
        return ct_new


if __name__ == '__main__':
    main()