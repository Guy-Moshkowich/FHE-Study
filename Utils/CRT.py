class CRT:
    def __init__(self, primes):
        self.primes = primes
        self.L = len(primes)
        self.qi_inv_mod_qj = [[0] * self.L] * self.L
        self.qi_hat = [0] * self.L
        self.qi_hat_inv_mod_qi = [0] * self.L
        self.init_q_hat()
        self.init_qi_hat_inv_mod_qi()
        self.init_qi_inv_mod_qj()
        self.Q = 1
        for i in range(self.L):
            self.Q = self.Q * self.primes[i]




    def init_qi_hat_inv_mod_qi(self):
        for i in range(self.L):
            qi_hat_mod_qi = self.qi_hat[i] % self.primes[i]
            self.qi_hat_inv_mod_qi[i] = inv_mod(qi_hat_mod_qi, self.primes[i])


    def init_qi_inv_mod_qj(self):
        for i in range(self.L):
            for j in range(self.L):
                if i != j:
                    self.qi_inv_mod_qj[i][j] = inv_mod(self.primes[i], self.primes[j])


    def init_q_hat(self):
        self.qi_hat = [0] * self.L
        for i in range(self.L):
            prod = 1
            for j in range(self.L):
                if i != j:
                    prod = prod * self.primes[j]
            self.qi_hat[i] = prod

    def rescale(self, crt):
        res =[]
        for i in range(len(crt)):
            diff = (crt[i]-crt[self.L-1]) % self.primes[i]
            r = (diff*self.qi_inv_mod_qj[self.L-1][i]) % self.primes[i]
            res.append(r)
        return res

    def crtToBigInt(self, crt, level):
        Q = 1
        for i in range(level):
            Q = Q * self.primes[i]
        s = 0
        for i in range(level):
            s = s + (crt[i] * self.qi_hat[i] * self.qi_hat_inv_mod_qi[i] % Q)
        t = s % Q
        return t if t <= ((Q - 1) / 2) else t-Q


def inv_mod(qi, qj):
    return pow_mod(qi % qj, qj - 2, qj)


def pow_mod( a,  b,  m):
    a = a % m
    res = 1
    while b > 0:
        if b & 1:
            res = res*a %m
        a = (a*a) % m
        b >>= 1
    return res




def main():
    primes = [9007199255019521, 9007199255347201, 18014398510645249, 18014398510661633]
    crt_module = CRT(primes)
    big_int2 = 5711057810600191723137071220892716933009785801822074336086157765
    crt2 = [big_int2 % primes[i] for i in range(crt_module.L)]
    assert big_int2 == crt_module.crtToBigInt(crt2, crt_module.L)
    print(90071992553308165/18014398510661633)
    print(crt_module.Q)
    print(crt_module.Q/2>3723092598124141037845717590833633670269211418086548928700099434)

if __name__ == '__main__':
    main()

