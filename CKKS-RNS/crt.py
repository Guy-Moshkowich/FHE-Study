from numpy.polynomial import Polynomial
from Utils import utils
import random

random.seed(0)
n = 8
qi = [97, 193]
pi = [101,103]
qipi = [97,193,101,103]
Q = 97*193
P = 101*103
cyc = Polynomial([1,0,0,0,0,0,0,0,1])


def gen_sk(primes):
    M = 1
    for p in primes:
        M = M * p
    return [random.choice([M-1, 0, 1]) for _ in range(n)]


def coef_to_crt(poly_coef, primes):
    crt = [0]*(n*len(primes))
    for j in range(len(primes)):
        for i in range(len(poly_coef)):
            idx = j*n+i
            crt[idx] = poly_coef[i] % primes[j]
    return crt


def crt_to_coef(poly_crt, primes):
    out = []
    for i in range(n):
        elm = []
        for j in range(len(primes)):
            idx = i+j*n
            elm.append(poly_crt[idx])
        out.append(crt_to_coef_elm(elm, primes))
    return out


def crt_to_coef_elm(crt_elm, primes):
    M = 1
    for p in primes:
        M = M * p
    out = 0
    for j in range(len(crt_elm)):
        out = out + crt_elm[j]*qi_hat(j,primes)*inv_qi_hat(j,primes)
    return out % M


def qi_hat(j, primes):
    out = 1
    for i in range(len(primes)):
        if i != j:
            out = out*primes[i]
    return out


def inv_qi_hat(j, primes):
    p = primes[j]
    out = pow(qi_hat(j, primes), p-2, p)
    return out


def add(poly1_crt, poly2_crt, primes):
    out_crt = [0]*n*len(primes)
    for j in range(len(primes)):
        for i in range(n):
            idx = j*n+i
            out_crt[idx] = (poly1_crt[idx]+poly2_crt[idx]) % primes[j]
    return out_crt


def mul(poly1_crt, poly2_crt, primes):
    M = 1
    for p in primes:
        M = M * p
    poly1 = Polynomial(crt_to_coef(poly1_crt, primes))
    poly2 = Polynomial(crt_to_coef(poly2_crt,primes))
    mul_mod_cyc = utils.modulo_polynomial(poly1*poly2, cyc).coef
    mul_mod_cyc_mod_q = [int(x) % M for x in mul_mod_cyc]
    return coef_to_crt(mul_mod_cyc_mod_q, primes)


def mul_scalar(scalar, poly_crt, primes):
    out = [0]*n*len(primes)
    for i in range(n):
        for j in range(len(primes)):
            idx = i+j*n
            out[idx] = (poly_crt[idx]*scalar) % primes[j]
    return out


def gen_rand_poly_crt(primes):
    M = 1
    for p in primes:
        M = M*p
    poly_rand_coef = [random.randint(0, M) for _ in range(n)]
    return coef_to_crt(poly_rand_coef, primes)


def main():
    sk = gen_sk(qi)
    print("sk",sk)
    sk_crt = coef_to_crt(sk,qi)
    print("sk_crt: ", sk_crt)
    sk_sqr_crt = mul(sk_crt, sk_crt, qi)
    print("sk_sqr: ", sk_sqr_crt)
    print(crt_to_coef(sk_crt,qi))

    tmp1 = coef_to_crt([2,0,0,0,0,0,0,0], qi)
    print("tmp1", tmp1)
    tmp2 = coef_to_crt([3,0,0,0,0,0,0,0], qi)
    print("tmp2", tmp2)
    print("tmp1*tmp2=", mul(tmp1, tmp2, qi))

    print("mul_scalar: ", mul_scalar(2, sk_crt, qi))
    sk_sqr_time_p_crt = mul_scalar(P, sk_sqr_crt,qi)
    ax_crt = gen_rand_poly_crt(qipi)
    ax_times_sk_crt = mul(sk_crt, ax_crt, qipi)




if __name__ == '__main__':
    main()
