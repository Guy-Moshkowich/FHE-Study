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

inv_P_qi = [0]*len(qi)
for i in range(len(qi)):
    P_qi = P % qi[i]
    inv_P = pow(P_qi, qi[i] - 2,  qi[i])
    inv_P_qi[i] = inv_P


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

def coef_to_crt_elm(elm_coef, primes):
    out = [0]*len(primes)
    for j in range(len(primes)):
        out[j] = elm_coef % primes[j]
    return out

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


def sub(poly1_crt, poly2_crt, primes):
    out_crt = [0]*n*len(primes)
    for j in range(len(primes)):
        for i in range(n):
            idx = j*n+i
            out_crt[idx] = (poly1_crt[idx]-poly2_crt[idx]) % primes[j]
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


def mod_up(poly_qi, qi, pi):
    poly_coef = crt_to_coef(poly_qi, qi)
    poly_pi = [0]*len(pi)*n
    for i in range(len(pi)):
        for j in range(n):
            idx = i*n+j
            poly_pi[idx] = poly_coef[j] % pi[i]
    poly_qi.extend(poly_pi)
    return poly_qi


def mod_down_elm(elm_qipi, qipi, qi):
    out = []
    pi = qipi[-(len(qipi)-len(qi)):]
    elm_pi = elm_qipi[-len(pi):]
    elm_P = crt_to_coef_elm(elm_pi, pi)
    a_qi = coef_to_crt_elm(elm_P, qi)
    for i in range(len(qi)):
        b_qi = elm_qipi[i]
        out.append(((b_qi - a_qi[i]) * inv_P_qi[i]) %qi[i])
    return out

def mod_down(poly_qipi, qipi, qi):
    out = [0]*n*len(qi)
    for i in range(n):
        elm_qipi = [0]*len(qipi)
        for j in range(len(qipi)):
            idx = j*n+i
            elm_qipi[j] = poly_qipi[idx]
        elm_qi = mod_down_elm(elm_qipi, qipi, qi)
        for j in range(len(qi)):
            idx = j * n + i
            out[idx] = elm_qi[j]
    return out