from numpy.polynomial import Polynomial
from Utils import utils
import random
from utils import *

random.seed(0)
n = 8
cyc = Polynomial([1,0,0,0,0,0,0,0,1])


class Evaluator:
    def __init__(self, n, primes):
        self.n = n
        self.primes = primes
        self.ntt_mat = {p:utils.ntt_matrix(n, p) for p in primes}
        self.inv_ntt_mat = {p: utils.inv_ntt_matrix(n, p) for p in primes}

    def mul(self, poly1_crt, poly2_crt, primes):
        out = []
        for i in range(len(primes)):
            x = utils.ntt_mat(primes[i], poly1_crt[i*n:i*n+n], self.ntt_mat[primes[i]])
            y = utils.ntt_mat(primes[i], poly2_crt[i*n:i*n+n], self.ntt_mat[primes[i]])
            ntt_res = [x[j] * y[j] % primes[i] for j in range(n)]
            res_i = [(x % primes[i]) for x in utils.inv_ntt_mat(primes[i], ntt_res, self.inv_ntt_mat[primes[i]])]
            out.extend(res_i)
        return out


def fast_base_conv(a, from_qi, to_pi):
    out = []
    for pi in to_pi:
        sum = 0
        for j in range(len(a)):
            tmp1 = (a[j] * inv_hat(j, from_qi)) % pi
            tmp2 = hat(j, from_qi) % pi
            sum += (tmp1*tmp2) % pi
        out.append(sum % pi)
    return out


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
    M = prod(primes)
    out = 0
    for j in range(len(crt_elm)):
        out = out + crt_elm[j] * hat(j, primes) * inv_hat(j, primes)
    return out % M

def coef_to_crt_elm(elm_coef, primes):
    out = [0]*len(primes)
    for j in range(len(primes)):
        out[j] = elm_coef % primes[j]
    return out

def hat(j, primes):
    out = 1
    for i in range(len(primes)):
        if i != j:
            out = out*primes[i]
    return out


def inv_hat(j, primes):
    p = primes[j]
    out = pow(hat(j, primes), p - 2, p)
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


def mul_scalar(scalar, poly_crt, primes):
    out = [0]*n*len(primes)
    for i in range(n):
        for j in range(len(primes)):
            idx = i+j*n
            out[idx] = (poly_crt[idx]*scalar) % primes[j]
    return out

def gen_noise_crt(max_noise, primes):
    poly_rand_coef = []
    for _ in range(n):
        if random.randint(0,10) > 8:
            poly_rand_coef.append(1)
            # if random.randint(0,1) ==1:
            #     poly_rand_coef.append(1)
            # else:
            #     poly_rand_coef.append(-1)
        else:
            poly_rand_coef.append(0)
    # poly_rand_coef = [random.randint(0, max_noise + 1) for _ in range(n)]
    return coef_to_crt(poly_rand_coef, primes)


def gen_rand_poly_crt(primes):
    M = prod(primes)
    # if debug:
    #     # rnd_poly_coef = [0]*n
    #     # rnd_poly_coef[0] = 30
    #     # rnd_poly = coef_to_crt(rnd_poly_coef,primes)
    #     rnd_poly_coef=[random.randint(0, M) for _ in range(n)]
    #     rnd_poly = coef_to_crt(rnd_poly_coef, primes)
    #     return rnd_poly

    poly_rand_coef = [random.randint(0, M) for _ in range(n)]
    return coef_to_crt(poly_rand_coef, primes)


def mod_up(poly_qi, qi, pi):
    poly_coef = crt_to_coef(poly_qi, qi)
    poly_pi = [0]*len(pi)*n
    for i in range(len(pi)):
        for j in range(n):
            idx = i*n+j
            poly_pi[idx] = poly_coef[j] % pi[i]
    out = poly_qi.copy()
    out.extend(poly_pi)
    return out


def mod_down_elm(elm_qipi, qipi, qi):
    pi = qipi[-(len(qipi)-len(qi)):]
    P = 1
    for p in pi:
        P = P*p

    inv_P_qi = [0] * len(qi)
    for i in range(len(qi)):
        P_qi = P % qi[i]
        inv_P = pow(P_qi, qi[i] - 2, qi[i])
        inv_P_qi[i] = inv_P

    out = []
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