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


# def fast_base_conv(a, from_qi, to_pi):
#     out = []
#     q = prod(from_qi)
#     for pi in to_pi:
#         s = 0
#         v = 0
#         for j in range(len(from_qi)):
#             tmp1 = ((a[j] * inv_hat(j, from_qi)) % from_qi[j]) % pi
#             tmp2 = hat(j, from_qi) % pi
#             s = (s + (tmp1*tmp2) % pi) % pi
#             v += ((a[j] * inv_hat(j, from_qi)) % from_qi[j])/from_qi[j]
#         v = round(v)
#         out.append(((s % pi) - (v % pi)*(q % pi)) % pi)
#     return out


def fast_base_conv_poly_kernel(poly_out_, poly_in_, idx, pi_idx,
                               from_qi_, to_pi_, q_mod_pi, inv_hat_qi_mod_qi_,
                                hat_qi_mod_pj_flat_, n):
    s = 0
    v = 0
    pi = to_pi_[pi_idx]
    for j in range(len(from_qi_)):
        a = poly_in_[idx + j * n]
        y_j = (a * inv_hat_qi_mod_qi_[j]) % from_qi_[j]
        tmp1 = y_j % pi
        tmp2 = hat_qi_mod_pj_flat_[pi_idx + j*len(from_qi_)]
        tmp3 = (tmp1 * tmp2) % pi
        s = (s + tmp3) % pi
        v += y_j / from_qi_[j]
    v = round(v) % pi
    tmp4 = (v * q_mod_pi[pi_idx]) % pi
    tmp5 = (s - tmp4) % pi
    #print(idx, ",", pi_idx, " :::: ", idx + pi_idx * n,",",tmp5)
    poly_out_[idx + pi_idx * n] = tmp5


def fast_base_conv_poly(poly_in, from_qi, to_pi, n):
    poly_out = [0] * n * len(to_pi)
    q = prod(from_qi)
    q_mod_pi = [(q % pi) for pi in to_pi]
    inv_hat_qi_mod_qi = [inv_hat(j, from_qi) for j in range(len(from_qi))]
    hat_qi = [hat(j, from_qi) for j in range(len(from_qi))]
    hat_qi_mod_pj_flat_ = []
    for i in range(len(from_qi)):
        for j in range(len(to_pi)):
            hat_qi_mod_pj_flat_.append(hat_qi[i] % to_pi[j])
    for coef_idx in range(n):
        for pi_idx in range(len(to_pi)):
            fast_base_conv_poly_kernel(poly_out, poly_in, coef_idx, pi_idx,
                                       from_qi, to_pi, q_mod_pi,inv_hat_qi_mod_qi,
                                       hat_qi_mod_pj_flat_,n)
    return poly_out


def coef_to_crt(poly_coef, primes):
    crt = [0]*(n*len(primes))
    for j in range(len(primes)):
        for i in range(len(poly_coef)):
            idx = j*n+i
            crt[idx] = poly_coef[i] % primes[j]
    return crt


def crt_to_coef(poly_crt, primes):
    out = []
    n = len(poly_crt)//len(primes)
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
    out = inv_mod(hat(j, primes) % p , p)
    return out


def inv_mod(x, q):
    return pow(x, q - 2, q)


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
    poly_rand_coef = [random.randint(0, M) for _ in range(n)]
    return coef_to_crt(poly_rand_coef, primes)


def mod_up(poly_qi, qi, pi, n):
    poly_pi = fast_base_conv_poly(poly_qi, qi, pi, n)
    out = poly_qi.copy()
    out.extend(poly_pi)
    return out

def mod_down_kernel(out_, idx_coef, crt_idx, diff_, inv_P_qi_,qi_):
    idx = crt_idx * n + idx_coef
    out_[idx] = (diff_[idx] * inv_P_qi_[crt_idx]) % qi_[crt_idx]


def mod_down(poly_qipi, qi,pi, inv_P_qi, n):
    out_ = [0] * n * len(qi)
    poly_pi = [0]* n * len(pi)
    poly_qi = [0]* n * len(pi)
    for idx_coef in range(n):
        # TODO: 14/7/24 Guy: check if it is faster to extract polynomials by primes in GPU rather in CPU.
        for j in range(len(pi)):
            idx_qipi = j * n + idx_coef + len(qi)*n
            idx_pi = j * n + idx_coef
            poly_pi[idx_pi] = poly_qipi[idx_qipi]
        for j in range(len(qi)):
            idx_qipi = j * n + idx_coef
            idx_qi = j * n + idx_coef
            poly_qi[idx_qi] = poly_qipi[idx_qipi]
    poly_pi_qi = fast_base_conv_poly(poly_pi, pi, qi, n)
    diff_ = sub(poly_qi, poly_pi_qi, qi)
    print('diff: ', diff_)
    for idx_coef in range(n):
        for crt_idx in range(len(qi)):
            mod_down_kernel(out_, idx_coef, crt_idx, diff_, inv_P_qi, qi)
    return out_

