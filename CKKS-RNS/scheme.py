from crt import *


def gen_sk(primes, debug=False):
    M = 1
    for p in primes:
        M = M * p
    if debug:
        print('M: ',M)
        sk_coef = [0] * len(primes)
        sk_coef[0] = 1
        sk = coef_to_crt(sk_coef, primes)
        print('sk_debug: ', sk)
        return sk
    sk_coef = [random.choice([M-1, 0, 1]) for _ in range(n)]
    return coef_to_crt(sk_coef, primes)

def he_mul(ax1, bx1, ax2, bx2, relin_key_ax_qipi, relin_key_bx_qipi, qi,qipi):
    d0_qi = mul(bx1, bx2, qi)
    d1_qi = add(mul(bx1, ax2, qi),mul(bx2, ax1, qi), qi)
    d2_qi = mul(ax1, ax2, qi)
    d2_qipi = mod_up(d2_qi, qi, qipi)
    d2_time_evk_ax_qipi = mul(d2_qipi, relin_key_ax_qipi, qipi)
    d2_time_evk_bx_qipi = mul(d2_qipi, relin_key_bx_qipi, qipi)
    d2_time_evk_ax_qi = mod_down(d2_time_evk_ax_qipi, qipi, qi)
    d2_time_evk_bx_qi = mod_down(d2_time_evk_bx_qipi, qipi, qi)
    ax_out = add(d2_time_evk_ax_qi, d1_qi, qi)
    bx_out = add(d2_time_evk_bx_qi, d0_qi, qi)
    return ax_out, bx_out


def he_add(ax1, bx1, ax2, bx2, qi):
    return add(ax1,ax2,qi), add(bx1,bx2,qi)


def encrypt(pt_qi, sk_qi, qi, debug=False):
    ax_qi = gen_rand_poly_crt(qi, debug)
    ax_time_sk_qi = mul(ax_qi, sk_qi, qi)
    bx_qi = add(ax_time_sk_qi, pt_qi, qi)
    return ax_qi, bx_qi


def decrypt(ax_qi, bx_qi, sk_qi, qi):
    ax_time_sk_qi = mul(ax_qi, sk_qi, qi)
    return sub(bx_qi,ax_time_sk_qi, qi)


def gen_relin_key(sk_qipi, qipi, pi, debug=False):
    P = 1
    for p in pi:
        P = P*p

    ax_qipi = gen_rand_poly_crt(qipi, debug)
    ax_time_sk_qipi = mul(ax_qipi, sk_qipi, qipi)
    sk_sqr_qipi = mul(sk_qipi, sk_qipi, qipi)
    sk_sqr_times_p_qipi = mul_scalar(P, sk_sqr_qipi, qipi)
    bx_qipi = add(ax_time_sk_qipi, sk_sqr_times_p_qipi, qipi)
    return ax_qipi, bx_qipi
