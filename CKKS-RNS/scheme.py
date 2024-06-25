from crt import *

# def he_mul(ax1, bx1, ax2, bx2):


def he_add(ax1, bx1, ax2, bx2, qi):
    return add(ax1,ax2,qi), add(bx1,bx2,qi)


def encrypt(pt_qi, sk_qi, qi):
    ax_qi = gen_rand_poly_crt(qi)
    ax_time_sk_qi = mul(ax_qi, sk_qi, qi)
    bx_qi = add(ax_time_sk_qi, pt_qi, qi)
    return ax_qi, bx_qi


def decrypt(ax_qi, bx_qi, sk_qi, qi):
    ax_time_sk_qi = mul(ax_qi, sk_qi, qi)
    return sub(bx_qi,ax_time_sk_qi, qi)


def gen_relin_key(sk, qi, pi):
    ax_qipi = gen_rand_poly_crt(qipi)
    sk_qipi = mod_up(sk, qi, pi)
    ax_time_sk_qipi = mul(ax_qipi, sk_qipi, qipi)
    sk_sqr_piqi = mul(sk_qipi, sk_qipi, qipi)
    sk_sqr_times_p_piqi = mul_scalar(P, sk_sqr_piqi, qipi)
    bx_qipi = add(ax_time_sk_qipi, sk_sqr_times_p_piqi, qipi)
    return ax_qipi, bx_qipi