from crt import *
from utils import *

def gen_sk(primes, debug=False):
    M = prod(primes)
    if debug:
        print('M: ',M)
        # sk_coef = [0] * n
        # sk_coef[0] =  M - 1
        # sk = coef_to_crt(sk_coef, primes)
        sk = coef_to_crt([random.choice([0, 1]) for _ in range(n)], primes)
        print('sk_debug: ', sk)
        return sk
    # sk_coef = [random.choice([M-1, 0, 1]) for _ in range(n)]
    sk_coef = [random.choice([0, 1]) for _ in range(n)]
    return coef_to_crt(sk_coef, primes)



def modulo(poly: Polynomial, M: int):
    mul_mod_cyc = utils.modulo_polynomial(poly, cyc)
    out = [int(x) % M for x in mul_mod_cyc]
    return out

def he_mul(ax1, bx1, ax2, bx2, relin_key_ax_qipi, relin_key_bx_qipi, qi,pi, sk_qi_dbg,m1_qi_dbg, m2_qi_dbg):
    P_dbg = prod(pi)
    Q_dbg = prod(qi)
    QP_dbg = Q_dbg*P_dbg
    qipi=qi.copy()
    qipi.extend(pi)
    print('qipi: ',qipi)
    d0_qi = mul(bx1, bx2, qi)
    d0_dbg = Polynomial(crt_to_coef(d0_qi,qi))
    ax1_dbg = Polynomial(crt_to_coef(ax1,qi))
    ax2_dbg = Polynomial(crt_to_coef(ax2,qi))
    bx1_dbg = Polynomial(crt_to_coef(bx1, qi))
    bx2_dbg = Polynomial(crt_to_coef(bx2, qi))
    sk_dbg = Polynomial(crt_to_coef(sk_qi_dbg,qi))
    m1_dbg = Polynomial(crt_to_coef(m1_qi_dbg,qi))
    m2_dbg = Polynomial(crt_to_coef(m2_qi_dbg,qi))
    print('m1: ', m1_dbg.coef)
    print('m2: ', m2_dbg.coef)
    print('m1*m2: ', modulo(m1_dbg*m2_dbg,Q_dbg))
    print('-----')
    print_coef('d0_dbg: ', d0_dbg)
    tmp_dbg = ax1_dbg*ax2_dbg*sk_dbg*sk_dbg + (ax1_dbg*m2_dbg+ax2_dbg*m1_dbg)*sk_dbg + m1_dbg*m2_dbg
    print_coef('a1*a2*s^2+(a1m2+a2m1)*s+m1*m2: ', modulo(tmp_dbg, Q_dbg))
    print('b1*b2: ', modulo(bx1_dbg * bx2_dbg, Q_dbg))
    tmp = bx1_dbg*bx2_dbg-ax1_dbg*ax2_dbg*sk_dbg*sk_dbg -sk_dbg*(ax1_dbg*m2_dbg+ax2_dbg*m1_dbg)
    print('m1*m2= b1*b2-a1*a2*sk^2-sk*(a1*m2+a2*m1): ', modulo(tmp, Q_dbg))
    print('-----')

    d1_qi = add(mul(bx1, ax2, qi),mul(bx2, ax1, qi), qi)
    print('d1: ', crt_to_coef(d1_qi,qi))
    tmp_dbg = ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg
    print('a1b2+a2b1: ', modulo(tmp_dbg, Q_dbg))
    print('-----')

    d2_qi = mul(ax1, ax2, qi)
    print('d2: ', crt_to_coef(d2_qi,qi))
    print('a1*a2: ', modulo(ax1_dbg*ax2_dbg, Q_dbg))
    print('a1 mod Q: ', modulo(ax1_dbg, Q_dbg))
    print('a2 mod Q: ', modulo(ax2_dbg, Q_dbg))
    print('-----')

    d2_qipi = mod_up(d2_qi, qi, pi)
    print('d2 mod QP: ', crt_to_coef(d2_qi,qi))
    print('-----')
    relin_ax_dbg = Polynomial(crt_to_coef(relin_key_ax_qipi, qipi))
    relin_bx_dbg = Polynomial(crt_to_coef(relin_key_bx_qipi, qipi))
    print('a4 mod QP: ', crt_to_coef(relin_key_ax_qipi, qipi))
    print('relin_b-reli_a*sk: ', modulo(relin_bx_dbg - relin_ax_dbg*sk_dbg,QP_dbg))
    print('P*sk^2: ', modulo(P_dbg*sk_dbg*sk_dbg, QP_dbg))
    print('-----')

    d2_time_evk_ax_qipi = mul(d2_qipi, relin_key_ax_qipi, qipi)
    print('d2_time_evk_ax_qipi: ', crt_to_coef(d2_time_evk_ax_qipi, qipi))
    d2_dbg = Polynomial(crt_to_coef(d2_qipi,qipi))
    print('d2*relin_a: ', modulo(d2_dbg*relin_ax_dbg,QP_dbg))
    print('-----')

    d2_time_evk_bx_qipi = mul(d2_qipi, relin_key_bx_qipi, qipi)
    print('d2_time_evk_bx_qipi: ', crt_to_coef(d2_time_evk_bx_qipi, qipi))
    print('d2*relin_b: ', modulo(d2_dbg*relin_bx_dbg, QP_dbg))

    print('-----')
    d2_time_evk_ax_qi = mod_down(d2_time_evk_ax_qipi, qipi, qi)
    print('d2*a4 mod QP: ', crt_to_coef(d2_time_evk_ax_qipi,qipi))
    print('d2*a4 mod QP / P: ', [int(x/P_dbg)for x in crt_to_coef(d2_time_evk_ax_qipi,qipi)])
    print('d2*a4 mod Q: ', crt_to_coef(d2_time_evk_ax_qi,qi))
    print('d2*a4 / P: ', [int(x/P_dbg) for x in modulo(d2_dbg*relin_ax_dbg, QP_dbg)])
    print('a1*a2*a4/P: ', modulo(ax1_dbg * ax2_dbg * relin_ax_dbg / P_dbg, Q_dbg))
    print('a1*a2*a4/P: ', [int(x / P_dbg) for x in modulo(ax1_dbg * ax2_dbg * relin_ax_dbg, QP_dbg)])
    print('-----')



    d2_time_evk_bx_qi = mod_down(d2_time_evk_bx_qipi, qipi, qi)
    print('d2*b4 mod QP: ', crt_to_coef(d2_time_evk_bx_qipi, qipi))
    print('d2*b4 mod QP / P: ', [int(x / P_dbg) for x in crt_to_coef(d2_time_evk_bx_qipi, qipi)])
    print('d2*b4 mod Q: ', crt_to_coef(d2_time_evk_bx_qi, qi))
    print('d2*b4 / P: ', [int(x/P_dbg) for x in modulo(d2_dbg*relin_bx_dbg, QP_dbg)])
    print('a1*a2*b4/P mod Q: ', modulo(ax1_dbg * ax2_dbg * relin_bx_dbg/P_dbg, Q_dbg))

    tmp = bx1_dbg * bx2_dbg + ax1_dbg*ax2_dbg * relin_bx_dbg/P_dbg - sk_dbg*(ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg + ax1_dbg*ax2_dbg*relin_ax_dbg/P_dbg)
    print('m1*m2 = b1*b2 + a1*a2*b4/P - sk*(a1*b2 + a2*b1 + a1*a2*a4/P) mod Q:', modulo(tmp,Q_dbg))
    print('-----')
    B = bx1_dbg * bx2_dbg + ax1_dbg*ax2_dbg * relin_bx_dbg/P_dbg
    print('B= b1*b2 + a1*a2*b4/P', modulo(B, Q_dbg))
    A = ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg + ax1_dbg*ax2_dbg*relin_ax_dbg/P_dbg
    print('A: ', A)
    print('A= a1*b2 + a2*b1 + a1*a2*a4/P', modulo(A, Q_dbg))
    print('B-A*sk: ', modulo(B - A*sk_dbg, Q_dbg) )
    print('A*sk: ', modulo(A * sk_dbg, Q_dbg))
    print('A*sk mod qi: ', coef_to_crt(modulo(A * sk_dbg, Q_dbg), qi))
    print('-----')
    B2 = add(d0_qi, mod_down(mul(d2_qipi,relin_key_bx_qipi,qipi),qipi,qi),qi)
    print('B2: ', crt_to_coef(B2, qi))
    A2 = add(add(mul(ax1,bx2,qi),mul(ax2,bx1,qi),qi),mod_down(mul(d2_qipi,relin_key_ax_qipi,qipi),qipi,qi),qi)
    print('A2=',crt_to_coef(A2,qi))
    print('B2_qi-A2_qi*sk_qi: ', crt_to_coef(sub(B2,mul(A2,sk_qi_dbg,qi), qi),qi))
    print('A2_qi*sk_qi: ', mul(A2, sk_qi_dbg, qi))
    print('A2*sk coef : ', crt_to_coef(mul(A2, sk_qi_dbg, qi),qi))

    print('-----')
    print('A.coef: ', A.coef)
    print('B.coef: ', B.coef)
    print('A2: ', A2)
    print('B2: ', B2)
    # ax_out = coef_to_crt(A.coef, qi)
    # bx_out = coef_to_crt(modulo(B, Q_dbg), qi)

    ax_out = A2
    bx_out = B2
    # ax_out = add(d2_time_evk_ax_qi, d1_qi, qi)
    # bx_out = add(d2_time_evk_bx_qi, d0_qi, qi)
    bx_out_dbg = Polynomial(crt_to_coef(bx_out, qi))
    ax_out_dbg = Polynomial(crt_to_coef(ax_out, qi))
    print('out_b-out_a*s mod Q: ', modulo(bx_out_dbg-ax_out_dbg*sk_dbg,Q_dbg))
    return ax_out, bx_out


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


def gen_relin_key(sk_qipi, qipi, pi):
    P = prod(pi)
    ax_qipi = gen_rand_poly_crt(qipi)
    ax_time_sk_qipi = mul(ax_qipi, sk_qipi, qipi)
    sk_sqr_qipi = mul(sk_qipi, sk_qipi, qipi)
    sk_sqr_times_p_qipi = mul_scalar(P, sk_sqr_qipi, qipi)
    bx_qipi = add(ax_time_sk_qipi, sk_sqr_times_p_qipi, qipi)
    return ax_qipi, bx_qipi
