from poly_rns import *
from utils import *

MAX_NOISE = 2

class Scheme:
    def __init__(self, n, qi, pi=[]):
        self.n = n
        self.qi = qi
        self.pi = pi
        self.Q = prod(qi)
        self.P = prod(pi)
        self.qipi = qi.copy()
        self.qipi.extend(pi)
        self.polyEval = Evaluator(n, self.qipi)
        self.sk_qi = self.gen_sk(Debug.POSITIVE_SK)
        self.relin_key_ax_qipi, self.relin_key_bx_qipi = self.gen_relin_key()

        self.inv_P_qi = [0] * len(qi)
        for i in range(len(qi)):
            self.inv_P_qi[i] = inv_mod(self.P, qi[i])

    def gen_sk(self, debug=0): # TODO: 11/7/24 Guy: sk with -1 values does not work. need to fix.
        if debug == Debug.POSITIVE_SK:
            # sk_coef = [0] * n
            # sk_coef[0] =  M - 1
            # sk = coef_to_crt(sk_coef, primes)
            sk = coef_to_crt([random.choice([0, 1]) for _ in range(n)], self.qi)
            print('sk_debug: ', sk)
            return sk
        sk_coef = [random.choice([0, -1, 1]) for _ in range(n)]
        return coef_to_crt(sk_coef, self.qi)



    def modulo(self, poly: Polynomial, M: int):
        mul_mod_cyc = utils.modulo_polynomial(poly, cyc)
        out = [int(x) % M for x in mul_mod_cyc]
        return out

    def he_mul(self, ax1, bx1, ax2, bx2, relin_key_ax_qipi, relin_key_bx_qipi, sk_qi_dbg,m1_qi_dbg, m2_qi_dbg):
        P_dbg = prod(self.pi)
        Q_dbg = prod(self.qi)
        QP_dbg = Q_dbg*P_dbg
        qipi=self.qi.copy()
        qipi.extend(self.pi)
        print('qipi: ',qipi)
        d0_qi = self.polyEval.mul(bx1, bx2, self.qi)
        d0_dbg = Polynomial(crt_to_coef(d0_qi,self.qi))
        ax1_dbg = Polynomial(crt_to_coef(ax1,self.qi))
        ax2_dbg = Polynomial(crt_to_coef(ax2,self.qi))
        bx1_dbg = Polynomial(crt_to_coef(bx1, self.qi))
        bx2_dbg = Polynomial(crt_to_coef(bx2, self.qi))
        sk_dbg = Polynomial(crt_to_coef(sk_qi_dbg,self.qi))
        m1_dbg = Polynomial(crt_to_coef(m1_qi_dbg,self.qi))
        m2_dbg = Polynomial(crt_to_coef(m2_qi_dbg,self.qi))
        print('m1: ', m1_dbg.coef)
        print('m2: ', m2_dbg.coef)
        print('m1*m2: ', self.modulo(m1_dbg*m2_dbg,Q_dbg))
        print('-----')
        print_coef('d0_dbg: ', d0_dbg)
        tmp_dbg = ax1_dbg*ax2_dbg*sk_dbg*sk_dbg + (ax1_dbg*m2_dbg+ax2_dbg*m1_dbg)*sk_dbg + m1_dbg*m2_dbg
        print_coef('a1*a2*s^2+(a1m2+a2m1)*s+m1*m2: ', self.modulo(tmp_dbg, Q_dbg))
        print('b1*b2: ', self.modulo(bx1_dbg * bx2_dbg, Q_dbg))
        tmp = bx1_dbg*bx2_dbg-ax1_dbg*ax2_dbg*sk_dbg*sk_dbg -sk_dbg*(ax1_dbg*m2_dbg+ax2_dbg*m1_dbg)
        print('m1*m2= b1*b2-a1*a2*sk^2-sk*(a1*m2+a2*m1): ', self.modulo(tmp, Q_dbg))
        print('-----')

        d1_qi = add(self.polyEval.mul(bx1, ax2, self.qi), self.polyEval.mul(bx2, ax1, self.qi), self.qi)
        print('d1: ', crt_to_coef(d1_qi,self.qi))
        tmp_dbg = ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg
        print('a1b2+a2b1: ', self.modulo(tmp_dbg, Q_dbg))
        print('-----')

        d2_qi = self.polyEval.mul(ax1, ax2, self.qi)
        print('d2_qi: ', d2_qi)
        print('d2_qi_mul: ', self.polyEval.mul(ax1, ax2, self.qi))
        print('d2_coef: ', crt_to_coef(d2_qi,self.qi))
        print('a1*a2: ', self.modulo(ax1_dbg*ax2_dbg, Q_dbg))
        print('a1 mod Q: ', self.modulo(ax1_dbg, Q_dbg))
        print('a2 mod Q: ', self.modulo(ax2_dbg, Q_dbg))
        print('-----')

        d2_qipi = mod_up(d2_qi, self.qi, self.pi, self.n)
        print('d2 mod QP: ', crt_to_coef(d2_qi,self.qi))
        print('-----')
        relin_ax_dbg = Polynomial(crt_to_coef(relin_key_ax_qipi, qipi))
        relin_bx_dbg = Polynomial(crt_to_coef(relin_key_bx_qipi, qipi))
        print('a4 mod QP: ', crt_to_coef(relin_key_ax_qipi, qipi))
        print('relin_b-reli_a*sk: ', self.modulo(relin_bx_dbg - relin_ax_dbg*sk_dbg,QP_dbg))
        print('P*sk^2: ', self.modulo(P_dbg*sk_dbg*sk_dbg, QP_dbg))
        print('-----')

        d2_time_evk_ax_qipi = self.polyEval.mul(d2_qipi, relin_key_ax_qipi, qipi)
        print('d2_time_evk_ax_qipi: ', crt_to_coef(d2_time_evk_ax_qipi, qipi))
        d2_dbg = Polynomial(crt_to_coef(d2_qipi,qipi))
        print('d2*relin_a: ', self.modulo(d2_dbg*relin_ax_dbg,QP_dbg))
        print('-----')
        print('d2*relin_b: ', self.modulo(d2_dbg*relin_bx_dbg, QP_dbg))

        print('-----')
        d2_time_evk_ax_qi = mod_down(d2_time_evk_ax_qipi, self.qi, self.pi, self.inv_P_qi, self.n)
        print('d2*a4 mod QP: ', crt_to_coef(d2_time_evk_ax_qipi,qipi))
        print('d2*a4 mod QP / P: ', [int(x/P_dbg)for x in crt_to_coef(d2_time_evk_ax_qipi,qipi)])
        print('d2*a4 mod Q: ', crt_to_coef(d2_time_evk_ax_qi,self.qi))
        print('d2*a4 / P: ', [int(x/P_dbg) for x in self.modulo(d2_dbg*relin_ax_dbg, QP_dbg)])
        print('a1*a2*a4/P: ', self.modulo(ax1_dbg * ax2_dbg * relin_ax_dbg / P_dbg, Q_dbg))
        print('a1*a2*a4/P: ', [int(x / P_dbg) for x in self.modulo(ax1_dbg * ax2_dbg * relin_ax_dbg, QP_dbg)])
        print('-----')

        tmp = bx1_dbg * bx2_dbg + ax1_dbg*ax2_dbg * relin_bx_dbg/P_dbg - sk_dbg*(ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg + ax1_dbg*ax2_dbg*relin_ax_dbg/P_dbg)
        print('m1*m2 = b1*b2 + a1*a2*b4/P - sk*(a1*b2 + a2*b1 + a1*a2*a4/P) mod Q:', self.modulo(tmp,Q_dbg))
        print('-----')
        B = bx1_dbg * bx2_dbg + ax1_dbg*ax2_dbg * relin_bx_dbg/P_dbg
        print('B= b1*b2 + a1*a2*b4/P', self.modulo(B, Q_dbg))
        A = ax1_dbg*bx2_dbg + ax2_dbg*bx1_dbg + ax1_dbg*ax2_dbg*relin_ax_dbg/P_dbg
        print('A: ', A)
        print('A= a1*b2 + a2*b1 + a1*a2*a4/P', self.modulo(A, Q_dbg))
        print('B-A*sk: ', self.modulo(B - A*sk_dbg, Q_dbg) )
        print('A*sk: ', self.modulo(A * sk_dbg, Q_dbg))
        print('A*sk mod qi: ', coef_to_crt(self.modulo(A * sk_dbg, Q_dbg), self.qi))
        print('-----')
        print('mul(d2_qipi,relin_key_bx_qipi,qipi): ', self.polyEval.mul(d2_qipi, relin_key_bx_qipi, qipi))

        # B2 = d0 + mod_down(mod_up(d2) * relin_key_bx mod QP) mod Q
        mul_tmp = self.polyEval.mul(d2_qipi, relin_key_bx_qipi, qipi)
        mod_down_tmp = mod_down(mul_tmp, self.qi, self.pi, self.inv_P_qi, self.n)
        B2 = add(d0_qi,mod_down_tmp , self.qi)

        print('B2: ', crt_to_coef(B2, self.qi))
        print('B2_qi: ', B2)
        # print('B2_qi_tmp: ', B2_tmp)

        mul_tmp = self.polyEval.mul(d2_qipi, relin_key_ax_qipi, qipi)
        mod_down_tmp = mod_down(mul_tmp, self.qi, self.pi, self.inv_P_qi, self.n)
        A2 = add(d1_qi, mod_down_tmp, self.qi)
        print('A2=',crt_to_coef(A2,self.qi))
        print('B2_qi-A2_qi*sk_qi: ', crt_to_coef(sub(B2, self.polyEval.mul(A2, sk_qi_dbg, self.qi), self.qi), self.qi))
        print('A2_qi*sk_qi: ', self.polyEval.mul(A2, sk_qi_dbg, self.qi))
        print('A2*sk coef : ', crt_to_coef(self.polyEval.mul(A2, sk_qi_dbg, self.qi), self.qi))

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
        bx_out_dbg = Polynomial(crt_to_coef(bx_out, self.qi))
        ax_out_dbg = Polynomial(crt_to_coef(ax_out, self.qi))
        print('out_b-out_a*s mod Q: ', self.modulo(bx_out_dbg-ax_out_dbg*sk_dbg,Q_dbg))
        return ax_out, bx_out

    def he_add(self, ax1, bx1, ax2, bx2):
        return add(ax1, ax2, self.qi), add(bx1, bx2, self.qi)

    def encrypt(self, pt_qi, debug=0):
        if debug == Debug.DISABLE_NOISE:
            noise = [0]*n*len(self.qi)
            print('noise: ', noise)
        else:
            noise = gen_noise_crt(MAX_NOISE, self.qi)
        ax_qi = gen_rand_poly_crt(self.qi)
        ax_time_sk_qi = self.polyEval.mul(ax_qi, self.sk_qi, self.qi)
        bx_qi = add(ax_time_sk_qi, pt_qi, self.qi)
        bx_qi = add(bx_qi, noise, self.qi)
        return ax_qi, bx_qi


    def decrypt(self, ax_qi, bx_qi):
        ax_time_sk_qi = self.polyEval.mul(ax_qi, self.sk_qi, self.qi)
        return sub(bx_qi,ax_time_sk_qi, self.qi)

    def gen_relin_key(self):
        P = prod(self.pi)
        ax_qipi = gen_rand_poly_crt(self.qipi)
        sk_qipi = mod_up(self.sk_qi, self.qi, self.pi, n)
        ax_time_sk_qipi = self.polyEval.mul(ax_qipi, sk_qipi, self.qipi)
        sk_sqr_qipi = self.polyEval.mul(sk_qipi, sk_qipi, self.qipi)
        sk_sqr_times_p_qipi = mul_scalar(P, sk_sqr_qipi, self.qipi)
        bx_qipi = add(ax_time_sk_qipi, sk_sqr_times_p_qipi, self.qipi)
        return ax_qipi, bx_qipi


