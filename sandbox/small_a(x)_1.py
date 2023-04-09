import math

import Utils.utils
from RLWE.ring_element import RingElement


def main():
    q = 128
    Q = 10000
    n = 2**5
    dnum = math.ceil(math.log2(q))
    a_swk = [RingElement.random(n, q) for i in range(dnum)]
    s = RingElement(Utils.utils.generate_ternary_polynomial(n//2, hamming_weight=64), n, q)
    phi_s = RingElement(Utils.utils.generate_ternary_polynomial(n//2, hamming_weight=64), n, q)
    I_swk = [RingElement.random_binary(n, q) for i in range(dnum)]
    e_swk = [RingElement.small_gauss(n, q) for i in range(dnum)]
    b_swk = [0]*dnum
    for j in range(dnum):
        b_swk[j] = a_swk[j]*s +  \
                   RingElement.const(2**j, n, q)*phi_s + \
                   RingElement.const(q, n, q)*I_swk[j] + \
                   e_swk[j]
    # assert correct decryption
    for j in range(dnum):
        msg = Utils.utils.modulo_int((b_swk[j] - a_swk[j]*s).poly, q)
        expected_msg = (RingElement.const(2**j, n, q)*phi_s).poly
        assert RingElement(msg - expected_msg, n, q).canonical_norm() < 3

    a = RingElement.random(n, q)
    a_decomp = Utils.utils.bit_decomp_poly(a.poly, dnum)

    sum = 0
    # RingElementa_decomp[j]*


    print('a=', a)
    print('a_decomp=', a_decomp)
    for j in range(dnum):
        print('a_swk[%d]=%s' % (j, a_swk[j]))
    for j in range(dnum):
        print('I_swk[%d]=%s' % (j, I_swk[j]))
    for j in range(dnum):
        print('b_swk[%d]=%s' % (j, b_swk[j]))

    print('s=', s)
    print('phi_s=', phi_s)




    # a_Q = RingElement.random(m, Q)
    # s_Q = RingElement.random_ternary(m, Q)
    # s_phi_Q = RingElement.random_ternary(m, Q)  #TODO: replace with rotated secret key
    # e = RingElement.small_gauss(m, Q)
    # b_Q = a_Q * s_Q + s_phi_Q + e
    #
    # a_q = a_Q.change_modulo(q)
    # b_q = b_Q.change_modulo(q)
    # print(b_Q-b_q)

    # s_q = s_Q.change_modulo(q)
    # print('a_q=', a_q)
    # print('a_Q=', a_Q)
    # print('b_q=', b_q)
    # print('b_Q=', b_Q)
    # print('|e|=', e.canonical_norm())
    # print('s_phi_Q=', s_phi_Q)
    # print('b_q-s_q*a_q', b_q-s_q*a_q)


if __name__ == '__main__':
    main()