from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial
from bgv import BGV

def main():
    bgv = BGV(m_power=4, q=10000, p=2, N=10, max_added_noise=10)
    plaintext = RingElement(Polynomial([1, 1]), bgv.m, 2)
    one = RingElement(Polynomial([1, 0, 0]), bgv.m, 2)
    zero = RingElement(Polynomial([0]), bgv.m, 2)
    ctx = bgv.encrypt(plaintext)
    ctx_one = bgv.encrypt(one)
    print(ctx.get_noise())

    result = 1
    while (bgv.decrypt(ctx, bgv.secret_key) == plaintext) and (plaintext != one) and (plaintext != zero):
        ctx = ctx + ctx_one
        plaintext = plaintext + one
        print(ctx.get_noise())
    #     # print(ctx)
    #     # print( bgv.decrypt(ctx, bgv.secret_key) )
    #     # print(plaintext)
    #     # print('-----')

if __name__ == '__main__':
    main()