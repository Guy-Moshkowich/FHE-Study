from RLWE.ring_element import RingElement


def noise_times_binary():
    log_n = 10
    n = 2 ** log_n
    q = 0x200000440001
    for i in range(10):
        e = RingElement.const(0, n, q)
        while e == RingElement.const(0, n, q):
            e = RingElement.small_gauss(n, q)
        a = RingElement.random_binary(n, q)
        print('a=', a)
        print('e=', e)
        print('|a|=', a.canonical_norm())
        print('|e|=', e.canonical_norm())
        print('|ae|=', (a * e).canonical_norm())
        print('----')


if __name__ == '__main__':
    noise_times_binary()
