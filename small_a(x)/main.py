from RLWE.ring_element import RingElement


def main():
    q = 1000
    Q = 10000
    m = 100
    a = RingElement.random(m, q)
    print(a)


if __name__ == '__main__':
    main()