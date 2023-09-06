import numpy as np
from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial

def main():
    # print(RingElement(Polynomial([0,0,0,0,1]),8,10))

    dim = 16*2
    q = 100
    mat = [[1,1,1],[2,2,2],[3,3,3]]
    vec = [ 1, 1, 1]
    result = np.dot(np.array(mat), np.array(vec))
    print("result: ", result)
    # print(Polynomial([1,1,1]))
    m1 = Polynomial([0,1,2,3 ,4,5,6,  7,8,9])
    v1_arr = [0]*(dim//2)
    v1_arr[-1] = 1
    v1_arr[-2] = 1
    v1_arr[-3] = 1
    v1 = (-1)*Polynomial(v1_arr)
    m = RingElement(m1, dim, q)
    v = RingElement(v1, dim, q)
    print("m*v: ", m*v)

if __name__ == '__main__':
    main()