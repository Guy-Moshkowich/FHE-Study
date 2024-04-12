from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial
import sympy
from sympy.abc import x
import numpy as np


def gen_poly1(mat, k, w):
    A = 0
    n = w*w
    print('w=', w)
    for i in range(w):
        for j in range(w):
            if ((i - k) * w + j) < 0:
                A = A + (-1 * mat[i, j] * (x ** (n + ((i - k) * w + j))))
            else:
                A = A + (mat[i, j] * (x ** ((i - k) * w + j)))
    print(A)
    return A


def gen_poly2(mat, k, w):
    B = 0
    n = w*w
    print('k=', k)
    for i in range(k):
        for j in range(k):
            if (w * k - (i * w + j)) < 0:
                B = B + (-1 * mat[i, j] * x ** (n + (w * k - (i * w + j))))
            else:
                B = B + (mat[i, j] * x ** (w * k - (i * w + j)))
    print(B)
    return B

def main():

    a = ["a_{"+str(i)+"}" for i in range(10)]
    print(a[0])
    I_mat =  sympy.symbols("a_{1\,1} a_{1\,2} a_{1\,3} a_{2\,1} a_{2\,2} a_{2\,3} a_{3\,1} a_{3\,2} a_{3\,3}")
    I_mat = np.array(I_mat).reshape(3, 3)
    print(I_mat)
    I = gen_poly1(I_mat, k=2, w=3)
    K_mat = sympy.symbols("b_{1\,1} b_{1\,2} b_{2\,1} b_{2\,2}")
    K_mat = np.array(K_mat).reshape(2, 2)
    K = gen_poly2(K_mat, k=2, w=3)
    print(sympy.collect(sympy.expand(I*K), x))

    # A = sympy.Matrix([[a_11,a_12],[a_21,a_22]])
    # B = sympy.Matrix([[b_11,b_12],[b_21,b_22]])
    # # print(A*B)




'''
    w = 3
    k = 2
    n = 25
    X = var('x')
    I = np.array(var('w', n=w ^ 2, latex_name='w')).reshape((w, w))
    K = np.array(var('k', n=k ^ 2, latex_name='k')).reshape((k, k))
    print(I)
    print(K)
    print(X)
    A = 0
    for i in range(w):
        for j in range(w):
            if ((i - k) * w + j) < 0:
                A = A + (-1 * I[i, j] * X ^ (n + ((i - k) * w + j)))
            else:
                A = A + (I[i, j] * X ^ ((i - k) * w + j))
    print(A)
    B = 0
    for i in range(k):
        for j in range(k):
            if (w * k - (i * w + j)) < 0:
                B = B + (-1 * K[i, j] * X ^ (n + (w * k - (i * w + j))))
            else:
                B = B + (K[i, j] * X ^ (w * k - (i * w + j)))
    print(B)
    R = (A * B).full_simplify()
    print(R)
    '''

def tmp():
    ''' I=
    1   2  3  4
    4   5  6  7
    8   9 10  11
    12 13 14  15

    K=
    1 1
    1 1
    '''
    w=4
    k=2
    N=w*w
    I_elm = [0]*N
    K_elm = [0]*N
    for i in range(w):
        for j in range(w):
            if (i-k)*w+j >= 0:
                I_elm[(i-k)*w+j] = w*i+j
            else:
                I_elm[((i-k)*w+j)] = -(w*i+j)
    print(I_elm)
    for i in range(k):
        for j in range(k):
            K_elm[w*k-i*w-j] = 1
    print(K_elm)

    I = RingElement(Polynomial(I_elm), N, 1000)
    K = RingElement(Polynomial(K_elm), N, 1000)
    print(I*K)



if __name__ == '__main__':
    main()