from RLWE.ring_element import RingElement
from numpy.polynomial import Polynomial
from Utils import utils

def load(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            # Convert each line to an integer and append to the list
            data.append(int(line.strip()))
    return data

# Q = 9007199255019521*9007199255347201*18014398510645249*18014398510661633
# n = 8192
# phi_m = utils.build_cyclotomic_poly(2*n)
# s = Polynomial(load("s_coeffs.txt"))
# print(s)
# ax_array = load("ax_coeffs.txt")
# # print(len(ax_array))
# ax = Polynomial(ax_array)
# # print(ax)
#
# bx = Polynomial(load("bx_coeffs.txt"))
# # print(bx)
#
# dec = bx-ax*s
# dec = utils.modulo_polynomial(dec, phi_m)
# dec = [x % Q for x in dec.coef]
# print(dec)

n = 16
Q = 97*193
phi_m = utils.build_cyclotomic_poly(2*n)
print('phi=', phi_m)
pt = Polynomial([965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965, 965])
ax = Polynomial([2840, 13237, 8632, 16720, 14636, 17279, 18101, 14467, 18161, 17081, 17500, 8432, 12012, 15280, 1806, 476])
bx = Polynomial([2841, 5487, 8635, 2006, 14641, 1449, 6, 4263, 8, 17091, 10, 8444, 12025, 15294, 16931, 15])
sk = Polynomial([1, 18720, 1, 18720, 1, 18720, 0, 18720, 0, 1, 0, 1, 1, 1, 18720, 0])
dec = utils.modulo_polynomial((pt+ax*sk), phi_m)
dec = [x % Q for x in dec.coef]
print(dec)
dec = bx-sk*ax
dec = utils.modulo_polynomial(dec, phi_m)
dec = [x % Q for x in dec.coef]
print(dec)


pk0 = [15881, 0, 8632, 0, 4085, 17279, 18101, 0, 0, 0, 17500, 0, 12012, 0, 0, 0]
pk1 = [2840, 13237, 8632, 16720, 14636, 17279, 18101, 14467, 18161, 17081, 17500, 8432, 12012, 15280, 1806, 476]
sk2 = [18720, 0, 1, 0, 18720, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0]
tmp = Polynomial(sk2)*Polynomial(pk1)
tmp = utils.modulo_polynomial(tmp, phi_m)
tmp = [x % Q for x in tmp.coef]
print('tmp=',tmp)


## x^{n-1} * x = x^n=-1
# x = [0]*16
# x[-1]=1
# y = [0]*16
# y[1]=1
# tmp = Polynomial(x)*Polynomial(y)
# print('tmp=', tmp)
# tmp = utils.modulo_polynomial(tmp, phi_m)
# print('tmp=', tmp)
# tmp = [x % Q for x in tmp.coef]
# print('tmp=', tmp)


