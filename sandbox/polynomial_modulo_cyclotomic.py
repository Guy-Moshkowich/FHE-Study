from numpy.polynomial import Polynomial
from Utils import utils

cyc = Polynomial([1,0,0,0,0,0,0,0,1])
Q = 97*193
d0 = Polynomial([18319, 11882, 1717, 18197, 18223, 6182, 15882, 9859, ])
b1 = Polynomial([18648, 13838, 575, 3144, 413, 8584, 13222, 4183, ])
b2 = Polynomial([18648, 7211, 1155, 3144, 413, 8584, 13222, 4183, ])

d2 = Polynomial([12789, 4219, 7330, 13141, 2779, 13508, 16245, 786, ])
a1 = Polynomial([2840, 13237, 8632, 16720, 14636, 17279, 18101, 14467, ])
a2 = Polynomial([2840, 13237, 8632, 16720, 14636, 17279, 18101, 14467, ])

d1 = Polynomial([8041, 1172, 2534, 14330, 6129, 1584, 6010, 3957, ])


tmp = utils.modulo_polynomial(a1*b2+a2*b1, cyc)
m = [x % Q for x in tmp.coef]
print(m)



