from Utils import utils

n = 8
p = 97
inv_ntt = utils.inv_ntt_matrix(n, p)
print(inv_ntt)
ntt = utils.ntt_matrix(n, p)
print(ntt)

n = 8
p = 193
inv_ntt = utils.inv_ntt_matrix(n, p)
print(inv_ntt)
ntt = utils.ntt_matrix(n, p)
print(ntt)
