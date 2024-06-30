import numpy as np

prime = 97
primitve_root = 19
assert 19**16 % 97 == 97-1

inv_primitve_root=46
assert (inv_primitve_root*primitve_root) % prime == 1

def ntt(x):
    N = len(x)
    mu = inv_primitve_root
    if N <= 1:
        return x
    even = ntt(x[0::2])
    odd = ntt(x[1::2])
    T = [(pow(mu, k) * odd[k] )% prime for k in range(N // 2)]
    return [(even[k] + T[k])%prime for k in range(N // 2)] + [(even[k] - T[k])%prime for k in range(N // 2)]

def intt(X):
    N = len(X)
    mu = primitve_root
    if N <= 1:
        return X
    even = intt(X[0::2])
    odd =  intt(X[1::2])
    combined = [0] * N
    for k in range(N // 2):
        t = (pow(mu, k) * odd[k])%prime
        combined[k] = (even[k] + t)%prime
        combined[k + N // 2] = (even[k] - t)%prime
    return [x / 2 for x in combined]  # Normalize by N here for true inverse

# Example usage:
x = [0]*16
x[1] = 2
y = [0]*16
y[2] = 3

ntt_x = ntt(x)
ntt_y = ntt(y)
result = np.multiply(ntt_x, ntt_y)
print(result)
print(" result:",intt(result))
