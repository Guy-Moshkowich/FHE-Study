import numpy as np


def find_primitive_root(p):
    """Find a primitive root modulo p."""
    # if not is_prime(p):
    #     raise ValueError("p must be prime")

    # Calculate the factors of p-1
    factors = []
    phi = p - 1
    temp = phi
    i = 2
    while i * i <= temp:
        if temp % i == 0:
            factors.append(i)
            while temp % i == 0:
                temp //= i
        i += 1
    if temp > 1:
        factors.append(temp)

    # Find the smallest primitive root
    for g in range(2, p):
        ok = True
        for f in factors:
            if pow(g, f, p) == 1:
                ok = False
                break
        if ok:
            return g
    return -1

def ntt(a, p, primitive_root):
    """Compute the Number Theoretic Transform of a."""
    n = len(a)
    if n == 1:
        return a

    # Calculate the twiddle factors
    wn = pow(primitive_root, (p - 1) // n, p)
    w = 1

    # Recursive NTT
    a_even = ntt([a[i] for i in range(0, n, 2)], p, primitive_root)
    a_odd = ntt([a[i] for i in range(1, n, 2)], p, primitive_root)

    # Combine the results
    out = [0] * n
    for i in range(n // 2):
        out[i] = (a_even[i] + w * a_odd[i]) % p
        out[i + n // 2] = (a_even[i] - w * a_odd[i]) % p
        w = (w * wn) % p
    return out


def inverse_ntt(a, p, primitive_root):
    """Compute the Inverse Number Theoretic Transform of a."""
    n = len(a)
    if n == 1:
        return a

    # Calculate the twiddle factors
    wn = pow(primitive_root, (p - 1) // n, p)
    w = 1

    # Recursive Inverse NTT
    a_even = inverse_ntt([a[i] for i in range(0, n, 2)], p, primitive_root)
    a_odd = inverse_ntt([a[i] for i in range(1, n, 2)], p, primitive_root)

    # Combine the results
    out = [0] * n
    for i in range(n // 2):
        out[i] = (a_even[i] + w * a_odd[i]) % p
        out[i + n // 2] = (a_even[i] - w * a_odd[i]) % p
        w = (w * wn) % p

    # Divide by n to get the final result
    return [x * pow(n, p - 2, p) % p for x in out]

p = 167
primitive_root = find_primitive_root(p)
a = [1, 2, 3, 4]
ntt_result = ntt(a, p, primitive_root)
print(ntt_result)  # Output: [10, 40, 70, 100]
inverse_ntt_result = inverse_ntt(ntt_result, p, primitive_root)
print(inverse_ntt_result)  # Output: [1, 2, 3, 4]



x = [0]*16
x[1] = 2
y = [0]*16
y[2] = 3
prime = 97
primitve_root = 19
inv_primitve_root=46
# print(find_primitive_root(97))



ntt_x = ntt(x, prime, primitve_root)
ntt_y = ntt(y, prime, primitve_root)
result = np.multiply(ntt_x, ntt_y)
print(" result:", inverse_ntt(result, prime, primitve_root))