import math

import sympy

def find_negacyclic_primes(order, limit):
    n = order * 2  # We are looking for a 32nd root of unity because g^16 should be -1 (hence g^32 = 1)
    primes_with_negacyclic_roots = []

    # Check primes p such that p-1 is divisible by 32
    for p in sympy.primerange(1, limit):
        if (p - 1) % n != 0:
            continue  # Skip primes that do not have a group large enough

        # Check for any generator g where g^order == -1 mod p
        found = False
        for g in range(2, p):
            if pow(g, order, p) == p - 1 and pow(g, n, p) == 1:
                primes_with_negacyclic_roots.append((p, g))
                found = True
                break

    return primes_with_negacyclic_roots

# Example: Find all primes with negacyclic roots of order 16 below 1000
primes_info = find_negacyclic_primes(16, 1000)
for prime, generator in primes_info:
    print(f"Prime: {prime}, Generator: {generator}")


prime = 17
n = 16
k=2
while k < prime:
    root = k;
    for i in range(0, int(math.log2(n))):
        root = (root **2) % prime
    print(root)
    if (root + 1) % prime == 0:
        print('root=',k)
    k = k+1



prime = 193
n = 16
root = 8;
for i in range(0, int(math.log2(n))):
    root = (root **2) % prime
print(root)
if (root + 1) % prime == 0:
    print('root=', root)



# print(root)
# 6071469**4096+1 == 1?
# print(math.log((6071469-1)/2,5))