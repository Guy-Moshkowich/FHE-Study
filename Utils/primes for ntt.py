import sympy


def find_primes_near_p_mod_2n(p, n, num_primes=4):
    primes = []
    start = max(p - 2 * n, 2)  # Start checking from a point around p, not less than 2
    end = p + 2 * n  # End checking a bit further than p

    candidate = start
    while len(primes) < num_primes:
        if sympy.isprime(candidate) and candidate % (2 * n) == 1:
            primes.append(candidate)
        candidate += 1

    return primes


# Example usage
p = 1000
n = 8  # 2n = 16, so looking for primes that satisfy p ≡ 1 (mod 16)
primes = find_primes_near_p_mod_2n(p, n)
print(f"Primes near {p} that satisfy p ≡ 1 (mod {2 * n}): {primes}")
