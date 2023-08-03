from typing import List
import random
import numpy as np


class LweParams:
    modulo: int
    dim: int

    def __init__(self, modulo: int, dim: int):
        self.modulo = modulo
        self.dim = dim

    def gen_secret_key(self):
        return  [int(random.uniform(0, self.modulo)) for _ in range(self.dim)]

    def encrypt(self, m: int, sk: List[int]) -> type['LweCiphertext']:
        a = [int(random.uniform(0, self.modulo)) for _ in range(self.dim)]
        e = int(random.gauss(mu=0, sigma=3)-1)
        b = np.dot(a, sk) + m + e
        return LweCiphertext(self, a, b)

    def decrypt(self, ct: type['LweCiphertext'], sk: List[int]):
        return ct.b-np.dot(ct.a, sk)

class LweCiphertext:
    a: List[int]
    b: int
    params: type['LweParams']

    def __init__(self, params: type['LweParams'], a: List[int], b: int):
        self.params = params
        self.a = a
        self.b = b


def main(): #TODO: move to tests
    params = LweParams(modulo=10, dim=3)



if __name__ == '__main__':
    main()