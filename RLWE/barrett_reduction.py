def almost_equal(a, b, tolerance=1e-9):
    return abs(a - b) < tolerance


class Barrett:
    def __init__(self, modulo: int):
        self.k = 64  # for 64 bits machine
        self.m = round(2**self.k/modulo)
        self.modulo = modulo
        assert almost_equal(self.m/(2**self.k), 1/modulo, 1/2**self.k)

    def mod(self, x: int) -> int:
        assert(x < self.modulo**2)
        r = x - ((x*self.m) >> self.k)*self.modulo
        if r >= self.modulo:
            r = r - self.modulo
        return r
