from numpy.polynomial import Polynomial


def prod(vec):
    out = 1
    for v in vec:
        out = out * v
    return out

def print_coef(title, poly_coef):
    print(title, [int(num) for num in poly_coef])