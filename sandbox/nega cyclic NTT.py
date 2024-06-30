import numpy
from numpy.fft import fft, ifft

def cyclic_polymul(p1, p2, N):
    """Multiply two integer polynomials modulo (x^N - 1).

    p1 and p2 are arrays of coefficients in degree-increasing order.
    """
    assert len(p1) == len(p2) == N
    product = fft(p1) * fft(p2)
    inverted = ifft(product)
    return numpy.round(numpy.real(inverted)).astype(numpy.int32)


def negacyclic_polymul_preimage_and_map_back(p1, p2):
    p1_preprocessed = numpy.concatenate([p1, -p1])
    p2_preprocessed = numpy.concatenate([p2, -p2])
    product = fft(p1_preprocessed) * fft(p2_preprocessed)
    inverted = ifft(product)
    rounded = numpy.round(numpy.real(inverted)).astype(p1.dtype)
    return (rounded[: p1.shape[0]] - rounded[p1.shape[0] :]) // 4

p1 = numpy.array([1,2])
p2 = numpy.array([4,5])
print(negacyclic_polymul_preimage_and_map_back(p1, p2))