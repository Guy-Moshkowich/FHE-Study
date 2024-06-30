import numpy as np

def fft(x):
    N = len(x)
    if N <= 1:
        return x
    even = fft(x[0::2])
    odd = fft(x[1::2])
    T = [np.exp(-2j * np.pi * k / N) * odd[k] for k in range(N // 2)]
    return [even[k] + T[k] for k in range(N // 2)] + [even[k] - T[k] for k in range(N // 2)]

def ifft(X):
    N = len(X)
    if N <= 1:
        return X
    even = ifft(X[0::2])
    odd =  ifft(X[1::2])
    combined = [0] * N
    for k in range(N // 2):
        t = np.exp(2j * np.pi * k / N) * odd[k]
        combined[k] = even[k] + t
        combined[k + N // 2] = even[k] - t
    return [x / 2 for x in combined]  # Normalize by N here for true inverse

# Example usage:
x = [0]*16
x[15]=1
y = [0]*16
y[1]=1

fft_x = fft(x)
fft_y = fft(y)
result = np.multiply(fft_x, fft_y)
print("FFT result:", [int(x.real) for x in ifft(result)])
