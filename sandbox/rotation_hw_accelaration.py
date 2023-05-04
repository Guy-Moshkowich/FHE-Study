import Utils.utils
from numpy.polynomial import Polynomial


def main():
    n = 4

    # define vector of slots
    slots = [20, 10, 10, 20]
    print("slots: ", slots)

    # encode by inv_special_FFT
    plaintext = Utils.utils.inv_special_fft(n, slots)
    plaintext = Polynomial([int(c.real) for c in plaintext.coef])
    print("plaintext: ", plaintext)
    sfft = Utils.utils.special_fft(n, plaintext)
    print("(sfft): ", sfft)

    # encode by NTT
    p = 73
    ntt = Utils.utils.ntt(n, p, list(plaintext.coef))
    print('ntt:', ntt)
    print('(inv_ntt): ', Utils.utils.inv_ntt(n, p, ntt))

    # rotate to the left by 1 step
    rot_ntt = [0]*n
    for i in range(0, n//2):
        rot_ntt[i] = ntt[i+1]
        rot_ntt[n//2-1] = ntt[0]

    for i in range(0, n//2):
        rot_ntt[n-i-1] = ntt[n-i-2]
        rot_ntt[n//2] = ntt[n-1]
    print('rot_ntt: ', rot_ntt)

    # decode by inv_NTT
    inv_rot_ntt = Utils.utils.inv_ntt(n, p, rot_ntt)
    print('inv_rot_ntt: ', inv_rot_ntt)
    # inv_rot_ntt = [15.0, -3.0, 0.0, 3.0] # remove. for testing.
    # decode by FFT
    fft_inv_rot_ntt = Utils.utils.special_fft(n, Polynomial(inv_rot_ntt))
    print('fft_inv_rot_ntt: ', fft_inv_rot_ntt)

    # verify that v was rotated to the left by 1 step

if __name__ == '__main__':
        main()