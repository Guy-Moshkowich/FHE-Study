import numpy as np


def square_element_wise(x):
    return


def main():

    # Coefficients of the polynomial (from highest to lowest degree)
    coefficients = [1, 0, 0, 0, 1, 1, 0, 1, 1]  # polynomial: x^8 +x^4+ x^3 +x+1

    # Find the roots of the polynomial
    roots = np.roots(coefficients)
    # print("roots: ", roots)

    # Vandermonde matrix
    nttMat = vandermonde(list(roots))
    invNttMat = np.linalg.inv(nttMat)


    # 254 = 32*7+30 = 32*7 + 8*3 + 4*1 + 2*1
    elm = [0,1,0,0,0,0,0,0] # some element in GF(256)
    print('elm original: ', elm)
    elm_ntt = np.dot(nttMat, elm)
    elm_pow_32 = expo_by_squaring(elm_ntt, 32)
    elm_pow_32_intt = np.dot(invNttMat, elm_pow_32)
    elm_pow_32_intt_real = np.real(elm_pow_32_intt)
    elm_pow_32_intt_int = np.round(elm_pow_32_intt_real)
    elm_pow_32_mod_2 = elm_pow_32_intt_int % 2
    print('elm_pow_32_mod_2: ',elm_pow_32_mod_2)

    elm_pow_32_fresh_ntt = np.dot(nttMat, elm_pow_32_mod_2)
    elm_pow_224_ntt = expo_by_squaring(elm_pow_32_fresh_ntt, 7)

    elm_pow_30_ntt = expo_by_squaring(elm_ntt, 30)
    elm_pow_30_intt = np.dot(invNttMat, elm_pow_30_ntt)
    elm_pow_30_intt_real = np.real(elm_pow_30_intt)
    elm_pow_30_intt_int = np.round(elm_pow_30_intt_real)
    elm_pow_30_mod_2 = elm_pow_30_intt_int % 2
    print('elm_pow_30_mod_2: ', elm_pow_30_mod_2)

    elm_pow_254_ntt = elm_pow_224_ntt * elm_pow_30_ntt
    one_ntt = elm_pow_254_ntt * elm_ntt
    one = np.dot(invNttMat, one_ntt)
    one_real = np.real(one)
    one_int = np.round(one_real)
    one_mod2 = one_int % 2
    print('one_mod2: ', one_mod2)
    # res = elm_ntt
    # for i in range(2,256):
    #     res = res * elm_ntt
    #     elm_power_new = np.dot(invNttMat, res)
    #     elm_power_new_real = np.real(elm_power_new)
    #     elm_power_new_int = np.round(elm_power_new_real)
    #     elm_power_new_mod_2 = elm_power_new_int % 2
    #     print('elm_power_new_mod_2: ', i, elm_power_new_mod_2)
    #     res = np.dot(nttMat, elm_power_new_mod_2)
    #
    # elm_power_new = np.dot(invNttMat, res)
    # print("elm_power_new:", elm_power_new)


def vandermonde(xi_list: list):
    mat = []
    for xi in xi_list:
        row = [xi**j for j in range(len(xi_list))]
        mat.append(row)
    return mat


def expo_by_squaring2(x, n, square_func):
    if n == 1:
        return x
    if (n % 2 == 0):
        return square_func(expo_by_squaring(x, n // 2,square_func))
        # return expo_by_squaring(x, n // 2) ** 2
    if (n % 2 == 1):
        return x * square_func(expo_by_squaring(x, (n - 1) // 2, square_func))
        # return x * expo_by_squaring(x, (n - 1) // 2) ** 2
    print(x,n)



def expo_by_squaring(x, n):
    if n == 1:
        return x
    if n == 0:
        return 1
    if (n % 2 == 0):
        x_new = expo_by_squaring(x, n // 2)
        return x_new*x_new
    if (n % 2 == 1):
        x_new = expo_by_squaring(x, (n - 1) // 2)
        return x * x_new*x_new


if __name__ == '__main__':
    main()