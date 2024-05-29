import numpy as np


def main():
    # Coefficients of the polynomial (from highest to lowest degree)
    coefficients = [1, 0, 0, 0, 1, 0, 0, 1, 1]  # polynomial: x^8 +x^4+x+1

    # Find the roots of the polynomial
    roots = np.roots(coefficients)
    print("roots: ", roots)

    # Vandermonde matrix
    nttMat = vandermonde(list(roots))
    invNttMat = np.linalg.inv(nttMat)


    elm = [1,1,1,0,0,0,0,0] # some element in GF(256)
    elm_in_ntt = np.dot(nttMat, elm)
    print("elm_in_ntt:", elm_in_ntt)

    elm_new = np.dot(invNttMat, elm_in_ntt)
    print("elm_new:", elm_new)

    elm_power = elm_in_ntt
    for i in range(254):
        elm_power = elm_power*elm_in_ntt
        elm_power_new = np.dot(invNttMat, elm_power)
        elm_power_new_real = np.real(elm_power_new)
        elm_power_new_int = np.round(elm_power_new_real)
        elm_power_new_mod_2 = elm_power_new_int % 2
        elm_power = np.dot(nttMat,elm_power_new_mod_2)
        print("elm_power_new_mod_2:",i," : ", elm_power_new_mod_2)

    elm_power_new = np.dot(invNttMat, elm_power)
    print("elm_power_new:", elm_power_new)


def vandermonde(xi_list: list):
    mat = []
    for xi in xi_list:
        row = [xi**j for j in range(len(xi_list))]
        mat.append(row)
    return mat

if __name__ == '__main__':
    main()