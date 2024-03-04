import numpy as np


def rot(v, k):
    n = len(v)
    k = k % n
    return np.concatenate((v[-k:], v[:-k]))

def main():
    n = 3
    values = np.array([10,20,30,11,21,31,12,22,32]) # order column first
    print(values)
    #   [1 2 3]
    #   [4 5 6]
    #   [7 8 9]
    diags = [[1,5,9],[2,6,7], [3,4,8]]
    print(diags)
    sum = np.zeros(values.size)
    for k in range(n):
        print('values=', values)
        print('diags=',diags)
        plain = np.zeros(values.size)
        for i in range(n): #diag number
            rot_diag_i = rot(diags[i], i)
            for j in range(n):
                plain[i+j*n] = rot_diag_i[j]
        print('plain=', plain)
        sum = sum + values * plain
        values = rot(values, n)
        diags = rot(diags, 1)
    matrix = np.reshape(sum, (n, n), order='F')
    print('sum=', matrix)
    matrix1 = np.array([[1,2,3], [4,5,6], [7, 8,9]])
    matrix2 = np.array([[10,20,30],[11,21,31],[12,22,32]])
    result_matrix = np.dot(matrix1, matrix2)
    print('expected=', result_matrix)



def main2():
    n = 2
    values = np.array([10,11,20,21])
    matrix = np.reshape(values, (n, n))
    print(matrix)
    diags = [[3,6],[4,5]]
    print(diags)
    sum = np.zeros(values.size)
    for k in range(n):
        print('values=', values)
        print('diags=',diags)
        plain = np.zeros(values.size)
        for i in range(n): #diag number
            rot_diag_i = rot(diags[i], i)
            for j in range(n):
                plain[i+j*n] = rot_diag_i[j]
        print('plain=', plain)
        sum = sum + values * plain
        values = rot(values, n)
        diags = rot(diags, 1)
    matrix = np.reshape(sum, (n, n), order='F')
    print('sum=', matrix)



if __name__ == '__main__':
    main()