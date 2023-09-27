import numpy as np
import unittest
import random


def mat_sum(A, B):
    n = len(A[0])
    C = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = A[i][j]+B[i][j]
    return C


def mat_mul(A, B):
    n = len(A[0])
    if n == 1:
        return [[A[0][0]*B[0][0]]]
    else:
        A_1_1 = [A[i][0:n // 2] for i in range(n // 2)]
        A_1_2 = [A[i][n // 2:n] for i in range(n // 2)]
        A_2_1 = [A[i][0:n // 2] for i in range(n // 2, n)]
        A_2_2 = [A[i][n // 2:n] for i in range(n // 2, n)]

        B_1_1 = [B[i][0:n // 2] for i in range(n // 2)]
        B_1_2 = [B[i][n // 2:n] for i in range(n // 2)]
        B_2_1 = [B[i][0:n // 2] for i in range(n // 2, n)]
        B_2_2 = [B[i][n // 2:n] for i in range(n // 2, n)]

        # A_1_1 A_1_2       B_1_1 B_1_2       A_1_1*B_1_1 + A_1_2*B_2_1
        # A_2_1 A_2_2   *   B_2_1 B_2_2   = ...

        res_1_1 = mat_sum(mat_mul(A_1_1, B_1_1), mat_mul(A_1_2, B_2_1))
        res_1_2 = mat_sum(mat_mul(A_1_1, B_1_2), mat_mul(A_1_2, B_2_2))
        res_2_1 = mat_sum(mat_mul(A_2_1, B_1_1), mat_mul(A_2_2, B_2_1))
        res_2_2 = mat_sum(mat_mul(A_2_1, B_1_2), mat_mul(A_2_2, B_2_2))
        res = []
        for i in range(n//2):
            res_1_1[i].extend(res_1_2[i])
            res.append(res_1_1[i])
        for i in range(n//2):
            res_2_1[i].extend(res_2_2[i])
            res.append(res_2_1[i])
        return res

    return mat_res


class TestMatrixMultiplication(unittest.TestCase):

    def test_sum_mat(self):
        # |1 2|   |5 6|     | 19 22 |
        # |3 4| * |7 8|  =  | 43 50 |
        mat1 = [[1, 2], [3, 4]]
        mat2 = [[5, 6], [7, 8]]
        np.testing.assert_array_equal([[6,8], [10, 12]], mat_sum(mat1, mat2))

    def test_matrix_multiplication_numpy(self):
        # |1 2|   |5 6|     | 19 22 |
        # |3 4| * |7 8|  =  | 43 50 |
        mat1 = [[1,2],[3,4]]
        mat2 = [[5,6],[7,8]]
        np.testing.assert_array_equal([[19,22],[43,50]], np.dot(mat1,mat2))

    def test_matrix_multiplication_recursion_one_step(self):
        # |1 2|   |5 6|     | 19 22 |
        # |3 4| * |7 8|  =  | 43 50 |
        mat1 = [[1,2],[3,4]]
        mat2 = [[5,6],[7,8]]

        np.testing.assert_array_equal([[19,22],[43,50]], mat_mul(mat1, mat2))


    def test_matrix_multiplication_recursion_two_step(self):
        A = np.zeros((8,8))
        B = np.zeros((8,8))

        for i in range(8):
            for j in range(8):
                A[i][j] = random.randint(1, 10)
                B[i][j] = random.randint(1, 10)

        np.testing.assert_array_equal(np.dot(A,B), mat_mul(A, B))

