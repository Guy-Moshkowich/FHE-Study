import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

def plot_matrix(matrix):
    cmap_colors = ['white', 'black']  # Define colors for 0 and 1 respectively
    cmap = ListedColormap(cmap_colors)

    # Convert the matrix to a numpy array for plotting
    matrix_array = np.array(matrix)

    plt.imshow(matrix_array, cmap=cmap, interpolation='nearest')
    plt.title('Matrix Plot (Custom Colors)')
    plt.show()

def print_matrix(matrix):
    for row in matrix:
        for element in row:
            print(element, end=" ")  # Using end=" " to print elements on the same line
        print()  # Move to the next line after printing a row

def print_square_matrix(matrix):
    if len(matrix) != len(matrix[0]):
        print("Not a square matrix")
        return

    max_length = max(len(str(element)) for row in matrix for element in row)
    for row in matrix:
        for element in row:
            print(f"{element:>{max_length}}", end=" ")
        print()

def switch_rows(matrix, row_index1, row_index2):
    # Check if the row indices are within the matrix size
    if row_index1 < 0 or row_index1 >= len(matrix) or row_index2 < 0 or row_index2 >= len(matrix):
        print("Invalid row indices")
        return

    # Perform row switch
    matrix[row_index1], matrix[row_index2] = matrix[row_index2], matrix[row_index1]
    return matrix


def rotate_matrix_by_column(matrix, column_index):
    if column_index < 0 or column_index >= len(matrix[0]):
        print("Invalid column index")
        return matrix

    # Transpose the matrix (rows become columns)
    transposed_matrix = list(zip(*matrix))

    # Rotate the specified column
    transposed_matrix[column_index] = list(reversed(transposed_matrix[column_index]))

    # Transpose the rotated matrix back (columns become rows)
    rotated_matrix = list(zip(*transposed_matrix))

    return rotated_matrix
def rotate_matrix_by_column(matrix, column_index):
    if column_index < 0 or column_index >= len(matrix[0]):
        print("Invalid column index")
        return matrix

    num_rows = len(matrix)
    rotated_matrix = [list(row) for row in matrix]  # Create a copy of the original matrix

    for i in range(num_rows):
        rotated_matrix[i][column_index] = matrix[(i + 1) % num_rows][column_index]

    return rotated_matrix



d = 4
n = d*d
mat = []
for i in range(n):
    row = []
    for j in range(n):
        if i == j:
            row.append(1)
        else:
            row.append(0)
    mat.append(row)

for i in range(d):
    switch_rows(mat,d*i,d*i+d-1)

for i in range(n):
    if (i %d) != 0 and (i+1) %d != 0:
        mat = rotate_matrix_by_column(mat,i)


# Calling the function to print the matrix
plot_matrix(mat)