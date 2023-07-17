import numpy as np

NUM_OF_BLANKS = 2 #6
NUM_OF_SQUARES = 4 #16

def main():
    squares = [[]]
    gen_all_squares(squares, 0)


    # 1 0 0 0
    # 1 0 0 1
    # 0 1 0 0
    # 1 0 0 0
    # square = [[1,0],
    #           [0,1]]
    for square in squares:
        if is_even_rows_cols(square):
            print(square)
        # print(is_even_rows_cols(square))


def gen_all_squares(squares, step):
    print(step)
    if step == NUM_OF_BLANKS:
        return
    # if step == 5: # optimization: if all rows and colms are even than adding one more will not.
    #    tmp_squares = squares[:] #copy
    #    print("before filtering: ", len(squares))
    #    for square in tmp_squares:
    #        if is_even_rows_cols(square):
    #            squares.remove(square)
    #    print("after filtering: ", len(squares))

    tmp_squares = squares[:]
    for k in range(len(tmp_squares)):
        squares.remove(tmp_squares[k])
        for i in range(NUM_OF_SQUARES):
            square = tmp_squares[k][:]  # copy
            max1 = -1
            if len(square)>0:
                max1=max(square)
            if  i > max1:
                square.append(i)
                squares.append(square)
    print('len(squares): ', len(squares))
    gen_all_squares(squares, step + 1)

def is_even_rows_cols(square):
    is_col = is_even_cols(square)
    is_row =  is_even_rows(square)
    return is_col and is_row

def is_even_cols(square):
    return is_even_rows(np.array(square).T.tolist())

def is_even_rows(square):
    s=[[0]*4]*4
    for i in range(4):
        for j in range(4):
            if (i+j*4) in square:
                s[i][j] = 0
            else:
                s[i][j] = 1
    for row in range(len(s)):
        cnt = 0
        for val in s[row]:
            if val == 1:
                cnt+=1
        if cnt % 2 != 0:
            return False
    return True

if __name__ == '__main__':
    main()

