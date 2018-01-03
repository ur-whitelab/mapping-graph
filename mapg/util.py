import numpy as np
import types

def binary_gauss_elim(A, b=None):
    '''
    GF(2) Gaussian elimination.
    '''
    if b is None:
        b = np.ones((A.shape[0], 1))
    Q = np.concatenate( (A, b), axis=1) == 1
    n,m = A.shape
    # iterate over columns
    for j in range(min(n,m)):
        # find first non-zero element to act as pivot row
        i = np.argmax(Q[j:, j]) + j
        # check if we found a non-zero
        if Q[i, j]:
            #swap rows
            tmp = np.copy(Q[j, :])
            Q[j, :] = Q[i, :]
            Q[i, :] = tmp
            # perform column addition
            for k in range(n):
                if k != j and Q[k,j]:
                    Q[k, :] ^= Q[j,:]
    return Q

def binary_is_row_echelon(A):
    n,m = A.shape
    # check zero rows are after non-zero rows
    nonzero = True
    for i in range(n):
        if sum(A[i, :]) == 0 and nonzero:
            nonzero = False
        elif sum(A[i, :]) > 0 and not nonzero:
            print('non-zero row', A)
            return False
    # check leading coefficients are to the right of previous
    leading_j = -1
    for i in range(n):
        # assumes 0/1
        j = np.argmax(A[i, :])
        # in zero rows?
        if A[i,j] == 0:
            break
        if j <= leading_j:
            print('leading coefficient', A)
            return False
        leading_j = j
    return True

def binary_solve(A):
    '''List solutions of given matrix where Ax = 1'''
    Q = binary_gauss_elim(A)
    #now substitute to get solutions
    N = Q.shape[1] - 1 # - 1 for solution column
    set_mask = np.zeros(N) == 1
    row_equations = [(np.zeros(N) == 1 for _ in range(1))] # these should yield a bit vector of the solution
    for k in range(Q.shape[0] - 1, -1, -1):
        #iterate over row equations backwards
        if np.sum(Q[k, :-1]):
            # non-zero
            #create generator of solutions for this row
            # need to bind these values so they are passed
            def row_solution(k, parent_equation, set_mask, free):
                print('you have called row solution on row', k, 'free', free, 'set_mask', set_mask)
                for s in parent_equation:
                    print('recievied solution', s * 1, 'in row', k)
                    # parity in this row will change based on set bits in solution
                    # get set solution bits which are used in this row
                    # the parity of those changes the parity of this row equation
                    # s & ~set_mask -> solutions which are set
                    # solutions which are set & Q[k, :-1] -> solutions which are set and in equation
                    # solutions which are set and in equation ^ Q[k,-1] -> compute change to parity
                    row_parity = bitarray_parity(s & set_mask & Q[k, :-1]) ^ Q[k, -1]
                    # enumerate free variables restricted by their parity
                    for row_value in enumerate_parity(free, row_parity):
                        # convert those to solutions
                        # expand list to be only unset solutions in equation
                        solution = expand_list(bitfield(row_value), ~set_mask & Q[k, :-1])
                        #combine our set solution variables with previous row
                        print('row', k, 'row value', row_value, 'parity', row_parity, 'unmod parity',
                        Q[k, -1], 'solution', np.array(solution) == 1, 'set', set_mask * 1, 'row', Q[k, :-1] * 1, 'free', (~set_mask & Q[k, :-1]) * 1)
                        print('combined solution for ', k, 's', s, 'combined', ((np.array(solution) == 1) | s) * 1)
                        yield (np.array(solution) == 1) | s


            #instantiate and append generator evaluated in this context
            # find number of free bits
            free = np.sum(~set_mask & Q[k, :-1])
            print('before append', free, Q[k,:-1], set_mask)
            if free > 0:
                row_equations.append(row_solution(k, row_equations[-1], np.copy(set_mask), free))
                #update set mask
                print(set_mask, Q[k, :-1])
                set_mask |= Q[k, :-1]
                print(set_mask, 'after')

    #generators are all chained together, so only outter one need be generated
    print(row_equations)
    solutions = [s for s in row_equations[-1]]
    return solutions

def expand_list(values, mask):
    result = [0 for _ in mask]
    j = 0
    for i,m in enumerate(mask):
        if m:
            if j == len(values):
                result[i] = 0
            else:
                result[i] = values[j]
                j += 1
    return result

def bitfield(n):
    '''convert an integer into a list of bits'''
    return [n >> i & 1 for i in range(n.bit_length() - 1,-1,-1)]

def enumerate_parity( nbits, set_parity ):
    '''I have no idea how to do this so I just brute force it'''
    value = 0
    while value < 2**nbits:
        if parity(value) == set_parity:
            yield value
        value += 1

def bitarray_parity( n ):
    return np.sum(n) % 2

def parity( n ):
    parity = 0
    while n:
        parity = 1 ^ parity
        n = n & (n - 1)
    return parity