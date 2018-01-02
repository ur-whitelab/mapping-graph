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
    print('you called this function')
    #now substitute to get solutions
    print(Q * 1)
    N = Q.shape[1] - 1 # - 1 for solution column
    set_mask = np.zeros(N) == 1
    row_equations = [(np.zeros(N) == 1 for _ in range(1))] # these should yield a bit vector of the solution
    for k in range(Q.shape[0] - 1, -1, -1):
        #iterate over row equations backwards
        if np.sum(Q[k, :-1]):
            # non-zero

            # find number of free bits
            free = np.sum(set_mask & Q[k, :-1])

            #create generator of solutions for this row
            def row_solution(parent_equation):
                for s in parent_equation:
                    print(s)
                    # parity in this row will change based on set bits in solution
                    # get set solution bits which are used in this row
                    # the parity of those changes the parity of this row equation
                    row_parity = bitarray_parity(s & set_mask & Q[k, :-1]) ^ Q[k, -1]
                    # enumerate free variables restricted by their parity
                    for row_value in enumerate_parity(free, row_parity):
                        # convert those to solutions
                        row_solution = expand_list(bitfield(row_value), set_mask & Q[k, :-1])
                        #combine our set solution variables with previous row
                        yield np.array(row_solution) == 1 | s
            #update set mask
            set_mask |= Q[k, :-1]

            #instantiate and append generator evaluated in this context
            if free > 0:
                row_equations.append(row_solution(row_equations[-1]))
    #generators are all chained together, so only outter one need be generated
    print(row_equations)
    solutions = [s for s in row_equations[-1]]
    print(solutions)
    return solutions

def expand_list(values, mask):
    result = [0 for _ in mask]
    j = 0
    for i,m in enumerate(mask):
        if m:
            result[i] = values[j]
            j += 1

def bitfield(n):
    '''convert an integer into a list of bits'''
    return [n >> i & 1 for i in range(n.bit_length() - 1,-1,-1)]

def enumerate_parity( nbits, set_parity ):
    '''I have no idea how to do this so I just brute force it'''
    value = 0
    if isinstance(set_parity, types.GeneratorType):
        sp = yield from set_parity
    else:
        sp = set_parity
    while value < 2**nbits:
        value += 1
        if parity(value) == sp:
            yield value

def bitarray_parity( n ):
    return np.sum(n) % 2

def parity( n ):
    parity = 0
    while n:
        parity = ~parity
        n = n & (n - 1)
    return parity