import numpy as np
import types

def integer_gauss_elim(A, b=None):
    '''
    Gaussian elmination for non-square.
    '''
    if b is None:
        b = np.ones((A.shape[0], 1))
    Q = np.concatenate( (A, b), axis=1).astype(np.int)
    m,n = Q.shape
    h,k = 0,0

    #gelim from wikipedia...hope we're lucky!
    while h < m and k < n:
        # find first non-zero element to act as pivot row
        p = np.argmax(np.abs(Q[h:, k])) + h
        # check if we found a non-zero
        if Q[p, k]:
            #swap rows
            tmp = np.copy(Q[h, :])
            Q[h, :] = Q[p, :]
            Q[p, :] = tmp
            #modify rows below pivot
            for i in range(h + 1, m):
                f = Q[i, k] // Q[h, k]
                Q[i, k] = 0
                for j in range(k + 1, n):
                    Q[i, j] -= Q[k, j] * f
            h += 1
            k += 1
        else:
            k += 1
    return Q

def integer_is_row_echelon(A):
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
        j = np.argmax(np.abs(A[i, :]))
        # in zero rows?
        if A[i,j] == 0:
            break
        if j <= leading_j:
            print('leading coefficient', A)
            return False
        leading_j = j
    return True

def integer_solve(A, debug=False):
    '''List solutions of given matrix where Ax = 1'''
    Q = integer_gauss_elim(A)
    if debug:
        print('Q = ', Q * 1)
    #now substitute to get solutions
    N = Q.shape[1] - 1 # - 1 for solution column
    set_mask = np.zeros(N + 1) == 1
    set_mask[-1] = True #answer is set
    row_equations = [(np.zeros(N) == 1 for _ in range(1))] # these should yield a bit vector of the solution
    for k in range(Q.shape[0] - 1, -1, -1):
        #iterate over row equations backwards
        if np.sum(np.abs(Q[k, :-1])):
            # non-zero
            #create generator of solutions for this row
            # need to bind these values so they are passed
            def row_solution(k, parent_equation, set_mask):
                if debug:
                    print('you have called row solution on row', k,'set_mask', set_mask)
                for s in parent_equation:
                    if debug:
                        print('recievied solution', s * 1, 'in row', Q[k, :], 'will have this subset after set', Q[k, ~set_mask])
                    #get set values and add to answer
                    row_answer = Q[k, -1] - s @ Q[k, :-1]
                    #get movable unset values
                    free = (Q[k, :] != 0) & (~set_mask)
                    if debug:
                        print('I believe these are free', free)
                    if debug:
                        print('can I make', Q[k, free],'@', '?', ' == ', row_answer)
                    #get degrees of freedom
                    for i in range(2**np.sum(free)):
                        proposed = bitfield(i,np.sum(free))
                        if debug:
                            print(f'Proposed ({i}) {proposed}')
                        if Q[k, free] @ proposed == row_answer:
                            if debug:
                                print(Q[k, free],'@', proposed, ' == ', row_answer)
                            solution = np.zeros(N)
                            solution[free[:-1]] = proposed
                            yield solution + s

            #instantiate and append generator evaluated in this context
            # find number of free bits
            if np.sum(set_mask == 0) > 0:
                if debug:
                    print('adding row solution generator ', k)
                row_equations.append(row_solution(k, row_equations[-1], np.copy(set_mask)))
                #update set mask
                set_mask[:-1] |= Q[k, :-1] != 0
            else:
                if debug:
                    print('Unable to find free var for row', k)

    #generators are all chained together, so only outter one need be generated
    if debug:
        print(row_equations)
    solutions = [s for s in row_equations[-1]]
    return solutions

def bitfield(n, size):
    '''convert an integer into a list of bits'''
    expanded = [1 if digit=='1' else 0 for digit in bin(n)[2:]]
    while len(expanded) < size:
        expanded.insert(0,0)
    return expanded