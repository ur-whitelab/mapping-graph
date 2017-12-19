import numpy as np

def binary_gauss_elim(A, b=None):
    '''
    Dennis Parkinson, Marvin Wunderlich, A compact algorithm for Gaussian elimination over GF(2) implemented on highly parallel computers, In Parallel Computing, Volume 1, Issue 1, 1984, Pages 65-73, ISSN 0167-8191, https://doi.org/10.1016/S0167-8191(84)90424-1.
    (http://www.sciencedirect.com/science/article/pii/S0167819184904241)
    Keywords: Factorization; Gaussian elimination; DAP; array processors
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
            Q[j,:], Q[i,:] = Q[i,:], Q[j,:]
            # perform column addition
            for k in range(n):
                if k != j and Q[k,j]:
                    Q[k, :] ^= Q[j,:]
    return Q[:,:-1]

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
