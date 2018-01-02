from mapg.util import *

def test_binary_row_echelon():
    A = np.zeros( (5, 5) )
    assert binary_is_row_echelon(A), 'Trivial zero matrix'

    A = np.identity(2)
    assert binary_is_row_echelon(A), '2 Identity matrix'

    A = np.zeros( (3, 5) )
    A[0,:] = 1
    A[1, 2:] = 1
    A[2, 4] = 1

    assert binary_is_row_echelon(A), 'non-square 3x5'

    A = np.zeros( (5, 5) )
    for i in range(5):
        A[i, i:] = 1

    #make some coefficients zero
    A[0, 2:] = 0
    A[2, 4:] = 0
    assert binary_is_row_echelon(A), 'square 5x5'

    A = np.zeros( (3, 5) )
    A[0,0] = 1
    A[1, 2] = 1
    A[2, 4] = 1
    print('naughty one!')
    assert binary_is_row_echelon(A), '3x5 with mostly zeros'

    A = np.ones( (5,5) )
    assert not binary_is_row_echelon(A), '5x5 ones'

    A = np.identity( 4 )
    A[3, 0] = 1
    assert not binary_is_row_echelon(A), '4x4 identity off-by-one'

    A = np.zeros( (4, 5) )
    A[3, 4] = 1
    assert not binary_is_row_echelon(A), '4x5 zeroes off-by-one'

    A = np.zeros( (2, 5) )
    A[0,1] = 1
    A[1,0] = 1

    assert not binary_is_row_echelon(A), '2x5 reversed leading coefficients'

def test_binary_guass_elim_examples():
    A = np.identity(5)
    assert binary_is_row_echelon(binary_gauss_elim(A)), '5 Identity matrix'

    A = np.zeros( (3, 5) )
    A[0,:] = 1
    A[1, 2:] = 1
    A[2, 4] = 1

    assert binary_is_row_echelon(binary_gauss_elim(A)), 'non-square 3x5 in row echelon form'

    A = np.zeros( (3, 5) )
    A[0, :] = 1
    A[1, 0] = 1
    A[2, 0] = 1
    A[2, 4] = 1
    assert binary_is_row_echelon(binary_gauss_elim(A)), 'non-square 3x5'

def test_binary_solve():
    #      0
    #     / \
    #    0   0
    #   / \ / \
    #  0   0   0
    A = np.array([
        [1,1,0,1,0,0],
        [1,1,1,0,1,0],
        [1,0,1,0,0,1]
    ])
    print(A * 1)
    solutions = binary_solve(A)
    for s in solutions:
        print(s)
    assert False