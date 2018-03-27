from mapg.util import *

def test_integer_row_echelon():
    A = np.zeros( (5, 5) )
    assert integer_is_row_echelon(A), 'Trivial zero matrix'

    A = np.identity(2)
    assert integer_is_row_echelon(A), '2 Identity matrix'

    A = np.zeros( (3, 5) )
    A[0,:] = 1
    A[1, 2:] = 1
    A[2, 4] = 1

    assert integer_is_row_echelon(A), 'non-square 3x5'

    A = np.zeros( (5, 5) )
    for i in range(5):
        A[i, i:] = 1

    #make some coefficients zero
    A[0, 2:] = 0
    A[2, 4:] = 0
    assert integer_is_row_echelon(A), 'square 5x5'

    A = np.zeros( (3, 5) )
    A[0,0] = 2
    A[1, 2] = 4
    A[2, 4] = 11
    assert integer_is_row_echelon(A), '3x5 with mostly zeros'

    A = np.ones( (5,5) )
    assert not integer_is_row_echelon(A), '5x5 ones'

    A = np.identity( 4 )
    A[3, 0] = 2
    assert not integer_is_row_echelon(A), '4x4 identity off-by-one'

    A = np.zeros( (4, 5) )
    A[3, 4] = 3
    assert not integer_is_row_echelon(A), '4x5 zeroes off-by-one'

    A = np.zeros( (2, 5) )
    A[0,1] = 6
    A[1,0] = 4

    assert not integer_is_row_echelon(A), '2x5 reversed leading coefficients'

def test_integer_guass_elim_examples():
    A = np.identity(5)
    assert integer_is_row_echelon(integer_gauss_elim(A)), '5 Identity matrix'

    A = np.zeros( (3, 5) )
    A[0,:] = 1
    A[1, 2:] = 1
    A[2, 4] = 1

    assert integer_is_row_echelon(integer_gauss_elim(A)), 'non-square 3x5 in row echelon form'

    A = np.zeros( (3, 5) )
    A[0, :] = 1
    A[1, 0] = 1
    A[2, 0] = 1
    A[2, 4] = 1
    assert integer_is_row_echelon(integer_gauss_elim(A)), 'non-square 3x5'

def test_integer_solve():
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
    solutions = integer_solve(A)
    for s in solutions:
        print(s * 1)
    assert False