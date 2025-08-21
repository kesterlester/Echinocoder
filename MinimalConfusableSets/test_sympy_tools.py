import sympy_tools as spt
import sympy as sp

def test_strip_zero_rows(): #(M: sp.Matrix) -> sp.Matrix:
    assert spt.strip_zero_rows(sp.Matrix([[0,1],[4,2]])) == sp.Matrix([[0,1],[4,2]])
    assert spt.strip_zero_rows(sp.Matrix([[0,0],[4,2]])) == sp.Matrix([[4,2]])
    assert spt.strip_zero_rows(sp.Matrix([[0,1],[0,0]])) == sp.Matrix([[0,1]])
    assert spt.strip_zero_rows(sp.Matrix([[0,0],[0,0]])) == sp.Matrix(0,2, [])

def test_max_pivot_positions_for_viable_stripped_RRE_matrix(): #(RRE_matrix_shape, k,):
    pass

def test_pivot_positions_are_all_viable_for_stripped_RRE_matrix(): #(RRE_matrix_shape, pivot_positions, k,):
    short = spt.pivot_positions_are_all_viable_for_stripped_RRE_matrix

    assert short(RRE_matrix_shape=(3,7), k=2, pivot_positions=(0,2,4))
                            # [1,0,0,0,0,0,0,],
                            # [0,0,1,0,0,0,0,],
                            # [0,0,0,0,1,0,0,],

    assert not short(RRE_matrix_shape=(3,6), k=2, pivot_positions=(0,2,4))
                            # [1,0,0,0,0,0,],
                            # [0,0,1,0,0,0,],
                            # [0,0,0,0,1,0,],

    assert short(RRE_matrix_shape=(3,7), k=2, pivot_positions=(0,2,3))
                            # [1,0,0,0,0,0,0,],
                            # [0,0,1,0,0,0,0,],
                            # [0,0,0,1,0,0,0,],

    assert not short(RRE_matrix_shape=(3,7), k=2, pivot_positions=(1,2,4))
                            # [0,1,0,0,0,0,0,],
                            # [0,0,1,0,0,0,0,],
                            # [0,0,0,0,1,0,0,],

def test_max_rows_for_viable_stripped_RRE_matrix(): #(M, k):
    short = spt.max_rows_for_viable_stripped_RRE_matrix

    assert short(M=6, k=2) == 2
                 # 01x0x.
                 # 0001xx
                 # -----------
                 # 123456

    assert short(M=7, k=2) == 3
                 # 1x0x0..
                 # 001x0x.
                 # 00001xx
                 # -----------
                 # 1234567

    assert short(M=8, k=2) == 3
                 # 01x0x0..
                 # 0001x0x.
                 # 000001xx
                 # -----------
                 # 12345678

    assert short(M=9, k=2) == 4
                 # 1x0x0.0..
                 # 001x0x0..
                 # 00001x0x.
                 # 0000001xx
                 # -----------
                 # 123456789

    assert short(M=10, k=2) == 4
                 # 01x0x0.0..
                 # 0001x0x0..
                 # 000001x0x.
                 # 00000001xx
                 # -----------
                 # 123456789a

    assert short(M=11, k=2) == 5
                 # 1x0x0.0.0..
                 # 001x0x0.0..
                 # 00001x0x0..
                 # 0000001x0x.
                 # 000000001xx
                 # -----------
                 # 123456789ab

    assert short(M=6, k=3) == 1
                 # 001xxx
                 # ------
                 # 123456

    assert short(M=7, k=3) == 2
                 # 1xx0x..
                 # 0001xxx
                 # -------
                 # 1234567

    assert short(M=8, k=3) == 2
                 # 01xx0x..
                 # 00001xxx
                 # --------
                 # 12345678

    assert short(M=9, k=3) == 2
                 # 001xx0x..
                 # 000001xxx
                 # ---------
                 # 123456789

    assert short(M=10, k=3) == 3
                 # 1xx0x.0...
                 # 0001xx0x..
                 # 0000001xxx
                 # ----------
                 # 123456789a

    assert short(M=11, k=3) == 3
                 # 01xx0x.0...
                 # 00001xx0x..
                 # 00000001xxx
                 # -----------
                 # 123456789ab

def test_some_row_causes_collapse(): #(mat: sp.Matrix, k: int):
    short = spt.some_row_causes_collapse

    assert not short(sp.Matrix([[0,2,3],[2,4,6]]), k=1) 
    assert not short(sp.Matrix([[0,2,3],[2,4,0]]), k=1) 
    assert not short(sp.Matrix([[3,2,3],[2,4,0]]), k=1) 
    assert not short(sp.Matrix([[0,2,3],[2,0,6]]), k=1) 
    assert not short(sp.Matrix([[2,2,3],[2,4,6]]), k=1) 

    assert     short(sp.Matrix([[0,2,3],[2,4,6]]), k=2) 
    assert     short(sp.Matrix([[0,2,3],[2,4,0]]), k=2) 
    assert     short(sp.Matrix([[3,2,3],[2,4,0]]), k=2) 
    assert     short(sp.Matrix([[0,2,3],[2,0,6]]), k=2) 
    assert not short(sp.Matrix([[2,2,3],[2,4,6]]), k=2) 

    assert     short(sp.Matrix([[0,2,3],[2,4,6]]), k=3) 
    assert     short(sp.Matrix([[0,2,3],[2,4,0]]), k=3) 
    assert     short(sp.Matrix([[3,2,3],[2,4,0]]), k=3) 
    assert     short(sp.Matrix([[0,2,3],[2,0,6]]), k=3) 
    assert     short(sp.Matrix([[2,2,3],[2,4,6]]), k=3) 


def test_normal_int_matrix():
    rows=3
    cols=6
    mat = spt.normal_int_matrix(rows, cols, seed=0)
    assert mat == sp.Matrix([
                           [ 0, 0,  1,  0, -1, 0],
                           [ 1, 1, -1, -1, -1, 0],
                           [-2, 0, -1, -1, -1, 0],
                        ])

    mat = spt.normal_int_matrix(rows, cols, seed=1, sigma=100)
    assert mat == sp.Matrix([
                           [ 35,  82,  33, -130, 91,  45],
                           [-54,  58,  36,   29,  3,  55],
                           [-74, -16, -48,   60,  4, -29],
                        ])

def test_rows_are_points_in_general_position():
    
    gen_pos = spt.rows_are_points_in_general_position

    assert True == gen_pos(sp.Matrix([
                      [0,0,0], # Nothing wrong with having a point at the origin.
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [2,3,4],
                      [3,4,5],
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [2,3,4], # R0
                      [4,6,8], # R1 is parallel to R0 -- but this is not what we are testing!
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [2,3,4], # R0
                      [2,3,4], # R1 is R0.
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [2,3,4], # R0
                      [1,3,4],
                      [2,9,4],
                      [2,9,8],
                      [2,3,4], # R4 is R0.
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [1,0,0,0], # R0
                      [0,1,0,0], # R1
                      [2,2,0,0], # R2 is in plane spanned by R0 and R1, but that is not what we are testing.
                      [0,0,1,0],
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [2,0,0,0], # R0
                      [0,2,0,0], # R1
                      [1,1,0,0], # R2 is on line between R0 and R1.
                      [0,0,1,0],
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [0,0,1,0],
                      [2,0,0,0], # R1
                      [0,2,0,0], # R2
                      [1,1,0,0], # R3 is on line between R1 and R2.
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [1,0,0,0],
                      [0,1,0,0],
                      [0,0,0,1],
                      [0,0,1,0],
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [3,0,0,1],
                      [0,3,0,1],
                      [0,0,3,1],
                      [1,1,1,1], # R3 is in triangle with corners at R0, R1 and R2
                   ]))

def test_rows_are_vectors_in_general_position():
    
    gen_pos = spt.rows_are_vectors_in_general_position

    assert False == gen_pos(sp.Matrix([
                      [0,0,0], # This is not a valid basis vector.
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [2,3,4],
                      [3,4,5],
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [2,3,4], # R0
                      [4,6,8], # R1 is parallel to R0.
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [1,0,0,0], # R0
                      [0,1,0,0], # R1
                      [2,3,0,0], # R2 is in plane spanned by R0 and R1.
                      [0,0,1,0],
                   ]))

    assert False == gen_pos(sp.Matrix([
                      [2,0,0,0], # R0
                      [0,2,0,0], # R1
                      [1,1,0,0], # R2 is in plane spanned by R0 and R1.
                      [0,0,1,0],
                   ]))

    assert True == gen_pos(sp.Matrix([
                      [1,0,0,0],
                      [0,1,0,0],
                      [0,0,0,1],
                      [0,0,1,0],
                   ]))

def test_int_tuple_without_zeros():

    print()

    dim=40
    vec = spt.int_tuple_without_zeros(dim) # Should have only ones and minus ones in.
    print(f"int vec without zeros {vec}")
    assert all( type(x)==int for x in vec)
    assert all( int(x)==x for x in vec)
    assert all( (x==+1 or x==-1) for x in vec)
    assert len(vec)==dim

    dim=40
    lam = 0
    vec = spt.int_tuple_without_zeros(dim, lam=lam) # Should have only ones and minus ones in.
    print(f"int vec without zeros {vec}")
    assert all( type(x)==int for x in vec)
    assert all( int(x)==x for x in vec)
    assert all( (x==+1 or x==-1) for x in vec)
    assert len(vec)==dim

    dim=40
    lam = 1
    vec = spt.int_tuple_without_zeros(dim, lam=lam) #  May have bigger magnitude numbers inside.
    print(f"int vec without zeros {vec}")
    assert all( type(x)==int for x in vec)
    assert all( int(x)==x for x in vec)
    assert all( (x>=1 or x<=-1) for x in vec)
    assert len(vec)==dim

    dim=40
    lam = 10
    vec = spt.int_tuple_without_zeros(dim, lam=lam) # Should have at least one number bigger than one, probabilistically.
    print(f"int vec without zeros {vec}")
    assert all( type(x)==int for x in vec)
    assert all( int(x)==x for x in vec)
    assert all( (x>=1 or x<=-1) for x in vec)
    assert any( (x>1 or x<-1) for x in vec) # Technically this could fail, but only with probability of order (1/10)^40 which we ought to be able to ignore!
    assert len(vec)==dim

