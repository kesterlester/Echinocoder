def test_M5_k2():

    import deciders
    from sympy import Matrix
    
    L_matrix = Matrix([
    [-1, 0, 0,  1, 1],
    [ 0, 0, 1, -1, 1]])
    bat_matrices = [Matrix([
    [-3, -2],
    [ 1, -1],
    [ 1,  2],
    [-3,  2],
    [-1,  4]]), Matrix([
    [ 2,  1],
    [-2,  2],
    [ 2,  0],
    [ 1, -3],
    [ 3,  4]]), Matrix([
    [ 0, -1],
    [ 1,  3],
    [-2, -1],
    [ 1,  1],
    [ 1, -1]]), Matrix([
    [ 2, -2],
    [-1,  0],
    [ 0, -1],
    [-1,  3],
    [-1, -1]])]
    decider = deciders.Rational_Decider(M=5, k=2, bat_matrices=bat_matrices)
    
    votes_for_collapse, null_spaces = decider.votes_for_collapse_and_null_spaces(L_matrix)
    null_space_dimensions = tuple(len(null_space) for null_space in null_spaces)

    print("votes_for_collapse are: ",votes_for_collapse)
    
    # Expect votes_for_collapse=(True, True, True, False) 
    
    has_True = True in votes_for_collapse
    has_False = True in votes_for_collapse
    
    assert  has_True and has_False 
    

def test_worse():
        ##################################
    import deciders
    from sympy import Matrix
    
    L_matrix = Matrix([
[-1, -1, -1, 0, 0],
[-1,  0,  1, 0, 1]])
    bat_matrices = [Matrix([
[-3, -2],
[ 1, -1],
[ 1,  2],
[-3,  2],
[-1,  4]]), Matrix([
[ 2,  1],
[-2,  2],
[ 2,  0],
[ 1, -3],
[ 3,  4]]), Matrix([
[ 0, -1],
[ 1,  3],
[-2, -1],
[ 1,  1],
[ 1, -1]]), Matrix([
[ 2, -2],
[-1,  0],
[ 0, -1],
[-1,  3],
[-1, -1]]), Matrix([
[-1, -3],
[ 0, -2],
[-1,  1],
[ 1,  0],
[-1, -1]]), Matrix([
[ 1,  0],
[ 0, -2],
[ 4,  4],
[-1,  2],
[ 2,  1]]), Matrix([
[-1, 2],
[-2, 3],
[-3, 1],
[-2, 0],
[ 1, 1]]), Matrix([
[ 3,  1],
[-1,  0],
[-1, -1],
[ 1,  2],
[ 2, -2]]), Matrix([
[ 2,  0],
[ 0, -1],
[ 2,  1],
[-1,  3],
[ 2, -1]]), Matrix([
[-2, -2],
[-2, -1],
[ 1, -1],
[ 3, -2],
[ 3, -4]]), Matrix([
[-1,  1],
[ 0, -2],
[-4,  2],
[ 1, -2],
[ 1, -3]]), Matrix([
[-1, -1],
[-1,  1],
[-1, -3],
[ 1,  2],
[ 0, -1]]), Matrix([
[-1,  0],
[ 1, -1],
[ 2,  2],
[ 1,  2],
[ 0, -1]]), Matrix([
[-4, -1],
[ 1,  0],
[ 3,  1],
[ 1, -1],
[ 1,  1]]), Matrix([
[ 1, -3],
[-1, -2],
[ 0, -1],
[ 1,  0],
[-2, -1]]), Matrix([
[-2,  1],
[-2,  5],
[-2,  0],
[-1,  1],
[ 0, -1]]), Matrix([
[ 1,  1],
[ 2, -3],
[-3,  1],
[-4, -3],
[-1,  0]]), Matrix([
[ 2, -2],
[-2,  1],
[ 2,  0],
[ 1,  1],
[ 3, -2]]), Matrix([
[-2,  1],
[-2, -2],
[-5, -1],
[-3, -1],
[-2, -1]]), Matrix([
[-3,  2],
[-2, -1],
[ 0, -3],
[ 1, -1],
[-2,  0]]), Matrix([
[ 3, -3],
[-2,  1],
[ 2, -3],
[-2, -1],
[ 0, -4]]), Matrix([
[-1,  0],
[-2, -1],
[-1,  1],
[ 2, -1],
[ 1,  2]]), Matrix([
[-2, -2],
[-2, -4],
[ 1,  0],
[ 2, -1],
[-1,  1]]), Matrix([
[ 3,  1],
[-2,  1],
[ 0, -1],
[ 1,  0],
[-1,  1]]), Matrix([
[ 0,  1],
[ 1, -2],
[ 1,  0],
[-2, -1],
[ 1,  1]]), Matrix([
[ 1,  2],
[ 0,  1],
[-2,  3],
[-3,  1],
[-2, -1]]), Matrix([
[-2, -1],
[ 1,  2],
[-2,  0],
[-3,  3],
[ 0, -1]]), Matrix([
[-2,  2],
[ 3,  1],
[ 0, -3],
[ 2,  3],
[ 3, -2]]), Matrix([
[-1, -2],
[-2,  0],
[-2, -1],
[-1, -1],
[ 0, -2]]), Matrix([
[-1, -4],
[-2, -4],
[ 2, -2],
[-3,  1],
[-1,  0]])]
    decider = deciders.Rational_Decider(M=5, k=2, bat_matrices=bat_matrices)
    
    votes_for_collapse, null_spaces = decider.votes_for_collapse_and_null_spaces(L_matrix)
    null_space_dimensions = tuple(len(null_space) for null_space in null_spaces)

    print("votes_for_collapse are: ",votes_for_collapse)
    # Expect votes_for_collapse=(True, True, True, True, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True) 

    print("null space lengths are: ",null_space_dimensions)
    # Expect null_space_dimensions=(1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) 
    
    print("null spaces are: ",null_spaces)
    """
Expect null_spaces=([Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]]), Matrix([
[ 1/2],
[  -1],
[-1/2],
[   0],
[   1]])], [Matrix([
[0],
[0],
[0],
[1],
[0]]), Matrix([
[   1],
[-1/2],
[-1/4],
[   0],
[   1]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])], [Matrix([
[0],
[0],
[0],
[1],
[0]])]) 
    """    
    has_True = True in votes_for_collapse
    has_False = True in votes_for_collapse

    null_spaces_have_different_dimensions = max(null_space_dimensions) != min(null_space_dimensions)

    assert not null_spaces_have_different_dimensions
    ##################################
    

