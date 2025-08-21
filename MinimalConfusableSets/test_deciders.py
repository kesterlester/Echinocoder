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
    
    votes_for_collapse = decider.votes_for_collapse(L_matrix)
    print("votes_for_collapse are: ",votes_for_collapse)
    
    # Expect votes_for_collapse=(True, True, True, False) 
    
    has_True = True in votes_for_collapse
    has_False = True in votes_for_collapse
    
    assert not ( has_True and has_False )
    
    # The above assertion is not strictly KNOWN to be true, but rather
    # reflects uncertainty as to how the algorithm is supposed to work.
    # I.e. I currently "feel" like it should be true, but see that it is not
    # always true, and don't understandy why it is not always true.


