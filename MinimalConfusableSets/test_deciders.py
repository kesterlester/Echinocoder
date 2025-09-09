
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
    
    has_True = False
    has_False = False

    for i, bat_matrix in enumerate(bat_matrices):
        decider = deciders.Rational_Decider(M=5, k=2, unscaled_bad_bat_matrix=bat_matrix)
        
        vote_for_collapse, null_space = decider.vote_for_collapse_and_null_space(L_matrix)

        print("vote_for_collapse is: ",vote_for_collapse)
        
        # Expect votes_for_collapse=(True, True, True, False) 
        expected = (True, True, True, False)
        assert vote_for_collapse == expected[i]
