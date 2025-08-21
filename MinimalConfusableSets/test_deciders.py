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
    has_True = True in votes_for_collapse
    has_False = True in votes_for_collapse
    assert not ( has_True and has_False )


