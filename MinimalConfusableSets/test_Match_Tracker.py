import Match_Tracker
import sympy as sp

Match_Tracker = Match_Tracker.Match_Tracker

def test():
    print()
    M = 5
    tracker = Match_Tracker(M)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)

    print()
    m0 = (0,0,0,-1,0)
    print("remembering {m0}:")
    tracker.remember_match(m0)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 0 # Even is Odd if an vector is zeroed!
    assert tracker.number_of_odd_vertices_present() == 0 # Even is Odd if an vector is zeroed!
    tracker.forget_match(m0)
    assert tracker.number_of_even_vertices_present() == 16

    print()
    m1 = (-1,0,1,0,1)
    z1 = sum(1 for i in m1 if i==0)
    print("remembering {m1}:")
    tracker.remember_match(m1)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 12
    assert tracker.number_of_even_vertices_present() == 16 - 2**z1

    print()
    print("forgetting {m1}:")
    tracker.forget_match(m1)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)

    print()
    m2 = (-1,-1,-1,-1,-1)
    z2 = sum(1 for i in m2 if i==0)
    print("remembering {m2}:")
    tracker.remember_match(m2)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16 - 2**z2

    print()
    print("forgetting {m2}:")
    tracker.forget_match(m2)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)


    print()
    m3 = (-1,-1,-1,0,0)
    z3 = sum(1 for i in m3 if i==0)
    print("remembering {m3}:")
    tracker.remember_match(m3)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 12
    assert tracker.number_of_even_vertices_present() == 16 - 2**z3

    print()
    print("forgetting {m3}:")
    tracker.forget_match(m3)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)

    print("remembering {m1, m2 and m3}:")
    assert tracker.number_of_even_vertices_present() == 16
    tracker.remember_match(m1)
    assert tracker.number_of_even_vertices_present() == 12
    tracker.remember_match(m2)
    tracker.remember_match(m3)
    print(tracker)

    print("forgetting {m1, m2 and m3}:")
    tracker.forget_match(m3)
    tracker.forget_match(m2)
    assert tracker.number_of_even_vertices_present() == 12
    tracker.forget_match(m1)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16

    print("remembering {m1, m2 and m3}:")
    assert tracker.number_of_even_vertices_present() == 16
    tracker.remember_match(m1)
    assert tracker.number_of_even_vertices_present() == 12
    tracker.remember_match(m2)
    tracker.remember_match(m3)
    print(tracker)

    print("forgetting {m1, m2 and m3} in a different order:")
    tracker.forget_match(m2)
    tracker.forget_match(m1)
    tracker.forget_match(m3)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16

def test2():
    M=6
    tracker=Match_Tracker(M=M)
    m = (-1, -1, -1, -1, -1, 0)
    tracker.remember_match(m)
    print(tracker) 

def test3():
    M=6
    mat = sp.Matrix([[-1, -1, -1, -1, -1, 0]])
    tracker=Match_Tracker(M=M, match_matrix=mat)
    print(tracker) 

def test4():
    M=6
    tracker=Match_Tracker(M=M)
    m1 = (-1, -1, -1, -1, -1, 0)
    m2 = (-1, -1, 0, 0, 0, -1)
    tracker.remember_match(m1)
    tracker.remember_match(m2)
    print(tracker) 

