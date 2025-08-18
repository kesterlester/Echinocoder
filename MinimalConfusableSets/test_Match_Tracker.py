import Match_Tracker

Match_Tracker = Match_Tracker.Match_Tracker

def test():
    print()
    M = 5
    tracker = Match_Tracker(M)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)

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
    tracker.remember_match(m1)
    tracker.remember_match(m2)
    tracker.remember_match(m3)
    print(tracker)

    print("forgetting {m1, m2 and m3}:")
    tracker.remember_match(m3)
    tracker.remember_match(m2)
    tracker.remember_match(m1)
    print(tracker)
    assert tracker.number_of_even_vertices_present() == 16
    assert tracker.number_of_even_vertices_present() == 2**(M-1)
