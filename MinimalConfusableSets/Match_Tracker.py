from itertools import combinations
import sympy as sp

class Match_Tracker:
    def __init__(self, M: int, match_matrix : sp.Matrix = None):
        """
        e_to_o is a dict for which the KEYS are tuples like (1,0,0,1,0) represeting and
        even vertex (even number of ones!) and the VALUES are a count of the number of times
        that vertex links to another vertex. Similarly o_to_e which is a check.
        """
        self.M = M
        self.e_present = self.o_present = 2**(M-1)
        self.e_to_o = { e:0 for e in Match_Tracker.even_or_odd_weight_tuples(M,even=True) }
        self.o_to_e = { o:0 for o in Match_Tracker.even_or_odd_weight_tuples(M,even=False) } # TODO: get rid of o_to_e in the long term, as e_to_o should be sufficient. Just here for check in short term.

        if match_matrix is not None:
            self.remember_match_matrix(match_matrix)

    def remember_match_matrix(self, match_matrix):
        n_rows, n_cols = sp.shape(match_matrix)
        assert n_cols == self.M
        for r in range(n_rows):
                row = tuple(x for x in match_matrix.row(r))
                for x in row:
                    assert int(x) == x
                self.remember_match(row)

    def get_e_tup_from(match):
        #print(f"Getting e_tup from {match}")
        return tuple(1 if x==+1 else 0 for x in match)

    def get_o_tup_from(match):
        #print(f"Getting o_tup from {match}")
        return tuple(1 if x==-1 else 0 for x in match)

    def get_z_positions_from(match):
        #print(f"Getting z_positions from {match}")
        return tuple(i for i,x in enumerate(match) if x ==0)

    def get_vertices_matched_by(match):
        """
        A match is a tuple comprising an n even number of 1s and an odd number of minus ones, others zero.

        This method yields (e_vtx, o_vtx) pairs ... for the vertices matched by 
        """
        #print(f"Getting vertices matched by {match}")
        match_count = 0

        e_tup = Match_Tracker.get_e_tup_from(match)
        o_tup = Match_Tracker.get_o_tup_from(match) # TODO: Get rid of o_tup long term. Was just a check. Not needed.
        z_positions = Match_Tracker.get_z_positions_from(match)

        #print(f"e_tup is {e_tup}")
        #print(f"o_tup is {o_tup}")
        #print(f"z_positions are {z_positions}")

        for even_zeros in True, False:
            for z_ons_and_offs in Match_Tracker.even_or_odd_weight_tuples(len(z_positions), even_zeros):

                if even_zeros:
                    e_vtx = list(e_tup)
                    o_vtx = list(o_tup)
                else:
                    e_vtx = list(o_tup)
                    o_vtx = list(e_tup)

                for i, z_on_off in enumerate(z_ons_and_offs):
                    if z_on_off:
                        assert e_vtx[ z_positions[i] ] == 0
                        assert o_vtx[ z_positions[i] ] == 0
                        e_vtx[ z_positions[i] ] = 1
                        o_vtx[ z_positions[i] ] = 1
                
                #print(f"e_vtx {e_vtx}")
                #print(f"o_vtx {o_vtx}")
                #print(f"Yielding e_vtx {tuple(e_vtx)}")
                #print(f"Yielding o_vtx {tuple(o_vtx)}")
                yield tuple(e_vtx), tuple(o_vtx)

                match_count += 1

        assert match_count == 2**len(z_positions)

    def remember_match(self, match):
        """
        A match is a tuple comprising an n even number of 1s and an odd number of minus ones, others zero.
        """
        #print("State is", self," and ", f"{self.e_present}, {self.o_present} when asked to remember {match}.")
        for e_vtx, o_vtx in Match_Tracker.get_vertices_matched_by(match):
            #print(f"   matching vertices: e,o = {e_vtx},{o_vtx}")
            if self.e_to_o[e_vtx] == 0: # No match yet exists for this vertex, so the match here will remove this vertex!
                #print("            NEW")
                self.e_present -= 1
            else:
                #print("       OLD")
                pass
       

            if self.o_to_e[o_vtx] == 0: # No match yet exists for this vertex, so the match here will remove this vertex! TODO: get rid of o_to_e long term.
                #print("            NEW")
                self.o_present -= 1
            else:
                #print("       OLD")
                pass

            self.e_to_o[e_vtx] += 1
            self.o_to_e[o_vtx] += 1 # TODO: Get rid of o_to_e long term as just check

        #print(f" {self.e_present}, {self.o_present}")
        assert self.e_present == self.o_present

    def forget_match(self, match):
        #print(f"FOG MATCH {match}")
        """
        A match is a tuple comprising an n even number of 1s and an odd number of minus ones, others zero.
        """
        for e_vtx, o_vtx in Match_Tracker.get_vertices_matched_by(match):
            assert self.e_to_o[e_vtx] > 0
            assert self.o_to_e[o_vtx] > 0

            self.e_to_o[e_vtx] -= 1
            self.o_to_e[o_vtx] -= 1 # TODO: Get rid of o_to_e long term as just check

            if self.e_to_o[e_vtx] == 0: # Last match for this vertex is removed, so this vertex returns to the fold!
                self.e_present += 1

            if self.o_to_e[o_vtx] == 0: # Last match for this vertex is removed, so this vertex returns to the fold!  get rid of o_to_e long term.
                self.o_present += 1

        assert self.e_present == self.o_present

    def even_or_odd_weight_tuples(M, even):
        """
        Generate all binary tuples of length M with an even (or odd) number of 1s.
        """
        result = []
        for k in range(0 if even else 1, M+1, 2):  # 1,3,5,... or 0,2,4,...
            for ones in combinations(range(M), k):
                tup = [0] * M
                for i in ones:
                    tup[i] = 1
                result.append(tuple(tup))
        return result

    def number_of_even_vertices_present(self):
        assert self.e_present == self.o_present
        assert 2**(self.M-1) - sum(1 for count in self.e_to_o.values() if count) == self.e_present
        assert 2**(self.M-1) - sum(1 for count in self.o_to_e.values() if count) == self.o_present
        return self.e_present

    def number_of_odd_vertices_present(self):
        assert self.e_present == self.o_present
        assert 2**(self.M-1) - sum(1 for count in self.e_to_o.values() if count) == self.e_present
        assert 2**(self.M-1) - sum(1 for count in self.o_to_e.values() if count) == self.o_present
        return self.o_present

    def __str__(self):
        ans =  (f"Match_Tracker(M={self.M},\n"
                f"e_present, o_present = {self.e_present, self.o_present},\n"
                f"e_to_o={self.e_to_o},\n"
                f"o_to_e={self.o_to_e}"
                )

        #for key,val in self.e_to_o.items():
        #    ans += f"E-to-O  {key} : {'   ' if val else ''}{val}\n"

        #for key,val in self.o_to_e.items():
        #    ans += f"O-to-E  {key} : {'   ' if val else ''}{val}\n"

        return ans

def demo():
    match_tracker = Match_Tracker(6)
    print(match_tracker)
    print("L")

if __name__ == "__main__":
    demo()
