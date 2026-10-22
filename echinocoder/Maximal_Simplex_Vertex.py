from dataclasses import dataclass,  field
from Eji import Eji

@dataclass
class Maximal_Simplex_Vertex:
    _vertex_set: set[Eji] = field(default_factory=set)

    def __len__(self) -> int:
        return len(self._vertex_set)

    def __iter__(self):
        return iter(self._vertex_set)

    def get_canonical_form(self):
        """Mod out by Sn for this single vertex, ignoring any others."""
        # Method: sort the Eji by the i index, then populate the j's in order.
        sorted_eji_list = sorted(list(self._vertex_set), key=lambda eji: eji.i)
        renumbered_eji_list = [ Eji(j=j, i=eji.i) for j,eji in enumerate(sorted_eji_list)]
        return Maximal_Simplex_Vertex(set(renumbered_eji_list))

    def check_valid(self):
        # every j index in the set must appear at most once
        j_vals = { eji.j for eji in self._vertex_set }
        assert len(j_vals) == len(self._vertex_set)

    def get_permuted_by(self, perm):
        return Maximal_Simplex_Vertex({Eji(perm[eji.j], eji.i) for eji in self._vertex_set})
