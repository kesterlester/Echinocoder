from __future__ import annotations
from symatom.atoms import Atom


def evaluate(atom: Atom, label_to_vec: dict) -> float:
    """
    Evaluate an atom numerically given a mapping from labels to vectors.

    The atom's sign is applied to the result of eval_fn, so the return
    value is:

        atom.sign * atom.operation.eval_fn([label_to_vec[lbl]
                                            for lbl in atom.labels])
    """
    vectors = [label_to_vec[lbl] for lbl in atom.labels]
    return atom.operation.eval_fn(vectors) * atom.sign
