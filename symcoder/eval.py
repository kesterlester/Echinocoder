from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
from symatom.atoms import ArgumentSymmetry, Operation, Atom


@dataclass(frozen=True)
class EvaluableOperation(Operation):
    """
    An Operation extended with a numerical evaluation function.

    eval_fn(vectors) -> float
        vectors : list of arrays, one per label, in atom.labels order.
        Return value : the numerical value of the operation on those vectors,
        BEFORE the atom's sign is applied.  The sign is applied by evaluate().

    eval_fn is excluded from equality and hashing so that two
    EvaluableOperations with identical symbolic properties (name, rank,
    parity, argument_symmetry) compare equal regardless of which callable
    was attached — consistent with how symatom treats Operations as purely
    symbolic objects.
    """
    eval_fn: Callable = field(compare=False, hash=False)

    def __repr__(self):
        par = f"+{self.parity}" if self.parity == 1 else str(self.parity)
        sym = self.argument_symmetry.name
        return f"EvaluableOperation('{self.name}', rank={self.rank}, parity={par}, {sym})"


def evaluate(atom: Atom, label_to_vec: dict) -> float:
    """
    Evaluate an atom numerically given a mapping from labels to vectors.

    atom.operation must be an EvaluableOperation.  The atom's sign is
    applied to the result of eval_fn, so the return value is:

        atom.sign * atom.operation.eval_fn([label_to_vec[lbl]
                                            for lbl in atom.labels])
    """
    if not isinstance(atom.operation, EvaluableOperation):
        raise TypeError(
            f"Cannot evaluate: operation '{atom.operation.name}' has no "
            f"eval_fn (it is not an EvaluableOperation)"
        )
    vectors = [label_to_vec[lbl] for lbl in atom.labels]
    return atom.operation.eval_fn(vectors) * atom.sign
