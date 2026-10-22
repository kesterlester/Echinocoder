"""
The evaluator: a torch.nn.Module mapping tensors of shape (..., n, k) to (..., D).

Its conventions follow IWPC's ``Encoding`` base class (https://github.com/jjhwilkinson/iwpc) without depending on
IWPC: ``input_shape`` and ``output_shape`` are integer buffers, and the module acts on tensors of shape
``(..., *input_shape)``.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import torch
from torch import nn, Tensor

from .bursar import bursar_coefficients, bursar_coefficients_complex, bursar_size
from .slots import Slot


class InvariantEmbedding(nn.Module):
    """
    Embeds n vectors in k dimensions, invariant under the group G and under rotations, and separating orbits.

    Pipeline: X -> w (the dots and determinants, in slot order) -> {g.w : g in G} -> Bursarial coefficients
    -> the kept subset. The signed-permutation table of G and the kept subset are fixed at build time, so evaluation
    enumerates nothing.

    Build instances with ``symembed.build_invariant_embedding`` rather than directly.
    """

    def __init__(self, n: int, k: int, slots: tuple, index: list, sign: list, dtype: torch.dtype,
                 keep: list | None = None, packing=None, metric: tuple | None = None, vector_w: bool = True):
        super().__init__()
        N, m = len(index), len(slots)
        width = m if packing is None else packing.r          # length of each polynomial q_g
        if keep is None:
            keep = list(range(bursar_size(N, width)))
        self.n, self.k, self.N, self.m = n, k, N, m

        # Raw quantities are computed as [all d_ij for i <= j (row-major)] + [all e for k-subsets (lexicographic)];
        # slot_source maps each slot to its position in that raw list. With vector_w=False, w is supplied by the
        # caller (see from_w), for example by the beam-frame front end, and none of this is needed.
        d_pairs = [(i, j) for i in range(n) for j in range(i, n)] if vector_w else []
        # Determinants are computed only if w has any (they are omitted, for example, when parity is in G).
        e_subsets = list(combinations(range(n), k)) if vector_w and any(s.kind == "e" for s in slots) else []
        raw_position = {Slot("d", p): r for r, p in enumerate(d_pairs)}
        raw_position.update({Slot("e", c): len(d_pairs) + r for r, c in enumerate(e_subsets)})
        if not vector_w:
            raw_position = {s: 0 for s in slots}

        self.register_buffer("input_shape", torch.tensor([n, k], dtype=torch.int64))
        self.register_buffer("output_shape", torch.tensor([len(keep)], dtype=torch.int64))
        self.register_buffer("d_rows", torch.tensor([p[0] for p in d_pairs], dtype=torch.int64))
        self.register_buffer("d_cols", torch.tensor([p[1] for p in d_pairs], dtype=torch.int64))
        self.register_buffer("e_subsets", torch.tensor(e_subsets, dtype=torch.int64).reshape(len(e_subsets), k))
        self.register_buffer("slot_source", torch.tensor([raw_position[s] for s in slots], dtype=torch.int64))
        self.register_buffer("g_index", torch.tensor(index, dtype=torch.int64))
        self.register_buffer("g_sign", torch.tensor(sign, dtype=dtype))
        self.register_buffer("metric_diag", torch.tensor(metric if metric is not None else (1,) * k, dtype=dtype))
        self.register_buffer("keep", torch.tensor(keep, dtype=torch.int64))
        self.keeps_all = len(keep) == bursar_size(N, width)
        self.packed = packing is not None
        if self.packed:
            self.register_buffer("pack_plus", torch.tensor(packing.plus, dtype=dtype).reshape(packing.r, m))
            self.register_buffer("pack_minus", torch.tensor(packing.minus, dtype=dtype).reshape(packing.r, m))

    @property
    def is_vector_input(self) -> bool:
        return False

    @property
    def is_vector_output(self) -> bool:
        return True

    def w(self, x: Tensor) -> Tensor:
        """The slots of w, shape (..., m)."""
        gram = (x * self.metric_diag) @ x.transpose(-1, -2)               # d_ij = v_i . eta . v_j
        d_vals = gram[..., self.d_rows, self.d_cols]
        if self.e_subsets.shape[0] > 0:
            e_vals = torch.linalg.det(x[..., self.e_subsets, :])
            raw = torch.cat([d_vals, e_vals], dim=-1)
        else:
            raw = d_vals
        return raw[..., self.slot_source]

    def forward(self, x) -> Tensor:
        """
        Parameters
        ----------
        x
            Tensor (or NumPy array) of shape (..., n, k). ``x[..., i, :]`` is the vector with the i-th label of the
            GroupSpec. NumPy input is converted to a tensor of the module's dtype.

        Returns
        -------
        Tensor
            Shape (..., D).
        """
        x = self._prepare(x)
        batch_shape = x.shape[:-2]
        out = self.from_w(self.w(x).reshape(-1, self.m))
        return out.reshape(*batch_shape, out.shape[-1])

    def from_w(self, w: Tensor) -> Tensor:
        """The outputs for given slot values w, shape (B, m) -> (B, D)."""
        gw = w[:, self.g_index] * self.g_sign                                # (B, N, m)
        if self.packed:
            # Involution packing: the coefficients are real, so only the real parts are output.
            out, _ = bursar_coefficients_complex(gw @ self.pack_plus.T, gw @ self.pack_minus.T)
        else:
            out = bursar_coefficients(gw)
        if not self.keeps_all:
            out = out[:, self.keep]
        return out

    def _prepare(self, x) -> Tensor:
        """Convert NumPy input to a tensor of the module's dtype, and check the trailing shape."""
        if isinstance(x, np.ndarray):
            x = torch.as_tensor(x, dtype=self.g_sign.dtype, device=self.g_sign.device)
        if tuple(x.shape[-2:]) != (self.n, self.k):
            raise ValueError(f"expected input of shape (..., {self.n}, {self.k}), got {tuple(x.shape)}")
        return x

    def flattened(self) -> "FlattenedInvariantEmbedding":
        """A view of this embedding taking input of shape (..., n*k), component c of vector i at position i*k + c."""
        return FlattenedInvariantEmbedding(self)


class FlattenedInvariantEmbedding(nn.Module):
    """
    Wraps an InvariantEmbedding so that it takes flat input of shape (..., n*k), as needed for IWPC-style
    concatenation of vector encodings.
    """

    def __init__(self, inner: InvariantEmbedding):
        super().__init__()
        self.inner = inner
        self.register_buffer("input_shape", torch.tensor([inner.n * inner.k], dtype=torch.int64))
        self.register_buffer("output_shape", inner.output_shape.clone())

    @property
    def is_vector_input(self) -> bool:
        return True

    @property
    def is_vector_output(self) -> bool:
        return True

    def forward(self, x) -> Tensor:
        if isinstance(x, np.ndarray):
            x = torch.as_tensor(x, dtype=self.inner.g_sign.dtype, device=self.inner.g_sign.device)
        if x.shape[-1] != self.inner.n * self.inner.k:
            raise ValueError(f"expected input of shape (..., {self.inner.n * self.inner.k}), got {tuple(x.shape)}")
        return self.inner(x.reshape(*x.shape[:-1], self.inner.n, self.inner.k))
