"""
Layer 2 (notes, Section 3): the Bursarial embedding, evaluated with torch tensor operations so that autograd works.

Given N polynomials q_g(x) = sum_j gw[g, j] x^j, each of degree m-1, form

    B(x, y) = prod_g (y + q_g(x)) = sum_c y^(N-c) E_c(x),     E_c(x) = sum_j E_{c,j} x^j,

and return the coefficients E_{c,j} for c = 1..N and j = 0..c(m-1), ordered by c then j (ascending powers of x,
as in Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset.py).
"""
from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import Tensor


def bursar_size(N: int, m: int) -> int:
    """The number of outputs, N + (m-1) N (N+1) / 2."""
    return N + (m - 1) * N * (N + 1) // 2


def bursar_output_labels(N: int, m: int) -> list[tuple[int, int]]:
    """The (c, j) of each output, in output order."""
    return [(c, j) for c in range(1, N + 1) for j in range(c * (m - 1) + 1)]


def _poly_mul_rows(rows: Tensor, q: Tensor) -> Tensor:
    """
    Multiply each row polynomial by that batch item's q.

    rows: (B, R, Lr) coefficients in ascending powers of x.  q: (B, m).  Returns (B, R, Lr + m - 1).
    """
    B, R, Lr = rows.shape
    m = q.shape[-1]
    weight = q.flip(-1).repeat_interleave(R, dim=0).unsqueeze(1)          # (B*R, 1, m)
    out = F.conv1d(rows.reshape(1, B * R, Lr), weight, padding=m - 1, groups=B * R)
    return out.reshape(B, R, Lr + m - 1)


def bursar_coefficients(gw: Tensor) -> Tensor:
    """
    Parameters
    ----------
    gw
        Tensor of shape (B, N, m): the coefficients of the N polynomials q_g(x), ascending powers of x.

    Returns
    -------
    Tensor
        Shape (B, D) with D = bursar_size(N, m): the E_{c,j}, ordered as in ``bursar_output_labels``.
    """
    B, N, m = gw.shape
    L = N * (m - 1) + 1
    # E[:, c, :] holds E_c(x). Start from the empty product: E_0 = 1.
    E = F.pad(torch.ones(B, 1, 1, dtype=gw.dtype, device=gw.device), (0, L - 1, 0, N))
    for g in range(N):
        # After g factors only E_0..E_g are non-zero, and E_c has degree c(m-1).
        Lg = g * (m - 1) + 1
        product = _poly_mul_rows(E[:, :g + 1, :Lg], gw[:, g, :])         # (B, g+1, Lg+m-1)
        # (y + q) * sum_c y^(N'-c) E_c  =>  E_c += q * E_{c-1}
        E = E + F.pad(product, (0, L - (Lg + m - 1), 1, N - g - 1))
    flat_index = torch.tensor([c * L + j for c, j in bursar_output_labels(N, m)], device=gw.device)
    return E.reshape(B, (N + 1) * L)[:, flat_index]
