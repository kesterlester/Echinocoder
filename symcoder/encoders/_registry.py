"""
symcoder.encoders._registry
============================
AtomOrbitEncoderRegistry: holds registered AtomOrbitEncoders and dispatches
orbit-encoding requests to them.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ._base import AtomOrbitEncoder, OrbitSpec, EncodingCapability, EncodingResult

if TYPE_CHECKING:
    from symatom.context import Plan


class AtomOrbitEncoderRegistry:
    """
    A collection of AtomOrbitEncoders that can be queried and dispatched to.

    Typical usage
    -------------
    Build a registry once, register all available encoders, then pass it to the
    high-level manager that will use it:

        registry = AtomOrbitEncoderRegistry()
        registry.register(SortEncoder())
        registry.register(PolyEncoder())

        spec = OrbitSpec.from_flavoured_operator(fo)

        # Primary path — let the manager choose:
        candidates = registry.query_all(spec, plan)
        # candidates: [(encoder, capability), ...] sorted by registration order,
        # filtered to only those that returned can_encode=True.
        # The high-level manager inspects capabilities and calls the chosen
        # encoder directly:
        #   encoder, cap = candidates[my_choice]
        #   result = encoder.encode(spec, event, plan)

        # Convenience path — use the highest-priority capable encoder:
        result = registry.encode(spec, event, plan)

    Registration order
    ------------------
    Encoders are queried in the order they were registered.  query_all() returns
    results in that same order.  When two capable encoders have equal priority,
    best_for() returns the one registered first (stable sort).

    Extending to higher-level objects
    -----------------------------------
    This class knows nothing about the specific types of encoders it holds.  To
    build a registry for higher-level objects (e.g. atom-pair orbits, molecules),
    replicate this class with AtomOrbitEncoder replaced by the corresponding
    base class.  No modifications here are needed.
    """

    def __init__(self) -> None:
        self._encoders: list[AtomOrbitEncoder] = []

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def register(self, encoder: AtomOrbitEncoder) -> None:
        """
        Add an encoder to the registry.

        Encoders are queried in registration order.  The same encoder instance
        may be registered more than once (though that is rarely useful).
        """
        self._encoders.append(encoder)

    def __len__(self) -> int:
        """Return the number of registered encoders."""
        return len(self._encoders)

    def __iter__(self):
        """Iterate over registered encoders in registration order."""
        return iter(self._encoders)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_all(
        self,
        spec: OrbitSpec,
        plan: Plan,
    ) -> list[tuple[AtomOrbitEncoder, EncodingCapability]]:
        """
        Ask every registered encoder whether it can handle spec.

        Returns a list of (encoder, capability) pairs — in registration order —
        for every encoder that returned can_encode=True.  Encoders that returned
        can_encode=False are silently excluded.

        This is the primary method for the high-level manager: it receives the
        full set of capable encoders plus their self-reported capability metadata,
        and can apply any selection strategy it chooses.
        """
        results = []
        for enc in self._encoders:
            cap = enc.assess(spec, plan)
            if cap.can_encode:
                results.append((enc, cap))
        return results

    def best_for(
        self,
        spec: OrbitSpec,
        plan: Plan,
    ) -> tuple[AtomOrbitEncoder, EncodingCapability] | None:
        """
        Return the (encoder, capability) pair with the highest priority among
        all capable encoders, or None if no encoder can handle spec.

        Ties are broken by registration order (earlier-registered wins).

        Primarily useful for unit tests and simple pipelines.  The high-level
        manager will typically call query_all() and apply richer selection logic.
        """
        candidates = self.query_all(spec, plan)
        if not candidates:
            return None
        # max() is stable on equal keys in Python 3.  We want registration order
        # to win on ties, so iterate from the end and take the last maximum seen
        # — equivalent to a stable max on the reversed list, then reversing back.
        # Simpler: just use max with a key; Python's max is stable (first maximum wins),
        # which is exactly registration-order preference on ties.
        return max(candidates, key=lambda pair: pair[1].priority)

    # ------------------------------------------------------------------
    # Convenience dispatch
    # ------------------------------------------------------------------

    def encode(
        self,
        spec: OrbitSpec,
        event: dict,
        plan: Plan,
    ) -> EncodingResult:
        """
        Encode spec using the highest-priority capable encoder.

        Raises RuntimeError if no registered encoder can handle spec.

        This is a convenience method for testing and simple pipelines.  For
        production use, prefer calling query_all() so the manager can apply
        its own encoder-selection logic.
        """
        best = self.best_for(spec, plan)
        if best is None:
            raise RuntimeError(
                f"No registered encoder can handle orbit spec {spec!r}. "
                f"({len(self._encoders)} encoder(s) are registered.)"
            )
        encoder, _ = best
        return encoder.encode(spec, event, plan)
