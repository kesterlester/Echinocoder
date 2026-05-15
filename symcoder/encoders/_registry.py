"""
symcoder.encoders._registry
============================
AtomOrbitEncoderRegistry: holds registered AtomOrbitEncoderFactory instances and
dispatches orbit-encoding requests to them.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ._base import AtomOrbitEncoderFactory, AtomOrbitEncoder, OrbitSpec, EncodingResult

if TYPE_CHECKING:
    from symatom.context import Plan


class AtomOrbitEncoderRegistry:
    """
    A collection of AtomOrbitEncoderFactory objects that can be queried and
    dispatched to.

    Typical usage
    -------------
        registry = AtomOrbitEncoderRegistry()
        registry.register(SortEncoderFactory())
        registry.register(PolyEncoderFactory())

        spec = OrbitSpec.from_flavoured_operator(fo)

        # Survey all factories and get back ready-to-use encoders:
        encoders = registry.query_all(spec, plan)
        # encoders: [AtomOrbitEncoder, ...] sorted by registration order,
        # flattened across all factories.  Each encoder can be inspected
        # (output_dim, priority, method_name) before committing:
        #   enc = encoders[my_choice]
        #   result = enc.encode(event)

        # Convenience: use the highest-priority encoder directly:
        result = registry.encode(spec, event, plan)

    Registration order
    ------------------
    Factories are queried in registration order.  query_all() returns results in
    that same order within each factory's output.  When two capable encoders have
    equal priority, best_for() returns the one that appears first in query_all()
    (stable sort).
    """

    def __init__(self) -> None:
        self._factories: list[AtomOrbitEncoderFactory] = []

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def register(self, factory: AtomOrbitEncoderFactory) -> None:
        """Add a factory to the registry (appended to the end)."""
        self._factories.append(factory)

    def __len__(self) -> int:
        """Return the number of registered factories."""
        return len(self._factories)

    def __iter__(self):
        """Iterate over registered factories in registration order."""
        return iter(self._factories)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_all(
        self,
        spec: OrbitSpec,
        plan: Plan,
    ) -> list[AtomOrbitEncoder]:
        """
        Ask every registered factory for encoders that can handle spec.

        Returns a flat list of AtomOrbitEncoder instances — in factory-registration
        order, with each factory's own encoders in the order it returned them.
        Factories that returned [] are silently skipped.

        Every returned encoder is ready to use: call enc.encode(event) directly.
        """
        results: list[AtomOrbitEncoder] = []
        for factory in self._factories:
            results.extend(factory.assess(spec, plan))
        return results

    def best_for(
        self,
        spec: OrbitSpec,
        plan: Plan,
    ) -> AtomOrbitEncoder | None:
        """
        Return the highest-priority encoder among all candidates, or None.

        Ties are broken by position in query_all() (earlier wins), which
        corresponds to factory registration order.
        """
        candidates = self.query_all(spec, plan)
        if not candidates:
            return None
        return max(candidates, key=lambda e: e.priority)

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
        Encode spec using the highest-priority available encoder.

        Raises RuntimeError if no registered factory can handle spec.
        For production pipelines prefer query_all() so the manager can apply
        its own selection logic.
        """
        enc = self.best_for(spec, plan)
        if enc is None:
            raise RuntimeError(
                f"No registered factory produced an encoder for {spec!r}. "
                f"({len(self._factories)} factory/factories registered.)"
            )
        return enc.encode(event)
