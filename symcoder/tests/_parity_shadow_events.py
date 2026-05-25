"""Auto-generated shadow event fixtures (parity-blindness regression).

Each event below has the property:

    encode(event) == encode({k: -v for k,v in event.items()})

even though those two events are NOT equivalent under SO(3) x S_n.
They ARE equivalent under O(3) x S_n.  The current Echinocoder Phase 1 +
Phase 2 encoder is therefore parity-blind on these structured inputs.

DO NOT EDIT BY HAND.  Regenerate with the recipe in the module's
docstring (see symcoder/tests/_parity_shadow_events.py source comments).
"""
from __future__ import annotations

# 31 shadow events, each a 5-particle event in R^3.

SHADOW_EVENTS = [
    # from family 'random_discrete_011'
    {
        'a': (-1.0, -1.0, -1.0),
        'b': (-1.0, -1.0, 0.0),
        'c': (1.0, -1.0, -1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (0.0, -1.0, 0.0),
    },
    # from family 'random_discrete_011'
    {
        'a': (0.0, -1.0, 1.0),
        'b': (-1.0, 1.0, 0.0),
        'c': (-1.0, -1.0, -1.0),
        'p': (0.0, -1.0, 0.0),
        'q': (0.0, 1.0, 0.0),
    },
    # from family 'random_discrete_011'
    {
        'a': (-1.0, 0.0, 1.0),
        'b': (0.0, 0.0, 1.0),
        'c': (1.0, -1.0, 1.0),
        'p': (0.0, 0.0, -1.0),
        'q': (0.0, 0.0, 1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, -1.0, 0.0),
        'b': (1.0, -1.0, 1.0),
        'c': (-1.0, -1.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, -1.0, 1.0),
        'b': (-1.0, 1.0, 1.0),
        'c': (1.0, 0.0, 0.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, -1.0),
        'b': (1.0, 1.0, 0.0),
        'c': (-1.0, -1.0, 0.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 1.0, -1.0),
        'b': (0.0, 1.0, 1.0),
        'c': (1.0, 1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, -1.0),
        'b': (1.0, 0.0, -1.0),
        'c': (0.0, 1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, -1.0, -1.0),
        'b': (0.0, 1.0, 1.0),
        'c': (1.0, -1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, -1.0),
        'b': (-1.0, -1.0, 0.0),
        'c': (1.0, -1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, 0.0, 0.0),
        'b': (0.0, 1.0, 0.0),
        'c': (0.0, -1.0, 1.0),
        'p': (1.0, 0.0, 0.0),
        'q': (-1.0, -0.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, -1.0),
        'b': (1.0, -1.0, 1.0),
        'c': (0.0, -1.0, -1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 0.0, -1.0),
        'b': (1.0, 0.0, -1.0),
        'c': (1.0, -1.0, -1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, 1.0),
        'b': (0.0, 1.0, 1.0),
        'c': (0.0, -1.0, 0.0),
        'p': (1.0, 0.0, 0.0),
        'q': (-1.0, -0.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, 1.0, -1.0),
        'b': (-1.0, 0.0, 0.0),
        'c': (1.0, 0.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, -1.0, 1.0),
        'b': (1.0, 0.0, 0.0),
        'c': (0.0, 1.0, -1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, -1.0, 0.0),
        'b': (0.0, -1.0, 1.0),
        'c': (1.0, -1.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, 0.0, -1.0),
        'b': (1.0, 1.0, 1.0),
        'c': (0.0, -1.0, 0.0),
        'p': (1.0, 0.0, 0.0),
        'q': (-1.0, -0.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, -1.0),
        'b': (1.0, 0.0, 0.0),
        'c': (0.0, 1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, 1.0, 1.0),
        'b': (1.0, 1.0, 0.0),
        'c': (1.0, 1.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, -1.0, -1.0),
        'b': (-1.0, 0.0, 0.0),
        'c': (1.0, 1.0, 0.0),
        'p': (1.0, 0.0, 0.0),
        'q': (-1.0, -0.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 0.0, 0.0),
        'b': (1.0, 1.0, -1.0),
        'c': (1.0, 1.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, 0.0),
        'b': (0.0, 1.0, 1.0),
        'c': (-1.0, -1.0, 1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, 1.0),
        'b': (1.0, 1.0, 1.0),
        'c': (1.0, 1.0, 1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, -1.0),
        'b': (0.0, -1.0, 1.0),
        'c': (1.0, 1.0, 0.0),
        'p': (1.0, 0.0, 0.0),
        'q': (-1.0, -0.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, 0.0),
        'b': (0.0, 0.0, -1.0),
        'c': (1.0, 1.0, 0.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (1.0, 1.0, -1.0),
        'b': (1.0, 1.0, -1.0),
        'c': (0.0, 1.0, -1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 1.0, -1.0),
        'b': (0.0, 1.0, -1.0),
        'c': (0.0, -1.0, 0.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, 1.0, -1.0),
        'b': (1.0, 1.0, 0.0),
        'c': (1.0, 1.0, 1.0),
        'p': (0.0, 1.0, 0.0),
        'q': (-0.0, -1.0, -0.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (-1.0, 0.0, -1.0),
        'b': (1.0, -1.0, -1.0),
        'c': (-1.0, 0.0, -1.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
    # from family 'random_discrete_antipodal'
    {
        'a': (0.0, -1.0, 1.0),
        'b': (-1.0, 1.0, 1.0),
        'c': (0.0, 0.0, 0.0),
        'p': (0.0, 0.0, 1.0),
        'q': (-0.0, -0.0, -1.0),
    },
]
