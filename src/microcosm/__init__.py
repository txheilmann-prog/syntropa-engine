"""microcosm -- a black-box netlist simulator for micro-ecosystems."""
from .component import Component, Port, load_component, SCHEMA_VERSION  # noqa: F401
from .engine import run_netlist, medium_grounding, composition_validity  # noqa: F401
from .transfer import register, build_transfer, fba_analyzer, FUNCTION_REGISTRY  # noqa: F401

# Stable, versioned engine-seam API -- the public surface external callers import against.
# Bump this version when a signature or semantics change; callers can guard on it.
API_VERSION = "0.1.0"

from .community import community_transform  # noqa: F401, E402
from .emergence import super_additivity  # noqa: F401, E402
from .thermo import route_dg  # noqa: F401, E402
try:                                     # optional higher-level orchestration entry point,
    from .pipeline import *  # noqa: F401,F403,E402  (its __all__ names the entry point; absent in the open build)
except ImportError:
    pass
