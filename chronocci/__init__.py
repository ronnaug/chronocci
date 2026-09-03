
import sys
import types

try:
    from anndata._io.specs.methods import H5Group
except ImportError:
    import anndata
    
  
    if 'anndata._io' not in sys.modules:
        sys.modules['anndata._io'] = types.ModuleType('anndata._io')
    if 'anndata._io.specs' not in sys.modules:
        sys.modules['anndata._io.specs'] = types.ModuleType('anndata._io.specs')
    

    methods_mod = types.ModuleType('anndata._io.specs.methods')
    sys.modules['anndata._io.specs.methods'] = methods_mod
    
    # Инициализируем registry
    registry_mod = types.ModuleType('anndata._io.specs.registry')
    sys.modules['anndata._io.specs.registry'] = registry_mod
    

    class _DummyGroup: pass
    methods_mod.H5Group = _DummyGroup
    methods_mod.ZarrGroup = _DummyGroup
    methods_mod.write_basic = lambda *args, **kwargs: None
    

    class _DummyRegistry:
        def register_write(self, *args, **kwargs):

            return lambda func: func
        def register_read(self, *args, **kwargs):
            return lambda func: func

    registry_mod._REGISTRY = _DummyRegistry() 
    
    class _IOSpec:
        def __init__(self, *args, **kwargs): pass
    registry_mod.IOSpec = _IOSpec

    anndata._io = sys.modules['anndata._io']
    anndata._io.specs = sys.modules['anndata._io.specs']
    anndata._io.specs.methods = methods_mod
    anndata._io.specs.registry = registry_mod

from .pipeline import (
    run_joint_chronological_cci_pipeline,
    snoop_py 
)
from .plots import (
    plot_cell_type_relay_timeline,
    plot_signaling_bifurcation, 
    plot_signaling_streamgraph
)

__all__ = [
    "run_joint_chronological_cci_pipeline",
    "snoop_py",
    "plot_cell_type_relay_timeline",
    "plot_signaling_bifurcation",
    "plot_signaling_streamgraph"
]
