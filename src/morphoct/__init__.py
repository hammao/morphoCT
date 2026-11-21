"""morphoCT: microCT analysis of rocks, core plugs and foraminifera.

Pipeline: io -> preprocess -> segment -> morphometrics -> rve -> properties ->
network -> report. Heavy lifting is delegated to morphoGen (measures, RVE statistics,
FFT and lattice Boltzmann solvers) and PoreSpy/OpenPNM (pore networks).
"""

__version__ = "0.1.0.dev0"

from morphoct.io import Volume, load_volume, save_volume  # noqa: E402,F401
