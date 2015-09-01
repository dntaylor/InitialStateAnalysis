PlottingUtils
=============
A comprehensize plotting framework for use with ntuples output from ISA. The style of the plots match
the latest recommendations from the [CMS Publication Committee](https://ghm.web.cern.ch/ghm/plots/).
The primary plotter,  [PlotterBase.py](./python/PlotterBase.py), accesses information about the cross section
([xsec.py](./python/xsec.py)) and plot styles ([tdrstyle.py](./python/tdrstyle.py), [CMS_lumi.py](./python/CMS_lumi.py),
and [dataStyles.py](./python/dataStyles.py)). Convenient access methods are available via [Plotter.py](./python/Plotter.py).

Limits can be plotted separately with [limits.py](./python/limits.py).
