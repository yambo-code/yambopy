from yambopy import *

#we try to use matplotlib, if not present we won't use it
try:
    # prevents crashes if no X server present
    import sys
    if "matplotlib" not in sys.modules:
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True
