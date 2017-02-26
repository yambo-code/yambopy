from yambopy import *

#we try to use matplotlib, if not present we won't use it
try:
    import matplotlib
    matplotlib.use('Agg') # prevents crashes if no X server present
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True
