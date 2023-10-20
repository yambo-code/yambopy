from yambopy import *
from yambopy.plot  import *

# Read the yambo_nl databases for all laser frequencies
NLDB=YamboNLDB()

#
# Analize the non-linear reponse using Eq.
# of Eq. 26 and 27 of Ref. Phys. Rev. B 88, 235113 (2013)
# results are stored in the files o.YamboPy-X_probe_order_x 
# where x goes from 0 to X_order
#
Harmonic_Analysis(NLDB,X_order=4)
