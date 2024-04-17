import matplotlib.pyplot as plt
from yambopy import *

y = YamboOut('bse',save_folder='bse/SAVE')

energy = y.files['o-yambo.eps_q1_diago_bse']['E/ev[1]']
im_eps = y.files['o-yambo.eps_q1_diago_bse']['EPS-Im[2]']

plt.plot(energy,im_eps)
#plt.savefig('bse.png')
plt.show()
