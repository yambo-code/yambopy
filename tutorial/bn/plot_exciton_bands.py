from yambopy import *


p = [ [0.0, 0.0, 0.0],
      [0.5, 0.0, 0.0],
      [1./3,1./3,0.0],
      [0.0, 0.0, 0.0]]


yw = YamboExcitonWeight('o-yambo.exc_weights_at_1')
print yw

yw.plot_exciton_bs(p,8)
#yw.plot_transitions()
