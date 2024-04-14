from yambopy import *
from yambopy.plot  import *

# read the yambo_nl databases
NLDB=YamboNLDB()

#store the polarization in the pol array (it contains the three Cartesian directions)
pol =NLDB.Polarization[0]

#the time series of the simulation
time=NLDB.IO_TIME_points

#time when the external field was turned on
t_initial=NLDB.Efield[0]["initial_time"]

# create an empty array for the damped polarization
pol_damped=np.empty_like(pol)

# damp the polarization using a Lorentzian function for all the three directions
for i_d in range(3):
    pol_damped[i_d,:]=damp_it(pol[i_d,:],time,t_initial,damp_type='LORENTZIAN',damp_factor=0.1/ha2ev)

# Plot the polarization
Plot_Pol(time=time, pol=pol_damped, xlim=[0,55], save_file='polarization.pdf')

# Fourier transfort the polarization, divide by the external field to get \xhi(\omega) = P(\omega)/E(\omega)
# and then get \epsilon=1 + 4 \pi \xhi
Linear_Response(time=time,pol=pol_damped,efield=NLDB.Efield[0],plot=False,plot_file='eps.pdf')
