#from numpy import arange, array
from pylab import *
from math import sqrt
import logging

from Physics.Radiation.Polarization import Signal
from Radio_Astronomy import dB

logging.basicConfig(level=logging.INFO)

def my_hybrid(e=0):
  """
  converts from linear to circular modes
  
  Its inverse converts from circular to linear modes.
  """
  return array([[(1+e)*(1+0j), (1-e)*(0+1j)],
                [(1+e)*(0+1j), (1-e)*(0-1j)]])/sqrt(2)
                  
mylogger = logging.getLogger()
mylogger.setLevel(logging.INFO)

imbalance = []
degree_linear = []
degree_circular = []
eccentricity = []

for e in arange(0.0, 0.7, 0.01):
  test_sig = Signal(mode="linear", E=(1,1))
  hybrid_out = my_hybrid(e=e/2).dot(test_sig.linear)
  imbalance.append(dB((1+e/2)/(1-e/2)))
  new_sig = Signal(mode="circular", E=hybrid_out)
  print "error=",e,"  Stokes:",new_sig.Stokes
  degree_linear.append(new_sig.ellipse['d_lin'])
  degree_circular.append(new_sig.ellipse['d_cir'])
  eccentricity.append(new_sig.ellipse['eccen'])

figure()
plot(imbalance,degree_linear,label="deg.lin.pol.")
plot(imbalance,degree_circular,label="deg.cir.pol")
plot(imbalance,eccentricity,label="eccentricity")
ylim(0,1)
xlabel("Power Imbalance (dB)")
legend(loc="center left")
grid()
#title()
show()
