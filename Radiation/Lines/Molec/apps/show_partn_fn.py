from pylab import *
import math as M
import logging

from Physics.Radiation.Lines.Molec.jpl import get_mol_metadata

tag = 17002

mylogger = logging.getLogger()
mylogger.setLevel(logging.INFO)
warnings.filterwarnings('error')

def show_partn_fn(tag):
  """
  Make a plot of the partition function as a function of temperature
  for the molecule specified by tag.
  """
  metadata = get_mol_metadata(tag)
  Q_data = metadata[2]
  Q_temps = Q_data.keys()
  Q_temps.sort()
  Q_values = []
  for temp in Q_temps:
    try:
      Q = float(Q_data[temp])
    except OverflowError:
      continue
    mylogger.debug("show_partn_fn: Q[%s] = %s", temp, Q)
    Q_values.append(Q)
  loglog(Q_temps[:len(Q_values)],Q_values)
  loglog(Q_temps[:len(Q_values)],Q_values,'.')
  #axis([9,300,30,10000])
  xlim(Q_temps[0],Q_temps[-1])
  title(metadata[0])
  xlabel('Temperature (K)',fontsize=12)
  ylabel('Partition Function')
  grid()
  show()


show_partn_fn(tag)
