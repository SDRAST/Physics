"""module Radiation
This module contains functions related to emission mechanisms. 
This module also defines the physical quantity for jansky and has some
handy conversions."""

import Physics as P
import math as m
import numpy as np
from Lines.Molec import jpl
from Lines.Molec import koln

# Some handy conversions

# MHz to 1/cm
MHz2invCm = (P.pq(1e6,'Hz')/P.pq(1,'c')).inUnitsOf('1/cm').value

# K per 1/cm
K_per_inv_cm = \
 (P.pq(1,'hplanck')*P.pq(1,'c')/P.pq(1,'k')).inUnitsOf('K*cm').value

def E_upper(E_lower,freq):
  """Given the lower state energy in 1/cm and the transition frequency
  in MHz, this returns the upper state energy in 1/cm."""
  return E_lower + freq*MHz2invCm

def column_density(freq,TaDv,g_up,Aul):
  """Given the line frequency in MHz, the integrated line intensity in
  K km/s, the upper state degeneracy, and the Einstein A, returns the
                      -2
  column density in cm  ."""
  MKS = m.log(8*m.pi*P.k*m.pow(freq*1e6,2)*TaDv/(P.h*P.c**3*g_up*Aul))
  cgs = MKS - m.log(1e4)
  return cgs

def brightness_temperature(source_function,optical_depth):
  """Given the source function of the medium, expressed using the Rayleigh
  Jeans approximation in K, and the optical depth through the medium,
  returns the brightness temperature assuming no background continuum
  emission."""
  return source_function*(1. - m.exp(-optical_depth))

def source_function_K(freq,population_ratio):
  """Given the ratio of the upper and lower populations per degenerate
  sublevel, i.e., (n_u/g_u)/(n_l/g_l), and the frequency of the transition
  in Hz, returns the Rayleigh-Jeans approximation of the source function
  in K."""
  return (P.h*freq/P.k)/(1./population_ratio - 1.)

def MKS_extinction_coefficient(n_u, pop_ratio, freq, Aul, dv):
  """Given the upper and lower level degeneracies, the upper and lower
  degenerate sublevel populations in m^-3, the transition frequency in Hz,
  the Einstein A in 1/s and the linewidth in m/s, returns the extinction
  coefficient."""
  factor = pow(P.c/freq,3)*Aul/8/m.pi
  pop_term = (1/pop_ratio - 1) * n_u
  return factor*pop_term/dv

def optical_depth(extinc_coef,path_length):
  """Given the extinction coefficient and the path length in the same
  system of units, returns the optical depth."""
  return extinc_coef * path_length

def get_partition_func(moltag,temps):
  """Returns log10 of the partition function at the temperatures given
  in list or array 'temps' for the molecule specified by 'moltag'.  If
  'temps' has length 1, numpy treats it as a scalar and so this function
  will return a scalar.
  The JPL catalog should be used because the Koln catalog does not give the
  partition function for methanol."""
  metadata = jpl.get_mol_metadata(moltag)
  log_Q = []
  T = metadata[2].keys()
  T.sort()
  log_T = np.log10(T)
  for key in T:
    log_Q.append(metadata[2][key])
  x=np.log10(temps)
  return np.interp(x,log_T,log_Q)

def n_upper_LTE(n_species,transition,temp):
  """Returns the population of the upper state of the line whose data
  are given in the dictionary 'transition', assuming LTE at the temperature
  'temp'."""
  moltag = abs(transition['tag'])
  # switch to the JPL catalog for methanol
  if moltag == 32504:
    moltag = 32003
  # get_partition_func() expects a list and returns a list.  If the input
  # list is length 1, it returns a scalar.
  temps = []
  temps.append(temp)
  Qs = get_partition_func(moltag,temp)
  Q = m.pow(10.,Qs)
  E_u = E_upper(transition['E lower'],transition['freq'])
  level_energy_K = K_per_inv_cm*E_u
  n_u = n_species*transition['g upper']* \
    m.exp(-level_energy_K/temp)/Q
  return n_u

def LTE_pop_ratio(transition,temp):
  """Given the transition data as a dictionary such as returned by
  methanol.get_A_state_transitions() and a temperature in K, returns the LTE
  population ratio for the transition levels."""
  g_upper = 2*int(transition['q upper'][0])+1
  g_lower = 2*int(transition['q lower'][0])+1
  freq = transition['freq']
  return (g_upper/g_lower)*m.exp(-P.h*freq*1e6/P.k/temp)

