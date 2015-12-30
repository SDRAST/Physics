"""
module Radiation

This module contains functions related to emission mechanisms. 
"""
import logging
import Physics as P
import math as m
import numpy as np
from Lines.Molec import jpl
from Lines.Molec import koln

logging.basicConfig(level=logging.WARNING)
module_logger = logging.getLogger(__name__)

# Some handy conversions

# MHz to 1/cm
MHz2invCm = (P.pq(1e6,'Hz')/P.pq(1,'c')).inUnitsOf('1/cm').value

# K per 1/cm
K_per_inv_cm = \
 (P.pq(1,'hplanck')*P.pq(1,'c')/P.pq(1,'k')).inUnitsOf('K*cm').value

def E_upper(E_lower, freq):
  """
  Energy of the upper level of a transition
  
  Given the lower state energy and the transition frequency, this returns 
  the upper state energy.
  
  @param E_lower : lower state energy in 1/cm
  @type  E_lower : float
  
  @param freq : frequency of the transition in MHz
  @type  freq : float
  
  @return: float
  """
  return E_lower + freq*MHz2invCm

def column_density(freq, TaDv, g_up, Aul):
  """
  Column density from an integrated spectral line
  
  Given the line frequency, the integrated line intensity, the upper state
  degeneracy, and the Einstein A, returns the column density in cm**-2.
  
  @param freq : line frequrncy in MHz
  @type  freq : float
  
  @param TaDv : integrated line intensity in K-km/s
  @type  TaDv : float
  
  @param g_up : upper state degeneracy
  @type  g_up : int
  
  @param Aul : Einstein A
  @type  Aul : float
  
  @return: float
  """
  MKS = m.log(8*m.pi*P.k*m.pow(freq*1e6,2)*TaDv/(P.h*P.c**3*g_up*Aul))
  cgs = MKS - m.log(1e4)
  return cgs

def brightness_temperature(source_function, optical_depth):
  """
  Brightness temperature of a gas cloud
  
  Given the source function of the medium, expressed using the Rayleigh
  Jeans approximation in K (see source_function_K), and the optical depth 
  through the medium, returns the brightness temperature assuming no 
  background continuum emission.
  
  @param source_function : in form of Rayleigh-Jeans approximation in K
  @type  source_function : float
  
  @param optical_depth :
  @type  optical_depth : float
  
  @return: float
  """
  return source_function*(1. - m.exp(-optical_depth))

def source_function_K(freq, population_ratio):
  """
  Source function for a given population ratio and frequency
  
  Given the ratio of the upper and lower populations per degenerate
  sublevel and the frequency of the transition, returns the Rayleigh-Jeans 
  approximation of the source function in K.
  
  Note that the population ratio here is not the same as returned by
  LTE_pop_ratio.
  
  @param freq : frequency of the transition in Hz
  @type  freq : float or int
  
  @param population_ratio : (n_u/g_u)/(n_l/g_l)
  @type  population_ratio : float or rational fraction
  
  @return: float
  """
  return (P.h*freq/P.k)/(1./population_ratio - 1.)

def extinction_coefficient_MKS(n_u, pop_ratio, freq, Aul, dv):
  """
  Extinction coefficient in MKS units
  
  Given the upper and lower level degeneracies, the upper and lower degenerate
  sublevel populations, the transition frequency, the Einstein A and the 
  linewidth, returns the extinction coefficient.
  
  @param n_u : upper degenerate sublevel population in cm^-3
  @type  n_u : float
  
  @param pop_ratio : n_u/n_l
  @type  pop_ratio : float
  
  @param freq : transition frequency in Hz
  @type  freq : float or int
  
  @param Aul : Einstein A in 1/s
  @type  Aul : float
  
  @param dv : linewidth in cm/s
  @type  dv : float
  
  @return: float (1/m)
  """
  factor = pow(P.c/freq,3)*Aul/8/m.pi
  pop_term = (1/pop_ratio - 1) * n_u
  return factor*pop_term/dv

def extinction_coefficient_CGS(n_u, pop_ratio, freq, Aul, dv):
  """
  Extinction coefficient in CGS units
  
  Given the upper and lower level degeneracies, the upper and lower degenerate
  sublevel populations, the transition frequency, the Einstein A and the 
  linewidth, returns the extinction coefficient.
  
  @param n_u : upper degenerate sublevel population in m^-3
  @type  n_u : float
  
  @param pop_ratio : n_u/n_l
  @type  pop_ratio : float
  
  @param freq : transition frequency in Hz
  @type  freq : float or int
  
  @param Aul : Einstein A in 1/s
  @type  Aul : float
  
  @param dv : linewidth in m/s
  @type  dv : float
  
  @return: float (1/m)
  """
  return extinction_coefficient_MKS(n_u*1e6, pop_ratio, freq, Aul, dv/100)/100
  
def optical_depth(extinc_coef,path_length):
  """
  Optical depth
  
  Given the extinction coefficient and the path length in the same
  system of units, returns the optical depth.
  
  @param extinc_coef : extinction coefficient in 1/cm
  @type  extinc_coef : float
  
  @param path_length : depth of gas cloud in cm
  @type  path_length : float
  
  @return: float
  """
  return extinc_coef * path_length

def get_partition_func(moltag, temps):
  """
  Partition function from JPL catalog
  
  Returns log10 of the partition function at the temperatures given
  in list or array 'temps' for the molecule specified by 'moltag'.  If
  'temps' has length 1, numpy treats it as a scalar and so this function
  will return a scalar.
  The JPL catalog should be used because the Koln catalog does not give the
  partition function for methanol.
  
  @param moltag : molecule ID in JPL database
  @type  moltag : int
  
  @param temps : list of temperatures in K
  @type  temps : list of int or float
  
  @return: nparray
  """
  metadata = jpl.get_mol_metadata(moltag)
  log_Q = []
  T = metadata[2].keys()
  T.sort()
  log_T = np.log10(T)
  for key in T:
    log_Q.append(metadata[2][key])
  x=np.log10(temps)
  return np.interp(x,log_T,log_Q)

def n_upper_LTE(n_species, transition, temp):
  """
  Upper state population of a transition
  
  Returns the population density of the upper state of the line whose data
  are given in the dictionary 'transition', assuming LTE at the temperature
  'temp'.
  
  @param n_species : column density of the species (cm^-3 or m^-3)
  @type  n_species : float
  
  @param transition : spectral line data
  @type  transition : dict
  
  @param temp : temperature in K
  @type  temp float or int
  
  @return: float (same units as n_species)
  """
  moltag = abs(transition['tag'])
  # switch to the JPL catalog for methanol
  if moltag == 32504:
    moltag = 32003
  # get_partition_func() expects a list and returns a list.  If the input
  # list is length 1, it returns a scalar.
  temps = []
  temps.append(temp)
  log10_Q = get_partition_func(moltag, temp)
  module_logger.info(" n_upper_LTE: log10(Q) = %f", log10_Q)
  Q = 10**log10_Q
  E_u = E_upper(transition['E lower'], transition['freq'])
  level_energy_K = K_per_inv_cm*E_u
  n_u = n_species*transition['g upper']* \
    m.exp(-level_energy_K/temp)/Q
  return n_u

def LTE_pop_ratio(transition, temp, degenerate=True):
  """
  Given the transition data as a dictionary and a temperature in K, returns the
  LTE population ratio for the transition levels.
  """
  g_upper = 2*int(transition['q upper'][0])+1
  g_lower = 2*int(transition['q lower'][0])+1
  freq = transition['freq']
  R = m.exp(-P.h*freq*1e6/P.k/temp)
  if degenerate:
    return (float(g_upper)/g_lower)*R
  else:
    return R

