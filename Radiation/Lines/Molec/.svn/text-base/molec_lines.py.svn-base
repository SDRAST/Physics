"""molec_lines - functions related to molecular line radiation

Tom Kuiper, 2008 Nov 10."""

from Physics import *
import math as M

DEBUG = False

def E_K(E_inv_cm):
  """Converts energy in 1/cm to K"""
  E_hz = E_inv_cm*c # (1/cm)*(cm/s)
  E_ergs = h*E_hz   # ergs
  return E_ergs/k   # K

def equilibrium_Boltzman_ratio(g_1,E_1,g_2,E_2,T):
  """For two states 1 and 2, given the degeneracies g_1 and g_2
  and energies E_1 and E_2 in 1/cm, returns the equilibrium
  population ratio of state 1 over state 2 at a temperature T."""
  delta_E = E_1-E_2
  if DEBUG:
    print "energy difference =",delta_E,"1/cm"
    print "                  =",c*delta_E,"hz"
    print "                  =",h*c*delta_E,"ergs"
    print "                  =",h*c*delta_E/k,"K"
  return (g_1/g_2)*M.exp(-delta_E/T)

