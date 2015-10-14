# -*- coding: utf-8 -*-
from Physics import eps_0, mu_0
from math import pi, sqrt

resistivity = {'copper': 1.7e-8, 'aluminum': 2.8e-8, 'iron': 1e-7,
               'steel-electrical': 4.6e-7, 'steel-stainless': 6.9e-7,
               'gold': 2.44e-8, 'silver': 1.68e-8,
               'graphite-min': 2.5e-6, 'graphite-max': 5e-6}

permeability = {'steel-electrical': 5e-3, 'steel-stainless': 1000*mu_0,
                'steel-carbon': 8.75e-4, 'copper': mu_0,
                'aluminum': mu_0}

permittivity = {'metal': eps_0}

def skin_depth(omega, rho, mu=mu_0, eps=eps_0):
  """
  Depth of the current layer in a conductor subject to AC fields::
   J = J  exp(-d/delta)
        S
   where J  is the surface current density and delta is the skin depth.
          S
  
  Resistivity is defined so that the resistance of a bulk conductor is::
       rho
   R = --- L
        A
  where A is the cross-sectional area and L is the length.
  
  @param omega : angular frequency (rad/s)
  @type  omega : float
  
  @param mu : magnetic permeability (H/m)
  @type  mu : float
  
  @param eps : electric permittivity (F/m)
  @type  eps : float
  
  @param rho : resistivity (ohm-m)
  @type  rho : float
  
  @return: m (float)
  """
  return 1/omega/sqrt( (mu*eps/2) * (sqrt(1+(1/(rho*omega*eps))**2) -1) )

def skin_resistance(freq, rho, diam):
  """
  Resistance in a 1-m thin wire.
  
  A metal wire is assumed.
  
  @param freq : Hz
  @type  freq : float
  
  @param rho : material resistivity, ohm-m
  @type  rho : float
  
  @param diam : diameter, m
  @type  diam : float
  
  @return: ohm/m
  """
  omega = 2*pi*freq
  delta = skin_depth(omega, rho)
  return rho/(pi*(diam-delta)*delta)