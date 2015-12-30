# -*- coding: utf-8 -*-
"""
module Physics

The module Scientific.Physics is pretty limited.  This package is intended
mainly to support astrophysical calculations, particularly as related to
radiation. This module defines some SI constants as globals and also extends
the unit objects defined in Scientific.Physics.PhysicalQuantities.

In addition, the module is the entry point for a package with various functions
used in astrophysical calculations. The package structure is::
  ConstantsExtra      - extra constants not globally defined in this module
  Radiation
    Continuum
    Lines             - functions and data for spectral lines
      recomb_lines    - functions and data for recombination line radiation
      Molec_Lines     - functions for molecular line radiation
        jpl_cat       - reading and manipulating data in the JPL Catalog
        quantum_label - decodes the quantum numbers in JPL Catalog format

The global physical constants are in the Systeme Internationale
based on Scientific.Physics.PhysicalQuantities
but extended with::
     ke     
These are defined as globals::
     alpha, fine structure constant, e^2/(4 pi eps_0 h_bar c)
     c,     speed of light in m/s
     e,     elementary charge in coulomb
     eps_0, vacuum permittivity (or electric constant), 1/mu_0 c^2 in farad/m
     eV,    electron volt, J
     G,     Newtonian constant, m^3 /kg / s^2
     h,     Planck constant, J s
     h_bar, Planck constant, J s
     k,     Boltzman constant, J/K
     k_e,   Coulomb constant, 1/(4 pi eps_0), N m^2/C^2
     k_m,   Ampere constant
     m_e,   electron mass in kg
     m_p,   proton mass
     m_P,   Planck mass
     m_u,   mass of an atomic mass unit, kg
     mu_0,  vacuum permittivity or magnetic constant, 4 pi 10^7 N A^-2
     Z_0,   vacuum impedance, sqrt(mu0/eps_0), ohm
The dictionary AMU provides the atomic mass for::
   H, He, C, N, O, Ne, Na, Xe, and Hg.

A suggested way to use this module is::
  from Physics import *
This defines 'PQ' as an alias for Scientific.Physics.PhysicalQuantities
and 'pq' as an alias for PQ.PhysicalQuantity. pq is an object with a value
and a combination of units. The physical quantity itself is a unit too.
See 'help(pq)' for the documentation.  An example::
 In [3]: pq(1,'c')
 Out[3]: PhysicalQuantity(1,'c')
 In [4]: pq(1,'c').inBaseUnits()
 Out[4]: PhysicalQuantity(299792458.0,'m/s')
 In [5]: pq(1,'c').inUnitsOf('cm/s').value
 Out[5]: 29979245800.0
 In [6]: pq(1,'c').inUnitsOf('mi/h').value
 Out[6]: 670616629.38439512

The functions 'sound_speed' and 'plasma_frequency' are also defined here.
"""

__author__    = "Tom Kuiper kuiper@jpl.nasa.gov"
__version__   = "$Revision: 1.1.1.1 $"
__date__      = "$Date: 2008/11/10 17:47 $"
__copyright__ = "Copyright (c) 2008 California Institute of Technology"
__license__   = ""

import math
from math import pi
import Scientific.Physics.PhysicalQuantities as PQ
pq = PQ.PhysicalQuantity

# Extension to physical quantity objects:

# For EM system conversions see
# http://physics.nist.gov/Pubs/SP811/appenB9.html#ELECTRICITY

# CGS-EMU
PQ._addUnit('abA',  'A*10.')      # abampere
PQ._addUnit('Bt',   'abA')        # biot
PQ._addUnit('abC',  'C*10.')      # abcoulomb
PQ._addUnit('abF',  'F*1.e-9')    # abfarad = 1.0e-9 A**2*s**4/kg/m**2
PQ._addUnit('abH',  'H*1.e9')     # abhenry = 1.0e+9 kg*m**2/A**2/s**2
PQ._addUnit('abohm','ohm*1.e9')   # abohm   = 1.0e+9 kg*m**2/A**2/s**3
PQ._addUnit('abV',  'V*1.e8')     # abvolt  = 1.0e+8 kg*m**2/A   /s**3
PQ._addUnit('G',    'T*1e4')      # gauss
PQ._addUnit('Mx',   'Wb*1e8')     # maxwell

# CGS-ESU
PQ._addUnit('statA',  'abA/29979245800')   # statampere = abA/c(CGS)
PQ._addUnit('statC',  'abC/29979245800')   # statcoulomb
PQ._addUnit('esu',    'statC')             # ESU unit charge
PQ._addUnit('Fr',     'statC')             # franklin
PQ._addUnit('statF',  'F*1.112650e-12')    # statfarad
PQ._addUnit('statH',  'H*8.987552e+11')    # stathenry
PQ._addUnit('statohm','ohm*8.987552e+11')  # statohm
PQ._addUnit('statV',  'V*2.997925e+02')    # statvolt

# More constants of proportionality

# Coulomb's Law (SI)
# ke    = (4*pi*pq(1,'eps0'))**-1
PQ._addUnit('ke',     '(4*pi*eps0)**-1')
# Ampere's Law (SI)  mu0/(4*pi)
km    = pq(1,'mu0')/(4*pi)
# Planck mass
mP = (pq(1,'hbar')*pq(1,'c')/pq(1,'Grav'))**0.5
# impedance of free space
Z0 = pq(1,'mu0')*pq(1,'c')
# debye
#D = 1e-10*pq(1,'statC')*pq(1,'cm')
PQ._addUnit('D',      '1e-10*statC*Ang')

# Dimensionless constants (no associated units)

# fine structure constant: 1/alpha = 137.03599074450398
alpha = pq(1,'e')*pq(1,'e')/pq(1,'hbar')/pq(1,'c')/pq(4*pi,'eps0')

# Globals:
  
# speed of zero mass particles
c     = pq(1,'c').inBaseUnits().value       # 299,792,458 = 3e8   m/s
# electron charge
e     = pq(1,'e').inBaseUnits().value       # 1.602176462e-19     coulomb
#     Gauss's Law
# This is the constant of proportionality between surface integral of the
# normal component of the electric field and the enclosed charge.
eps_0 = pq(1,'eps0').inBaseUnits().value    # 8.854187817e-12     farad/m
# Energy gained by an electron moving 1 cm in an electric potential of 1 V
eV    = pq(1,'eV').inBaseUnits().value      # 1.60217733e-19      J
# Law of Gravitation
G     = pq(1,'Grav').inBaseUnits().value    # 6.673e-11           m^3 /kg /s^2
h     = pq(1,'hplanck').inBaseUnits().value # 6.62606876e-34      J s
h_bar = pq(1,'hbar').inBaseUnits().value    # 1.054571596e-34     J s
k     = pq(1,'k').inBaseUnits().value       # 1.3806503e-23       J/K
#k_e   = ke.inBaseUnits().value              # 8,987,551,788       N m^2 C^-2
k_e   = pq(1,'ke').inBaseUnits().value      # 8,987,551,788       N m^2 C^-2
k_m = km.inBaseUnits().value                # 1.0000000e-07       kg*m/A**2/s**2
m_e   = pq(1,'me').inBaseUnits().value      # 9.10938188e-31      kg
m_p   = pq(1,'mp').inBaseUnits().value                # 2.1767141e-08       kg
m_P   = mP.inBaseUnits().value
m_u   = pq(1,'amu').inBaseUnits().value  # 1.66053873e-27         kg
mu_0  = pq(1,'mu0').inBaseUnits().value  # 1.2566370614359173e-06 N/A**2
Z_0   = Z0.inBaseUnits().value           # 376.7303               ohm


# Nuclear masses
AMU   = {'H':1.0081,'He':4.0028,'C':12.01161,'N':14.0067,'O':15.9994,\
         'Ne':20.179,'Na':22.98977,'Xe':131.3,'Hg':200.59}

# Proton number
Z = {'H':1,'He':2,'C':6,'N':7,'O':8,'Ne':10,'Na':11,'Xe':54,'Hg':80}

# Some general functions

def sound_speed(T,M):
    """
    Sound speed in a gas
    
    @param T : the temperature in K
    @type  T : float
    @param M : the mean molecular weight in AMU
    @type  M : float
    @return:  the sound speed in m/s (float).
    """
    return math.sqrt(k*T/(M*m_u))

def plasma_frequency(electron_density):
    """
    Radial plasma frequency (rads/sec)

    To get cyclic plasma frequency (Hz) divide by 2*pi

    Tests
    =====
    In [1]: import Physics as P
    In [2]: P.plasma_frequency(1.)/(2*pi)
    Out[2]: 8978.6637622453945
    
    If electron density is in electrons/m^3, divide by 1e6.
    
    In [3]: P.plasma_frequency(1.e-6)/(2*pi)
    Out[3]: 8.9786637622453949
    Reference
    =========
    http://www.carnicom.com/plasma1.htm

    @param electron_density : electron density in cm^{-3}
    @type  electron_density : float
    
    @return: plasma frequency in radians/sec (float)
    """
    n = electron_density*1e6 # convert from cm^-3 to m^-3
    const = e/math.sqrt(m_e*eps_0)
    return const*math.sqrt(n)


def wavelength (wn):
  """
  wavelength in m given the wavenumber in 1/m
  """
  return 1./wn

def frequency (wn):
  """
  Frequency in Hz given the wavenumber in 1/m.
  """
  return wn*Physics.c

def wavenumber(wvln):
  """
  Wave number given the wavelength

  @return: cycles per unit distance (m)
  """
  return 1./wvln


def test():
  print "Coulomb's Law in ESU:"
  print "  (ke*pq(1,'statC')**2/pq(1,'cm')**2).inUnitsOf('dyn') = ",
  print (ke*pq(1,'statC')**2/pq(1,'cm')**2).inUnitsOf('dyn')
  print "Ampere's Law in EMU:"
  print "  (km*(pq(1,'cm')*pq(1,'abA'))**2/pq(1,'cm')**2).inUnitsOf('dyn') = ",
  print (km*(pq(1,'cm')*pq(1,'abA'))**2/pq(1,'cm')**2).inUnitsOf('dyn')

if __name__ == "__main__":
  test()
