"""
Basic physical constants and formula

There is a module Scientific.Physics.PhysicalQuantities which,
according to the documentation, should provide the physical
constants with units.  However, I have not been able to find the
constants, only tools for keeping track of units.
   
So, from the 1998 compilation by the Committee on Data for
Science and Engineering (CODATA) of the International Council of
Scientific Unions (ICSU), this package provides
Universal Constants::
     c,     speed of light in m/s
     eV,    electron volt, J
     mu_0,  vacuum permittivity or magnetic constant, 4 pi 10^7 N A^-2
     eps_0, vacuum permittivity (or electric constant), 1/mu_0 c^2 in farad/m
     Z_0,   vacuum impedance, sqrt(mu0/eps_0), ohm
     G,     Newtonian constant, m^3 /kg / s^2
     k_e,   Coulomb constant, 1/(4 pi eps_0), N m^2/C^2
     k,     Boltzman constant, J/K
Empirical Constants::
     m_e,   electron mass in kg
     m_u,   mass of an atomic mass unit, kg
     e,     elementary charge in coulomb
     h,     Planck constant, J s
     h_bar, Planck constant, J s
     m_p,   proton mass
     alpha, fine structure constant, e^2/(4 pi eps_0 h_bar c)
The dictionary AMU provides the atomic mass for H, He, C, N, O,
Ne, Na, Xe, and Hg.
In addition, from Allen, Astrophysical Quantities::
     pc, parsec in cm

The modules available here are::
  continuum    - functions relating to continuum emission,
  lines        - general functions for line emission
  recomb_lines - function relating to recombination line emission
  jpl_cat      - functions for accessing and decoding the JPL Spectral Line
                 Catalog
  molec_lines  - functions for molecular line emission
The module 'radiation' is obsolete and will be removed.
"""

__author__    = "Tom Kuiper kuiper@jpl.nasa.gov"
__version__   = "$Revision: 1.1.1.1 $"
__date__      = "$Date: 2008/11/10 17:47 $"
__copyright__ = "Copyright (c) 2007 California Institute of Technology"
__license__   = ""

import math
#import jpl_cat
#import continuum
#import lines
#import radiation
#import recomb_lines
#import molec_lines

c     = 299792458       # 29,792,458 = 3e8 m/s
mu_0  = 12.566370614e7  # N/A^2
eps_0 = 8.854187817e-12 # farad/m
Z_0   = 376.730313461   # ohm
G     = 6.673e-11       # m^3 /kg / s^2
k_e   = 8987551788.0    # 8,987,551,788 = 9e9 N m^2 C^-2
k     = 1.3806503e-23   # J/K

m_e   = 9.10938188e-31  # kg
m_u   = 1.66053873e-27  # kg
e     = 1.602176462e-19 # coulomb
h     = 6.62606876e-34  # J s
h_bar = 1.054571596e-34 # J s
m_p   = 2.1767e-8       # m
alpha = 7.297352533e-3  # dimensionless

AMU   = {'H':1.0081,'He':4.0028,'C':12.01161,'N':14.0067,'O':15.9994,\
         'Ne':20.179,'Na':22.98977,'Xe':131.3,'Hg':200.59}

eV    = e               # J


def sound_speed(T,M):
    """
    given the temperature in K and the mean molecular weight
    in AMU, returns the sound speed in cm/s.
    """
    return math.sqrt(k*T/(M*M_amu))

def plasma_frequency(electron_density):
    """ 
    returns plasma frequency in radians/sec
    given electron density in cm^{-3} 
    """
    # verified: http://www.carnicom.com/plasma1.htm
    n = electron_density*1e6 # convert from cm^-3 to m^-3
    const = e/math.sqrt(m_e*eps_0)
    return const*math.sqrt(n)

