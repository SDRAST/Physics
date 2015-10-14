# -*- coding: utf-8 -*-
"""
Functions relation to recombination line emission
by Tom Kuiper 2008 Nov 10"""

import math
import Physics

# Rydberg constant, R_inf, in 1/m
R_inf = Physics.m_e \
        *math.pow(Physics.k_e,2) \
	       *math.pow(Physics.e,4)/2 \
	       /math.pow(Physics.h_bar,2) \
        /Physics.h/Physics.c

def recomb_freq(n, delta_n, Z, species):
    """
    Recombination line frequency Gordon, eq. 3-20::
     n -       lower state quantum number
     delta_n - change in the principal quantum number
     Z -       effective nuclear charge, 1 for singly
               ionized, 2 for doubly ionized, etc.
     species - as a string: 'H', 'He', 'C'
    Verified Lilley and Palmer, Ap.J.Suppl. 16, 143 (1968)
    """
    m_0 = Physics.m_u * Physics.AMU[species]
    R = R_inf/(1 + (Physics.m_e/m_0))
    # print "R=",R
    wn = R*Z**2*(1./n**2 - 1./(n+delta_n)**2)
    # print "k=",wn
    result = wn*Physics.c
    return result

def E_n(Z,n):
    """
    Given the effective nuclear charge and the orbital
    quantum number n, returns the energy of the level in
    eV
    """
    return R_inf*math.pow(Z,2)* \
           Physics.h*Physics.c/Physics.eV
    
def wavenumber (Z, n_f, n_i):
    """
    Given exposed nuclear charge, the final n and the initial n, returns the 
    wavenumber in 1/m.
    
    DO NOT USE FOR RADIO FREQUENCIES. It's accurate enough for optical
    recombination lines. For example           ::
       >>> Physics.wavelength(wavenumber(1,2,3))*1e9
       656.112
       which is the wavelength of H alpha in nm.
       Similarly, for H beta
       >>> Physics.wavelength(wavenumber(1,2,4))*1e9
       486.009
    """
    rydberg_constant = R_inf*math.pow(Z,2)
    paren = 1./math.pow(n_f,2) - 1./math.pow(n_i,2)
    return rydberg_constant*paren

def nearest_recomb_line(f,delta_n,Z):
    """
    Given a frequency in Hz, returns the upper orbital
    quantum number of the nearest hydrogen recombination
    line. In the radio range for delta_n << n, Gordon, 
    eq. 3.21.
    """
    m_0 = Physics.m_u*Physics.AMU['H']
    R = R_inf/(1 + (Physics.m_e/m_0))
    n_cubed = 2*R*Physics.c*Z**2*delta_n/f
    # The nearest integer
    result =\
      int(math.floor(math.pow(n_cubed,0.33333)+0.5))
    return result

def recomb_freq_interval(n,delta_n,Z):
    """
    Frequency spacing between recombination lines.
    In the radio domain; Gordon, eq. 3.23.
    """
    m_0 = Physics.m_u*Physics.AMU['H']
    R = R_inf/(1 + (Physics.m_e/m_0))
    result = 6*R*Physics.c*Z**2*delta_n/math.pow(n,4)
    return result


def F_hydro(R):
    """
    fraction of all free-free emissions due to
    hydrogen. R is the ratio of the number of heliums to
    the number of hydrogens, which is approximately
    0.08
    """
    result = 1 - R
    return result


def H_osc_str(n,delta_n):
    """
    Oscillator strengths for hydrogen for n=50-150,
    delta_n=1-5. Taken from Menzel, ApJSuppl 161, 221
    (1970)
    """
    # The intercept for the formulas below are fit quite well with
    # 10*n**-2.85
    if delta_n == 1:
        result = 9.8 + 9.6*(n-50)/50.
    elif delta_n == 2:
        result = 1.4 + 1.31*(n-50)/50.
    elif delta_n == 3:
        result = 0.44 + 0.41*(n-50)/50.
    elif delta_n == 4:
        result = 0.195 + 0.175*(n-50)/50.
    elif delta_n == 5:
        result = 0.105 + 0.090*(n-50)/50.
    else:
        result = -1
    return result


def H_line_to_cont_ratio(n,delta_n,Te):
    """
    Hydrogen recombination line to continuum ratio.
    Gordon, eq. 3.49.
    """
    fnn = H_osc_str(n,delta_n)
    freq = recomb_freq(n,1,1,'H')/1e9
    F = F_hydro(0.08)
    result = \
       1.299e5*delta_n \
              *(fnn/n) \
              *math.pow(freq,2.1) \
              *math.pow(Te,-1.15)/F \
              *math.exp(1.579e5/Te/n**2)
    return result

def electron_temperature(TaL,TaC,n,delta_V):
    """
    Given the peak antenna temperature in the line, TaL, the
    continuum temperature, TaC, the quantum number of the lower
    level, n, and the line FWHM, delta_V in km/s, returns the
    electron temperature in K.  Lockman and Brown, Ap.J. 207,
    436 (1976), eq. 2
    """
    # print TaL, TaC, n, delta_V
    return 10000*math.pow( (1.29/delta_V) \
                          *(TaC/TaL)      \
                          *pow(n/100.,-3.3), 0.87)

def velocity_width(Te,species,V_turb):
    """
    Given the electron temperature in K, the atomic species,
    and the turbulent velocity in km/s, returns the linewidth in m/s
    """
    v_squared = 2 * \
      Physics.k*Te/(Physics.m_u*Physics.AMU[species]) + \
       2*math.pow(V_turb,2)/3
    delta_v = 2*math.sqrt(math.log(2)*v_squared)/1000
    return delta_v

