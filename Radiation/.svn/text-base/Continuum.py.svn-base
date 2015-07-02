"""Functions relation to emission mechanisms
   by Tom Kuiper 2008 Nov 10"""

import math as M
import Physics as P

def free_free_absorp_coef(n_e,n_i,T,f):
  """Absorption coefficient for free-free radiation.
  n_e is the electron density,
  n_i the ion density,
  T is the electron temperature in K, and
  f the frequency in GHz.
  The result in 1/pc."""
  result = \
      0.08235*n_e*n_i/M.pow(f,2.1)/M.pow(T,1.35)
  return result

def free_free_absorp_coefPQ(n_e,n_i,T,f):
  """Returns a physical quantity for the free-free absorption coefficient
  given the electron density, ion density, kinetic temperature and frequency
  as physical quantities. From Shklovsky (1960) as quoted by Kraus (1966)."""
  value = 9.8e-13 * n_e.inBaseUnits().value * n_i.inBaseUnits().value \
          * M.pow(T.inBaseUnits().value,-1.5) * M.pow(f.inBaseUnits().value,-2) \
          * (19.8 + M.log(M.pow(T.inBaseUnits().value,1.5)/f.inBaseUnits().value))
  return P.pq(value,'1/m')
  
def BB_intensity(T,f):
    """Blackbody radiation intensity, W/m^2/Hz/rad^2,
    f is in Hz, and T is in K."""
    left_term = 2*P.h*M.pow(f,3)/M.pow(P.c,2)
    exponent = P.h*f/P.k/T
    if exponent < 700:
      right_term = M.pow(M.e,exponent) - 1
      result = left_term/right_term
    else:
      result = 1e-304
    return result

def BB_intensityPQ(T,f):
    """Physical quantity object for Blackbody radiation intensity,
    in W/m^2/Hz/rad^2, given physical quantity objects for f and T."""
    left_term = 2*P.pq(1,'hplanck')*f**3/P.pq(1,'c')**2
    print "2 h f^3 / c^2 =", left_term.inBaseUnits()
    print "h =", P.pq(1,'hplanck').inBaseUnits()
    print "h f =",(P.pq(1,'hplanck')*f).inBaseUnits()
    print "k =",P.pq(1,'k').inBaseUnits()
    print "k T =", (P.pq(1,'k')*T).inBaseUnits()
    exponent = (P.pq(1,'hplanck')*f)/(P.pq(1,'k')*T)
    print "hf/kT =",exponent
    if exponent < 700:
      right_term = M.pow(M.e,exponent) - 1
      result = left_term/right_term
    else:
      result = 1e-304
    return result

