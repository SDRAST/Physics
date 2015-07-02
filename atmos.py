"""
Physics.atmos - functions for atmospheric opacity, etc.
"""

atm = 10.1325            #  newton/cm^2 = 1013.25 mb
newton = 100000          #  dyn = g cm / s^2
amu = 1.66E-24           #  atomic mass unit in g
k = 1.38046E-16          #  erg/K                 Boltzmann constant

def rel_humidity(temp, dewpoint):
  """
  Relative humidity

  input air temperature and dew point temperature in deg C
  returns relative humidity in percent
  """
  rel_humidity = 100*critical_humidity(dewpoint)/critical_humidity(temp)
  return rel_humidity

def vapor_pressure(temp):
  """ 
  Vapor pressure

  input temp in deg C; returns vapor pressure of water in mm Hg.
  Formula from quadratic fit of log10(vapor pressure) vs.
  1/(temp in Kelvin) for data in the range 0 < deg C < 220
  taken from Table 19-5 in
  Sears and Zemansky, COLLEGE PHYSICS, Addison-Wesley (1957)
  and adjusted for the range -20 < deg C < 0 to give the
  values of critical humidity listed for -10 and -20 C in
  METEOROLOGY FOR AIR NAVIGATORS, Can. Dept. of Transport,
  Form No. 2261-2 (1943)
  """
  theta=1/(273.15 + temp)
  if temp > 0:
    vapor_pressure= 10**(7.86932-theta*1567.76*(1+69.906*theta))
  else:
    vapor_pressure= 10**(7.86932-theta*1350.0*(1+125.0*theta))
  return vapor_pressure

def critical_humidity(temp):
  """
  Water vapor density of a fully saturated gas

  input temp in deg C
  returns humidity of fully saturated gas in gm/m^3
  """
  vap_pres_mm = vapor_pressure(temp)
  vap_pres_atm  = vap_pres_mm/760.
  newton_per_sq_cm = vap_pres_atm * atm
  dyn_per_sq_cm = newton_per_sq_cm * newton
  kelvin = 273.15 + temp
  number_per_cc = number_density(dyn_per_sq_cm, kelvin)
  number_per_m3 = number_per_cc*1000000.0
  critical_humidity = number_per_m3 * 18 * amu
  return critical_humidity

def number_density(pres, temp):
  """
  Number density of gas molecules

  input pressure in dyne/cm^2, temp in deg K
  returns number of molecules per cm^3
  """
  number_density = pres / (k * temp)
  return number_density

def abs_humidity(rel_humidity, temp):
  """
  Water vapor density

  input percent relative humidity and temp in deg C
  returns absolute humidity in g/m^3
  """
  abs_humidity = rel_humidity * critical_humidity(temp)/100
  return abs_humidity

