# -*- coding: latin-1 -*-
"""
Provides class, data and methods for computing with polarized signals

The IAU defines polarization in terms of a right-handed coordinate system with
the wave propagating towards the observer as the Z axis. The other coordinates
are Right Ascension (X) and Declination (Y). The position angle of the electric
vector maximum (for linear or elliptical polarization) is measured from North
through East.

Note that the polarization angle theta_P is related to the usual polar angle
[theta = cotan(Y/X)] by theta_P = pi/2 - theta, or theta_P = cotan(X/Y).

A left circularly polarized wave is one in which the electric vector rotates 
counterclockwise when viewed along the axis of propagation, that is, increasing
polar angle. From the point of view of an observer, the vector rotates
clockwise. To drive the point home, an incoming left circularly polarized wave
rotates to the right.

An LCP wave can be obtained from two linearly polarized waves in which the wave
parallel to Y (X=0 plane) is pi/2 ahead of the wave which is parallel to X
(Y=0 plane)::
  E    = E  sin(2 pi f t) + E  cos(2 pi f t)
   LCP    H                  V
   
       = E  sin(2 pi f t) + E  sin(2 pi x/lambda + pi/2)
          H                       V

Expressing that using complex notation (which accounts for the phase as well as
the amplitude of each component::
             
  E    = (E  + i E )/sqrt(2)
   LCP   ( H      V)
   
  E    = (i E  + E )/sqrt(2)
   RCP   (   H    V)

In matrix notation::
  |E   |   | 1   i | |E |
  | LCP|   |       | | H|
  |    | = |       | |  |
  |E   |   | i   1 | |E |
  | RCP|   |       | | V|
"""
from numpy import abs, array
from numpy.linalg import inv
from math import atan, sqrt
import logging

module_logger = logging.getLogger(__name__)

class PolarizationError(Exception):
  """
  Class for reporting exceptions in module Polarization
  """
  def __init__(self,value,message):
    """
    Initialize a class instance
    """
    ValueError.__init__(self,value,message)
  
  
class Signal(object):
  """
  Polarized signal class
  
  Currently does not include a non-polarized signal component
  
  Attributes::
    circular - orthogonal circular components of the signal
    ellipse  - dict with properties of the polarization ellipse
    I        - Stokes I
    linear   - orthogonal linear components of the signal
    logger   - logging.Logger
    Q        - Stokes Q
    Stokes   - (I,Q,U,V)
    U        - Stokes U
    V        - Stokes V
  """
  def __init__(self, mode="linear", E=(1+0j, 0+0j)):
    """
    Initialize a polarized signal
    
    @param mode : circular or linear
    @type  mode : str
    
    @param E : initialize a signal, default linear along X
    @type  E : tuple of complex (X,Y)
    """
    object.__init__(self)
    self.logger = logging.getLogger(module_logger.name+".Signal")
    if mode[0].lower() == 'c':
      self.from_circular(E)
    elif mode[0].lower() == 'l':
      self.from_linear(E)
    else:
      self.logger.error("__init__: invalid mode %s", mode)
      raise PolarizationError()
    self.polarization_ellipse()
    
  def from_linear(self, E=(1+0j, 0+0j)):
    """
    Create a signal with linear orthogonal components
    
    @param E : initialize a signal, default linear along X
    @type  E : tuple of complex (X,Y)
    """
    self.linear = array(E).reshape((2,1))
    self.circular = self.to_circular()
    self.Stokes_from_linear()
    return self
    
  def from_circular(self, E=(1+0j, 0+0j)):
    """
    Create a signal with circular orthogonal components
    
    @param E : initialize a signal, default left circular
    @type  E : tuple of complex (L,R)
    """
    self.circular = array(E).reshape((2,1))
    self.linear = self.to_linear()
    self.Stokes_from_circular()
    return self
    
  def to_circular(self):
    """
    Converts signal with orthogonal linear components to circular components
    """
    return quad_hybrid.dot(self.linear)
  
  def to_linear(self):
    """
    Converts signals with orthogonal circular components to linear components
    """
    return inv(quad_hybrid).dot(self.circular)
  
  def Stokes_from_linear(self):
    """
    Computes Stokes parameters from orthogonal linear components
    """
    (Ex,Ey) = self.linear
    self.logger.debug("Stokes_from_linear: (Ex, Ey)  = %s", (Ex,Ey))
    (Exc,Eyc) = self.linear.conj()
    self.logger.debug("Stokes_from_linear: (Ex*,Ey*) = %s", (Exc,Eyc))
    (Sxx,Syy) = abs(self.linear*self.linear.conj())
    self.logger.debug("Stokes_from_linear: Sxx, Syy  = %s", (Sxx,Syy))
    Sxy = Ex*Eyc
    Syx = Ey*Exc
    self.logger.debug("Stokes_from_linear: Sxy, Syx  = %s", (Sxy,Syx))
    self.I = float(Sxx+Syy)
    self.Q = float(Sxx-Syy)
    self.U = float((Sxy+Syx).real)
    self.V = float(((0-1j)*(Sxy-Syx)).real)
    self.Stokes = self.I,self.Q,self.U,self.V
    return self.Stokes
    
  def Stokes_from_circular(self):
    """
    Computes Stokes parameters from orthogonal circular components
    """
    (El,Er) = self.circular
    self.logger.debug("Stokes_from_circular: (El, Er)  = %s", (El,Er))
    (Elc,Erc) = self.circular.conj()
    self.logger.debug("Stokes_from_circular: (El*,Er*) = %s", (Elc,Erc))
    (Sll,Srr) = abs(self.circular*self.circular.conj())
    self.logger.debug("Stokes_from_circular: Sll, Srr  = %s", (Sll,Srr))
    Slr = El*Erc
    Srl = Er*Elc
    self.logger.debug("Stokes_from_circular: Slr, Srl  = %s", (Slr, Srl))
    self.I = float(Sll+Srr)
    self.Q = float((Slr+Srl).real)
    self.U = float(((0-1j)*(Slr-Srl)).real)
    self.V = float(Sll-Srr)
    self.Stokes = self.I,self.Q,self.U,self.V
    self.logger.debug("Stokes_from_circular: Stokes: %s", self.Stokes)
    return self.Stokes

  def polarization_ellipse(self):
    """
    Computes parameters of the polarization ellipse.
    """
    self.ellipse = {}
    self.ellipse['d_lin'] = sqrt(self.Q**2 + self.U**2)/self.I
    self.ellipse['d_cir'] = abs(self.V)/self.I
    self.ellipse['d'] = sqrt(self.Q**2 + self.U**2 + self.V**2)/self.I
    if self.Q:
      self.ellipse['angle'] = 0.5*atan(self.U/self.Q)
    else:
      self.ellipse['angle'] = float('NaN')
    if (self.Q**2 + self.U**2):
      self.ellipse['eccen'] = self.V/sqrt(self.Q**2 + self.U**2)
    else:
      self.ellipse['eccen'] = float('inf')

# --------------------------------------- DATA --------------------------------

# This converts from linear to circular modes.  Its inverse converts from
# circular to linear modes.
quad_hybrid = array([[1+0j, 0+1j],
                     [0+1j, 1+0j]])/sqrt(2)

if __name__ == "__main__":
  """
  Test cases from 
  https://en.wikipedia.org/wiki/Stokes_parameters#Relation_to_the_polarization_ellipse
  |I|   | 1|
  |Q|   | 1| 
  |U| = | 0| 	Linearly polarized (horizontal)
  |V|   | 0|

      | 1| 
      |-1| 
      | 0| 	Linearly polarized (vertical)
      | 0|
      
      | 1|
      | 0| 
      | 1|  Linearly polarized (+45°)
      | 0|

      | 1| 
      | 0|
      |-1| 	Linearly polarized (−45°)
      | 0|
      
      | 1| 
      | 0| 
      | 0|  Right-hand circularly polarized
      | 1|
      
      | 1| 
      | 0| 
      | 0| 	Left-hand circularly polarized
      |-1|
      
      | 1|
      | 0| 
      | 0| 	Unpolarized
      | 0|
  The tests below agree with Kraus, Radio Astronomy, section 4-4 (1966). The 
  confusion is clarified in 
  https://en.wikipedia.org/wiki/Circular_polarization#Left.2Fright_handedness_conventions
  """
  X = Signal()
  print "\nlinear horizontal"
  print "X,Y electric field =", X.linear.flatten()
  print "Stokes =", X.Stokes

  Y = Signal(E=(0,1))
  print "\nlinear vertical"
  print "X,Y electric field =", Y.linear.flatten()
  print "Stokes =", Y.Stokes
  
  XYp45 = Signal(E=(sqrt(2)/2,sqrt(2)/2))
  print "\nlinear +45 deg"
  print "X,Y electric field =", XYp45.linear.flatten()
  print "Stokes =", XYp45.Stokes
  
  XYm45 = Signal(E=(-sqrt(2)/2,sqrt(2)/2))
  print "\nlinear -45 deg"
  print "X,Y electric field =", XYm45.linear.flatten()
  print "Stokes =", XYm45.Stokes
  
  LC = Signal(mode='circular')
  print "\nleft circular"
  print "L,R electric fields =", LC.circular.flatten()
  print "Stokes =", LC.Stokes
  
  RC = Signal(mode='circular',E=(0,1))
  print "\nright circular"
  print "L,R electric fields =", RC.circular.flatten()
  print "Stokes =", RC.Stokes
  
