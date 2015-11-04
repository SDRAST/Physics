"""
Provides class, data and methods for computing with polarized signals
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
    if self.V:
      self.ellipse['eccen'] = sqrt(self.Q**2 + self.U*2)/self.V
    else:
      self.ellipse['eccen'] = float('inf')

# --------------------------------------- DATA --------------------------------

quad_hybrid = array([[1+0j, 0+1j],
                     [1+0j, 0-1j]])/sqrt(2)
