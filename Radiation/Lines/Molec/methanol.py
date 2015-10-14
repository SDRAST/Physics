"""
Functions to support operations on methanol laboratory line data, like
selecting A or E state data and plotting A and E state levels and
transitions.
"""

#from Physics import *
import Physics as P
from Physics.Radiation.Lines.Molec import koln
from support.Graphics import colormap as C
#from pylab import *
import pylab as py
from math import log10

# converts MHz to 1/cm.
MHz2invCm = (P.pq(1e6,'Hz')/P.pq(1,'c')).inUnitsOf('1/cm').value

def get_A_state_transitions(line_data,v):
  """
  Extracts the A state transitions from the Koln catalog for the v=0 or
  v=1 vibrational states.  These are characterized by having the fourth quantum
  number as +0 or -0 for v=0 and +1 or -1 for v=1.
  The data are obtained like this::
    cat = open_molecule(32504)
    line_data = cat.readlines()
    A_transitions = get_A_state_transitions(line_data,0)
  """
  A_state_transitions = []
  for line in line_data:
    line_dict = koln.parse_catalog_line(line)
    if v==0:
      if ((line_dict['q upper'][3] == '+0') or (line_dict['q upper'][3] == '-0')) \
        and \
         ((line_dict['q lower'][3] == '+0') or (line_dict['q lower'][3] == '-0')):
        A_state_transitions.append(line_dict)
    elif v==1:
      if ((line_dict['q upper'][3] == '+1') or (line_dict['q upper'][3] == '-1')) \
        and \
         ((line_dict['q lower'][3] == '+1') or (line_dict['q lower'][3] == '-1')):
        A_state_transitions.append(line_dict)
  return A_state_transitions

def get_E_state_transitions(line_data,v):
  """
  Extracts the E state transitions from the Koln catalog for the vibrational
  state v=0 or v=1.  These are characterized by having the fourth quantum
  number as 0 or 1. The data are obtained like this::
    cat = open_molecule(32504)
    line_data = cat.readlines()
    E_transitions = get_E_state_transitions(line_data,0)
  """
  E_state_transitions = []
  for line in line_data:
    line_dict = koln.parse_catalog_line(line)
    if len(line_dict['q lower'][0]) > 2:
      # This takes care of cases where the quantum numbers run together
      line_dict['q lower'] = [line_dict['q lower'][0][0:-3],\
      line_dict['q lower'][0][-2:],line_dict['q lower'][1],line_dict['q lower'][2]]
    if v == 0:
      if (line_dict['q upper'][3] == '0') and (line_dict['q lower'][3] == '0'):
        E_state_transitions.append(line_dict)
    elif v == 1:
      if (line_dict['q upper'][3] == '1') and (line_dict['q lower'][3] == '1'):
        E_state_transitions.append(line_dict)
  return E_state_transitions

def plot_A_levels(A_transitions, QNindex, E_index, Emax, kmax):
  """
  Plots the energy levels of the A-states up to Emax (1/cm) and K=kmax.
  Use QNindex (e.g. 'q upper' or 'q lower') as the basis for extracting the
  quantum number data from the transition dictionaries.  Use E_index (typically
  'E lower') to extract the energy data. Because there are two A states, except
  for K=0, they are plotted side by side with different line styles.
  """
  levels = []
  for level in A_transitions:
    try:
      # See if the state was already found.
      levels.index(level[QNindex])
    except ValueError:
      # Add the state to the list of states
      levels.append(level[QNindex])
      j = int(level[QNindex][0])
      k = int(level[QNindex][1])
      state = level[QNindex][3]
      Elo = float(level["E lower"])
      if E_index == "E lower":
        E = Elo
      elif E_index == "E upper":
        E = Elo + float(level['freq'])*MHz2invCm
      if state == '+0':
        # A+ state
        if k == 0:
          # for J=0 plot full width black line
          color = 'k-'
          left  = k - 0.4
          right = k + 0.4
          if E < Emax and k <= kmax:
            py.text(k+0.4,E,level[QNindex][0])
        else:
          # for J>0, plot on left as blue line
          color = 'k-'
          left = k - 0.4
          right = k
      elif (state == '-0') and (k != 0):
        # A- state J>0, plot on right as green line
        color = 'k:'
        left  = k
        right = k + 0.4
        if E < Emax and k <= kmax:
          py.text(k+0.4,E,level[QNindex][0])
      else:
        # something is wrong
        color = 'r'
        left  = k - 0.4
        right = k + 0.4
      lw = 2
      py.plot([left,right],[E,E],color,lw=lw)
  # put state sign inside bottom of graph
  for K in range(1,kmax):
    py.text(K-0.25,1,'+')
    py.text(K+0.2,1,'-')

def plot_E_levels(E_transitions, QNindex, E_index, Emax, kmin, kmax):
  """
  Plots the energy levels of the A-states up to Emax (1/cm) and
  kmin <= K< = kmax. Use QNindex (e.g. 'q upper' or 'q lower') as the basis for
  extracting the quantum number data from the transition dictionaries.  Use
  E_index (typically 'E lower') to extract the energy data. Because there are
  two A states, except for K=0, they are plotted side by side with different
  line styles.
  """
  lower_E = []
  levels = []
  for level in E_transitions:
    try:
      levels.index(level[QNindex])
    except ValueError:
      levels.append(level[QNindex])
      j = int(level[QNindex][0])
      k_signed = level[QNindex][1]
      k = int(k_signed)
      Elo = float(level[E_index])
      if level[QNindex][3] == '0':
        # v=0 state
        left = k - 0.4
        color = 'k'
        right = k + 0.4
      else:
        # something is wrong
        color = 'r'
        left  = k - 0.4
        right = k + 0.4
      lw = 1
      py.plot([left,right],[Elo,Elo],color,lw=lw)
      if Elo < Emax and k <= kmax and k>= kmin:
        py.text(k+0.4,Elo,level[QNindex][0])

def plot_A_transitions(A_transitions,Amin,Amax,color=True,logA=True):
  """
  Plot the radiative decays in the list A_transitions. Amin and Amax are the
  rates for the bottom and top of the color bar. By default the logarithm of
  the rates are printed as colors.  If color is false, the rate or logarithm
  of the rate is plotted as a black line with the linewidth proportional to
  the rate.
  """
  for trans in A_transitions:
    Elo = float(trans['E lower'])
    f = float(trans['freq'])
    Eup = Elo + f*MHz2invCm
    A = koln.einstein_a(f,float(trans['str']),int(trans['g upper']))
    if logA == True:
      rate = py.log10(A)
    else:
      rate = A
    lineColor = C.strRgb(rate,Amin,Amax)
    # Plot separated lines for A+ and A- states
    # lower level
    if trans['q lower'][1] == '0':
      # for K=0 there is no distinction. Plot from/to level centers
      klo = int(trans['q lower'][1])
    elif trans['q lower'][3] == '+0':
      # A+ on the left
      klo = int(trans['q lower'][1]) - 0.2
    elif trans['q lower'][3] == '-0':
      # A- on the right
      klo = int(trans['q lower'][1]) + 0.2
    # upper level
    if trans['q upper'][1] == '0':
      # K=0, plot between line centers
      kup = int(trans['q upper'][1])
    elif trans['q upper'][3] == '+0':
      # A- on the left
      kup = int(trans['q upper'][1]) - 0.2
    elif trans['q upper'][3] == '-0':
      # A+ on the right
      kup = int(trans['q upper'][1]) + 0.2
    if color == True:
      if rate > Amin:
        # double width lines above the minimum value
        py.plot([kup,klo],[Eup,Elo],'-',lw=2,c=lineColor)
      else:
        # single width lines below the minimum value
        py.plot([kup,klo],[Eup,Elo],'-',lw=1,c=C.strRgb(Amin,Amin,Amax))
    else:
      if rate > Amin:
        linew = 3*(rate-Amin)/(Amax-Amin)
        py.plot([kup,klo],[Eup,Elo],'-',lw=linew, c='k')

def plot_E_transitions(E_transitions,Amin,Amax):
  """
  Plot the radiative decays in the list A_transitions. Amin and
  Amax are the rates for the bottom and top of the color bar.
  """
  for trans in E_transitions:
    Elo = float(trans['E lower'])
    f = float(trans['freq'])
    Eup = Elo + f*MHz2invCm
    A = koln.einstein_a(f,float(trans['str']),int(trans['g upper']))
    logA = log10(A)
    lineColor = C.strRgb(logA,Amin,Amax)
    klo = int(trans['q lower'][1])
    kup = int(trans['q upper'][1])
#    if abs(int(trans['q upper'][1])-int(trans['q lower'][1])) >2:
#      print trans['q upper'],'->', trans['q lower'],':', \
#        trans['freq']/1000.,', str=',trans['str']
    if logA > Amin:
      # double width lines above the minimum value
      py.plot([kup,klo],[Eup,Elo],'-',lw=2,c=lineColor)
    elif logA > Amin - 2:
      # single width lines below the minimum value
      py.plot([kup,klo],[Eup,Elo],'-',lw=1,c=C.strRgb(logA,Amin,Amax))
