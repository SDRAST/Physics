"""
Functions related to calculating the rotational energy of asymmetric
molecules.  Townes and Schawlow, Ch. 4
"""
from pylab import poly1d

def asymmetry (A,B,C):
  """
  Ray's asymmetry parameter for molecular rotation.
  For a prolate symmetric top (B = C), kappa = -1.
  For an oblate symmetric top (B = A), kappa = +1.
  See Townes and Schawlow, Ch. 4.
  """
  return (2.*B - A - C)/(A - C)

def b_prolate(kappa):
  """
  0 = prolate <= b_P <= -1 = oblate
  Townes and Schawlow, Ch. 4
  """
  return (kappa+1.)/(kappa-3.)

def b_oblate(kappa):
  """
  -1 = oblate <= b_O <= 0 = prolate
  Townes and Schawlow, Ch. 4
  """
  return (kappa-1.)/(kappa+3.)

def asym_quantum_factor(J,b):
  """
  This takes the places of K^2 in calculating the energy levels
  for asymmetric rotators. Townes and Schawlow, Ch. 4.  For
  J > 6 this returns an empty tuple.  Note that it doesn't matter which version of
  b is used since b_prolate(kappa) = b_oblate(-kappa) and the equations are
  symmetric in b or depend on b**2.
  """
  roots = ()
  if J == 0:
    roots = (0,)
  elif J == 1:
    roots =  (0., 1+b, 1-b)
  elif J == 2:
    roots = ( 4., 1-3*b, 1+3*b)
    p = poly1d([1, -4, -12*b**2])
    roots = roots + tuple(p.r)
  elif J == 3:
    roots = (4.,)
    p = poly1d([1, -4, -60*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -10+6*b, 9-54*b-15*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -10-6*b, 9+54*b-15*b**2])
    roots = roots + tuple(p.r)
  elif J == 4:
    p = poly1d([1, -10*(1-b), 9-90*b-63*b**2])
    roots = tuple(p.r)
    p = poly1d([1, -10*(1+b), 9+90*b-63*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -20, 64-28*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -20, 64-208*b**2, 2880*b**2])
    roots = roots + tuple(p.r)
  elif J == 5:
    p = poly1d([1, -20, 64-108*b**2])
    roots = tuple(p.r)
    p = poly1d([1, -20, 64-528*b**2,6720*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -35+15*b, 259-510*b-213*b**2, -225+3375*b+4245*b**2-675*b**3])
    roots = roots + tuple(p.r)
    p = poly1d([1, -35-15*b, 259+510*b-213*b**2, -225-3375*b+4245*b**2+675*b**3])
    roots = roots + tuple(p.r)
  elif J == 6:
    p = poly1d([1, -35+21*b, 259-714*b-525*b**2, -225+4725*b+9165*b**2-3465*b**3])
    roots = tuple(p.r)
    p = poly1d([1, -35-21*b, 259+714*b-525*b**2, -225-4725*b+9165*b**2+3465*b**3])
    roots = roots + tuple(p.r)
    p = poly1d([1, -56, 784-336*b**2, -2304+9984*b**2])
    roots = roots + tuple(p.r)
    p = poly1d([1, -56, 784-1176*b**2, -2304+53664*b**2, -483840*b**2+55440*b**4])
    roots = roots + tuple(p.r)
  else:
    roots = ()
  return roots

def Ejk(A,B,C,J):
  """
  Rotational energy of an asymmetric molecule using Eq. 4-4 and 4-5 in
  Townes and Schawlow.  Returns energies in units used for A,B and C
  """
  kappa = asymmetry(A,B,C)
  # assume prolate form
  b_P = b_prolate(kappa)
  ws = asym_quantum_factor(J,b_P)
  # print "w's=",ws
  result = []
  for w in ws:
    result.append( (B+C)*J*(J+1)/2. + (A - (B+C)/2. )*w )
  return result

def plot_E_vs_kappa(Amax,C,maxJ):
  """
  Plots diagram showing how energy of an asymmetric rotor depends on its
  asymmetry as it varies between prolate and oblate, and how the J(-K,+K)
  labelling arises.  Townes and Schawlow , Ch. 4.
  This assumes that the volume of the moment of inertia ellipsoid is a
  constant::
   Ia*Ib*Ic = 3*V/(2*pi)
  or::
   (h/8 pi**2)**3(ABC) = 3*V/(2*pi)
  or::
   A*B*C = K, a constant
  The ellipsoid's minimum semi-axis C is also a constant. So in the prolate case,
  B=C and K = Amax*C*C. In the oblate case, Amin=B and K = A*A*C.
  The constraints are then::
   2*B = (1+kappa)*A + (1-kappa)*C
  and::
   A = Amax*C/B
  """
  n_kappas = 21
  kappas = linspace(-1,1,n_kappas)
  for J in range(maxJ):
    n_columns = 2*J+1
    # create a matrix on n_kappas rows and n_
    E = zeros((n_kappas,n_columns),float)
    for i in range(n_kappas):
      kappa = kappas[i]
      p = poly1d([2,(kappa-1)*C,-(kappa+1)*Amax*C])
      if p.r[0] > 0.:
        B = p.r[0]
      else:
        B = p.r[1]
      print B
      A = Amax*C/B
      # This should yield n_columns of energy values for this kappa
      Es = Ejk(A,B,C,J)
      E[i,:] = Es
    # Now we have a 2D array of energies to plot
    for k in range(n_columns):
      # select a line style and plot
      if J%3 == 0:
        ls = "-"
      elif J%3 == 1:
        ls = "--"
      elif J%3 == 2:
        ls = ":"
      else:
        ls = "-."
      plot(kappas,E[:,k],label=r"$J_{\tau}="+str(J)+"_{"+str(k-J)+"}$",ls=ls,lw=2)
    # label the lines
    for K in range(J+1):
      # For prolate, B=C
      E_prolate = C*J*(J+1)+(Amax-C)*K**2
      # For oblate, B=A, using the last value of A
      E_oblate =  A*J*(J+1)+(C-A)*K**2
      text(-0.98+0.07*(J%2),E_prolate,r"$"+str(J)+"_{"+str(K)+"}$")
      text( 0.93-0.07*(J%2),E_oblate,r"$"+str(J)+"_{"+str(K)+"}$")

def test(A,B,C,maxJ):
  """
  Checks the calculation of energy levels.  For example, for water::
   >>> test(835.83910,435.347353,278.139826,5)
   0 [0.0]
   1 [713.48717899999997, 1113.978926, 1271.186453]
   2 [4056.8435790000003, 2855.3683380000002, 2383.7457570000001,
      4094.7813912145548, 2102.5237247854457]
   3 [6197.3051160000005, 6374.3863667070382, 4103.8418232929625,
      8620.1335561107444, 5204.2902778892558, 8614.2071531433267, 4266.9715188566752]
   4 [11569.54077149916, 8277.1955485008384, 11529.522215260266, 6745.1388347397333,
      14830.335414585266, 9021.318375414734, 14831.098979756451, 9490.7431163792135,
      6664.6834838643363]
  Divide by 29.997 to convert GHz to 1/cm
  """
  for J in range(maxJ):
    print J,Ejk(A,B,C,J)

if __name__ == "__main__":
  C=25
  A=125
  maxK = 4

  plot_E_vs_kappa(A,C,maxK)
  a = axis()
  b = [a[0], a[1], -10, a[3]]
  axis(b)

  rcParams.update({'legend.fontsize': 10})
  legend(loc=9)
  title(r"$\mathrm{Asymmetric~rotor,~A}_{\mathrm{max}}=$"+str(A)+"$,~\mathrm{C}=$"+str(C))
  xlabel(r"$\mathrm{Asymmetry}$")
  ylabel(r'$E(J,\tau)/h~(GHz)$')
  show()