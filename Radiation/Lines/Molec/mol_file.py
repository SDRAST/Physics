"""
Functions which use the atomic position data in MOL files such as those
obtained from
http://webbook.nist.gov/chemistry/name-ser.html
or
http://www.ebi.ac.uk/chebi/
"""

from re import search
import Physics as P
from numpy import zeros, linalg

def get_position_data(mol_data):
  """
  This obtains data from a MOL file such as may be obtained from the
  Chemical Entities of Biological Interest (ChEBI) website, or
  the NIST website.
  """
  global n_atoms
  x = []
  y = []
  z = []
  ID = []
  index = 0
  for line in mol_data:
    if search('V2000',line):
      metadata = line.strip().split()
      n_atoms = int(metadata[0])
      n_bonds = int(metadata[1])
      for i in range(n_atoms):
        # print mol_data[index+i+1].strip()
        atom_data = mol_data[index+i+1].strip().split()
        x.append( atom_data[0] )
        y.append( atom_data[1] )
        z.append( atom_data[2] )
        ID.append( atom_data[3] )
      break
    index += 1
  return ID,x,y,z

def format_LaTeX_table(ID,x,y,z):
  """
  Makes a LaTeX formatted table of the atomic position data in a MOL file.
  """
  print "\\begin{tabular}{c|ccc}"
  print "{\\bfseries Atom} & {\\bfseries X} & {\\bfseries Y} & {\\bfseries Z} \\\\"
  print "{\\bfseries ID} & $\stackrel{\\circ}{\\mathrm{A}}$ & $\stackrel{\\circ}{\\mathrm{A}}$ & $\stackrel{\\circ}{\\mathrm{A}}$ \\\\"
  for i in range(len(ID)):
    print ID[i]," & ", x[i], " & ", y[i]," & ", z[i]," \\\\"
  print "\\end{tabular}"

def inertia_matrix(ID,x,y,z,diag=False):
  """
  Computes the rotation constants from MOL file data.  Returns a tuple
  with the inertia tensor, whether it is diagonalized or not, and the
  atom coordinates relative to the center of mass. For the coordinates the
  first index is the axis, the second the atom. 
  Diagnostics print out if diag = 1.
  """
  # Create arrays with masses and coordinates for all atoms
  n_atoms = len(ID)
  m = zeros((n_atoms),float)       # masses of nuclei
  r = zeros((3,n_atoms),float)     # x,y,z positions of nuclei
  if diag == True:
    print "Position relative to origin"
    print " m      x        y       z"
  for i in range(n_atoms):
    m[i] = P.AMU[ID[i]]
    r[0,i] = x[i]
    r[1,i] = y[i]
    r[2,i] = z[i]
  if diag == True:
    for nucleus in range(n_atoms):
      print "%5.2f  %6.3f  %6.3f  %6.3f" % \
          (m[nucleus], r[0,nucleus],r[1,nucleus],r[2,nucleus])
  # ------------- Find center of mass
  r_cm = zeros((3))
  m_tot = 0
  for axis in range(0,3):
    r_cm[axis] = 0.0
  for nucleus in range(0,n_atoms):
    for axis in range(0,3):
        r_cm[axis] += r[axis,nucleus]*m[nucleus]
    m_tot += m[nucleus]
  for axis in range(0,3):
    r_cm[axis] = r_cm[axis]/m_tot
  if diag == True:
    print "Position of center of mass:"
    print "X=%6.3f, Y=%6.3f, Z=%6.3f" % (r_cm[0],r_cm[1],r_cm[2])
  # ------------- change coordinates to center of mass system
  dist = zeros((n_atoms),float)
  for nucleus in range(0,n_atoms):
    dist[nucleus] = 0
    for axis in range(0,3):
        r[axis,nucleus] -= r_cm[axis]
        dist[nucleus] += pow(r[axis,nucleus],2)
    dist[nucleus] = sqrt(dist[nucleus])
  if diag == True:
    print "Positions relative to center of mass"
    print "  m      x       y       z       r"
    for nucleus in range(0,n_atoms):
      print "%5.2f  %6.3f  %6.3f  %6.3f  %6.3f" % \
	       (m[nucleus],r[0,nucleus],r[1,nucleus],r[2,nucleus],dist[nucleus])
  # ------------------------------------------ compute the inertia tensor
  # see Kittel et al., eqn. 14, p.240
  inertia= zeros((3,3),float)
  for i in range(0,3):
    for j in range(0,3):
      inertia[i,j] = 0
      for nucleus in range(0,n_atoms):
        if i==j:
          inertia[i,j] += \
            m[nucleus]*(pow(dist[nucleus],2) - pow(r[i,nucleus],2))
        else:
          inertia[i,j] -= m[nucleus]*r[i,nucleus]*r[j,nucleus]
  if diag == True:
    print "Rotational inertia tensor (in AMU*angstroms^2):"
  # Check if the tensor is diagonalized
  diagonalized = True
  for j in range(0,3):
    for i in range(0,3):
      if diag == True:
        print  "%6.3f" % inertia[i,j],
      if i<>j and abs(inertia[i,j])>1.0E-15:
        diagonalized = False
    if diag == 1:
      print
  return inertia, diagonalized, r

def diagonalize_matrix(inertia, diag=False):
  """
  Diagonalize the inertia tensor.  Print out diagnostics if diag = 1.
  """
  diagonal=zeros((3),float)
  off_diagonal=zeros((3),float)

  diagonal,eig_vectors = linalg.eig(inertia)
  if diag == True:
    print "Diagonal:"
    print diagonal
    print "Eigen vectors:"
    print eig_vectors
  return diagonal,eig_vectors

def rotation_constants(diagonal,diag=False):
  if diag == True:
    print "Eigenvalues:        Rotational constants:"
    print "   (AMU*A^2  (kg m^2)   (GHz)"
  rc=zeros((3))
  result = []
  for i in range(0,3):
    if diag == True:
      print chr(65+i),               # axis: A, B, C
      print "%8.5f" % (diagonal[i]), # diagonal matrix element
    # moment of inertia in MKS
    temp = diagonal[i]*P.m_u * pow(P.pq(1,'Ang').inBaseUnits().value,2)
    if diag:
      print "%9.2e" % temp,

    if diagonal[i] <> 0 :
      rot_const = P.h_bar/temp/4/pi/1e9
      result.append(rot_const)
    if diag == True:
      print "%9.4f" % rot_const,
      print
  return result

if __name__ == "__main__":
  # get the data
  mol_fd = open('methanol.mol','r')
  mol_data = mol_fd.readlines()
  mol_fd.close()
  ID,x,y,z = get_position_data(mol_data)

  # compute the inertia matrix
  inertia,diagonalized = inertia_matrix(ID,x,y,z,True)
  # diagonalize the inertia matrix
  if diagonalized == False:
    diagonal = diagonalize_matrix(inertia,True)
  # compute the rotation constants
  print rotation_constants(diagonal)
  