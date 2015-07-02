"""
module jpl - Routines to access the JPL Line Catalog
"""

import urllib2
from Physics import *
import Math as M
import operator
import text
import re
from os.path import exists

diag = False

jpl_url = 'http://spec.jpl.nasa.gov/ftp/pub/catalog/'
local_path = '/usr/local/line_cat/jpl/'

# This keeps track of where to find the molecular data
version = {}

def catalog_versions(code):
  global version
  # First find the available versions
  fd = urllib2.urlopen(jpl_url)
  lines = fd.readlines()
  fd.close()
  versions = []
  for line in lines:
    txt = text.remove_html_tags(line.strip())
    if txt != '' and txt[0] == 'v':
      parsed = txt.split('/')
      versions.append(parsed[0])
  # Now check each version's subdirectory, last first
  versions.sort(reverse=True) # should not be needed
  print versions
  return versions

def latest_version(code):
  code_string = ('%06d' % code)
  filename = 'c'+code_string+'.cat'
  for v in catalog_versions(code):
    print "Trying",jpl_url+v+'/'+filename
    try:
      mol_fd = urllib2.urlopen(jpl_url+v+'/'+filename)
    except urllib2.HTTPError:
      continue
    if not version.has_key(code):
      version[code] = v
    mol_fd.close()
    return v
  # Now try original catalog
  print "Trying",jpl_url+filename
  try:
    mol_fd = urllib2.urlopen(jpl_url+filename)
  except urllib2.HTTPError:
    return None
  if not version.has_key(code):
    version[code] = v
  mol_fd.close()
  return v
  
  
def open_molecule(code):
  """
  Returns a file descriptor for the data file for the molecule
  with the specified code.

  Catalog Version
  ===============
  It looks in the latest catalog version, then the one before, etc.
  until the molecule is found.

  @type code : int
  @param code : molecule identifier
  """
  code_string = ('%06d' % code)
  filename = 'c'+code_string+'.cat'
  v = latest_version(code)
  if v:
    catpath = jpl_url+v+'/'+filename
  else:
    catpath = jpl_url+filename
  print "Opening", catpath
  try:
    mol_fd = urllib2.urlopen(catpath)
    return mol_fd
  except urllib2.HTTPError:
    return None

def parse_catalog_line(line):
  """
  Parses a line of data extracted from a molecular data file.
  
  For the encoding details see
  http://spec.jpl.nasa.gov/ftp/pub/catalog/README

  @type line : str
  @param line : a line of data from a JPL spec data file

  @return: dictionary with fairly obvious keys::
    'freq'     - frequency in MHz
    'freq err' - frequency formal error in MHz
    'int'      - log10(integrated intensity in units of nm^2 MHz) at 300K
    'deg free' - Degrees of freedom in the rotational partition function
                 0 for atoms,
                 2 for linear molecules, and
                 3 for nonlinear molecules.
    'E lower'  - energy of the lower level in cm^{-1}
    'g upper'  - degeneracy of the upper state
    'tag'      - molecule code
    'qn format'- see documentation for quantum_label()
  """
  data = {}
  data['freq'] = line[0:13].strip()
  data['freq err'] = line[13:21].strip()
  data['int'] = line[21:29].strip()
  data['deg free'] = line[29:31].strip()
  data['E lower'] = line[31:41].strip()
  data['g upper'] = line[41:44].strip()
  data['tag'] = line[45:51]
  data['qn format'] = line[51:55]
  q_low = []
  q_up = []
  for i in range(6):
    ptr_up = 55+2*i
    ptr_low = 55+12+2*i
    q_up.append(line[ptr_up:ptr_up+2].strip())
    q_low.append(line[ptr_low:ptr_low+2].strip())
  data['q upper'] = q_up
  data['q lower'] = q_low
  return(data)

def extract_range(moltag,f_low,f_up):
  """
  Extract all the lines in the specified frequency range

  @type moltag : int
  @param moltag : molecule ID

  @type f_low : float
  @param f_low : lower end of the frequency range in MHz
  
  @type f_up : float
  @param f_up : upper end of the frequency range in MHz

  @return: list of dictionaries with parsed line data
  """
  mol_fd = find_latest_version(moltag)
  mol_data = mol_fd.readlines()
  mol_fd.close()
  result = []
  for line in mol_data:
    line_data = parse_catalog_line(line)
    print line_data
    if ((float(line_data['freq']) > float(f_low)) & \
        (float(line_data['freq']) < float(f_up) ) & \
        (float(line_data['int'])  > -4          ) ):
            result.append(line_data)
  return(result)

def quantum_label(tag,q_nums,q_fmt,deg_free):
  """
  Decode the quantum numbers in the JPL catalog.

  References
  ==========
  Refer to section V of the JPL Submm Catalog Supplement
  and Townes and Schawlow, Dover edition 1975, p.176.
  This version follows the JPL convention in use on 2008 Jun 15

  Quantum code
  ============
  The quantum code is a four-digit number.

  Molecule Type
  -------------
  The first two digits (most significant digits) give the molecule type::
    0 - atom (deg_free = 0)
    1 - linear, sigma (deg_free = 2)
    2 - linear case b (deg_free = 2) or symmetric rotor (deg_free = 3)
    3 - linear case a (2S+1 odd; deg_free = 2) or
        asymmetric rotor (deg_free = 3)
    8 - linear case a (2S+1 even; deg_free = 2)
   13 - symmetric rotor with vibration (deg_free = 3)
   14 - asymmetric rotor with vibration (deg_free = 3)
   73 - ammonia with hyperfine splitting

  Number of Primary Quantum Numbers
  ---------------------------------
  The molecule type modulo 5 gives the number of primary quantum numbers

  Number of Quantum Numbers
  -------------------------
  This is given by the last digit of the quantum code. If the number of
  quantum numbers  > number of primary quantum numbers, then the degeneracy
  is derived from the last quantum number.

  Half Integer Bits
  -----------------
  The third digit of the quantum code is a 3-bit binary number that
  indicates whether the last three quantum numbers have half-integer
  values.  If so, they have been rounded up, e.g., 2 => 1.5

  @type tag : int
  @param tag : molecule ID

  @type q_nums : list
  @param q_nums : up to six quantum numbers

  @type q_fmt : int
  @param q_fmt : quantum number format code

  @return: tuple::
  1) a list of quantum number labels such as
     ['J', 'K_a', 'K_c', 'v_t', '', ''] and
  2) a TeX-able string for the quantum assignment of the state with
     the quantum numbers 'q_nums'.
  Some of the labels may depend on the values of some of the
  quantum numbers.
  
  """
  if diag:
    print q_nums
  quantum_code = str(q_fmt)
  mol_type = int(quantum_code[0:2])
  if diag:
    print "molecule type =",mol_type
  # mol_type MOD 5 is the number of primary quantum numbers
  num_primary_qn = mol_type % 5
  if diag:
    print "number of primary quantum numbers =",num_primary_qn
  # half integer code
  half_int_code = int(quantum_code[2])
  if diag:
    print "half integer code =",half_int_code
  # number of quantum numbers for each state.
  num_qn = int(quantum_code[3])
  
  qn_id = ["","","","","",""]
  qn_str = ["","","","","",""]
  if mol_type == 0:
    # atoms
    if diag:
      print "Processing atom"
    #    test half integer code
    # defaults
    qn_id[0] = "J"; qn_id[1] = "F"
    if half_int_code & 1:
      qn_id[0] = "J+0.5"
      qn_str[0] = str(2*q_nums[0]-1)+"/2"
    #    test first quantum number
    if q_nums[0] == 0:
      qn_id[1] = ""
    elif q_nums[0] == 1:
      qn_id[1] = "F+.5"
      qn_str[1] = "(F="+str(2*q_nums[1]-1)+"/2)"
    elif q_num[0] == 2:
      qn_id[1] = "F"
      qn_str[1] = "(F="+str(q_nums[1])+")"
  elif mol_type == 1 or mol_type == 2 or mol_type == 3 or mol_type == 8:
    if diag:
      print "Processing linear or diatomic molecule"
    # linear and diatomic molecules
    #    test molecular type code
    if mol_type == 1:
      # sigma state
      # defaults
      qn_id[0] = "N"
      qn_id[1] = "J"
      qn_id[2] = "F1"
      qn_id[3] = "F2"
      qn_id[4] = "F"
    if num_qn == 1:
      qn_id[0] = "J"
      qn_str="(J="+str(q_num[0])+")"
    elif mol_type == 2:
      # linear case b
      if deg_free == 2:
        qn_id[0] = "N"
        qn_id[1] = "\Lambda"
        qn_id[2] = "F1"
        qn_id[3] = "F2"
        qn_id[4] = "F"
      elif deg_free == 3:
        # symmetric rotor
        qn_id[0] = "N"
        qn_id[1] = "K"
        qn_id[2] = "J"
        qn_id[3] = "F1"
        qn_id[4] = "F2"
        qn_id[5] = "F"
        if num_qn == 2:
          qn_id[0] = "J"
    elif mol_type == 3:
      if deg_free == 2:
        # case a (2S+1 odd)
        qn_id[0] = "J"
        qn_id[1] = "\Omega"
        qn_id[2] = "\Lambda"
        qn_id[3] = "F1"
        qn_id[4] = "F2"
        qn_id[5] = "F"
      elif deg_free == 3:
        # asymmetric rotor
        qn_id[0] = "N"
        qn_id[1] = "K_{-1}"
        qn_id[2] = "K_{+1}"
        qn_id[3] = "J"
        qn_id[4] = "F1"
        qn_id[5] = "F"
        if num_qn == 3:
          qn_id[0] = "J"
    elif mol_type == 8:
      # linear case a (2S+1 even)
      qn_id[0] = "J+\frac{1}{2}"
      qn_id[1] = "\Omega+\frac{1}{2}"
      qn_id[2] = "\Lambda"
      qn_id[3] = "F1"
      qn_id[4] = "F2"
      qn_id[5] = "F"
  elif mol_type == 13:
    # symmetric rotor with vibration
    if tag == 32003: # methanol
      J = int(q_nums[0])
      K = int(q_nums[1])
      if K < 0:
        state = "A^-"
      elif K > 0:
        state = "A^+"
      else:
        state = "E"
      qn_str = str(J)+"_{"+str(abs(K))+"} "+state
      qn_id[0] = "J"
      qn_id[1] = "K_a"
      qn_id[2] = "K_c"
      qn_id[3] = "v_t"
    else:
      if diag:
        print "Processing symmetric rotor with vibration"
      qn_id[0] = "N"
      qn_id[1] = "K"
      qn_id[2] = "v"
      qn_id[3] = "J"
      qn_id[4] = "F1"
      qn_id[5] ="F"
      if num_qn == 3:
        qn_id[0] = "J"
  elif mol_type == 73:
    # ammonia with hyperfine splitting
    qn_id[0] = "J"
    qn_str[0] = str(q_nums[0])
    qn_id[1] = "K"
    qn_str[1] = str(q_nums[1])
    qn_id[2] = "V"
    if q_nums[2] == '0':
      qn_str[2] = "$0^+$"
    elif q_nums[2] == '1':
      qn_str[2] = "$0^-$"
    elif q_nums[2] == '2':
      qn_str[2] = "$\nu_2^+$"
    elif q_nums[2] == '3':
      qn_str[2] = "$\nu_2^-$"
    else:
      qn_str[2] = "??"
    if half_int_code & 4:
      qn_id[3] = "F_1+0.5"
      qn_str[3] = str(2*q_nums[3]-1)+"/2"
    else:
      qn_id[3] = "F_1"
      qn_str[3] = str(q_nums[3])
    if half_int_code & 2:
      qn_id[4] = "I_{tot}+0.5"
      qn_str[4] = str(2*int(q_nums[4])-1)+"/2"
    else:
      qn_id[4] = "I_{tot}"
      qn_str[4] = str(q_nums[4])
    if half_int_code & 1:
      qn_id[5] = "F+0.5"
      qn_str[5] = str(2*int(q_nums[5])-1)+"/2"
    else:
      qn_id[5] = "F"
      qn_str[5] = q_nums[5]
    qn_id.append("F2")
    if int(q_nums[1]) % 3 == 0:
      q_nums.append(str(0))
    elif int(q_nums[1]) % 3 == 1:
      q_nums.append(str(4))
    else:
      q_nums.append(str(2))
  elif mol_type ==14:
    # asymmetric rotor with vibration
    if diag:
      print "Processing asymmetric rotor with vibration"
    qn_id[0] = "N"
    qn_id[1] = "K-"
    qn_id[2] = "K+"
    qn_id[3] = "v"
    qn_id[4] = "J"
    qn_id[5] = "F"
    if num_qn == 4:
      qn_id[0] = "J"
  return qn_id, qn_str

def get_mol_metadata(tag):
  """
  Get the metadata for a molecule from the doc file
  """
  code_string = ('%06d' % tag)
  filename = 'd'+code_string+'.cat'
  # Use local version if available
  if exists(local_path+filename):
    fd = open(local_path+filename,'r')
  else:
    if not version.has_key(tag):
      v = latest_version(tag)+'/'
    elif version[tag] != '':
      v = version[tag]+''
    else:
      v = ''
    doc_path = jpl_url+v+'doc/'+filename
    print "Looking for", doc_path
    # Now process the doc file
    try:
      fd = urllib2.urlopen(doc_path)
    except urllib2.HTTPError, details:
      print details
      return
  lines = fd.readlines()
  fd.close()
  part_fn = {}
  for line in lines[:14]:
    data = line.lstrip('\\').strip().split('&')
    print data
    if re.search('Species',data[0]):
      name = data[-1].strip()
    if re.search('Lines Listed',data[0]):
      n_lines = int(data[1])
    if re.search('Q\(',data[2]):
      temp = float(data[2][3:8])
      Q = float(data[3])
      part_fn[temp] = Q
  return name, n_lines, part_fn
    
def get_mol_metadata2(tag):
  """
  Get the metadata for a molecule from the directory
  
  Get the general data for the molecule specified by tag,
  returning a sequence consisting of the name, the number of
  lines in the data file, and
  """
  q_temps = [300, 225, 150, 75, 37.5, 18.25, 9.375]
  doc_file = urllib2.urlopen(jpl_url+'catdir.cat','r')
  doc_text = doc_file.readlines()
  doc_file.close()
  name = ''
  for line in doc_text:
    mol_data = line.strip().split()
    id = int(mol_data[0])
    if id == tag:
      name = mol_data[1]
      n_lines = int(mol_data[2])
      partn_fn = {}
      for i in range(7):
        partn_fn[q_temps[i]] = mol_data[3+i]
      break
  if name != '':
    #print name, n_lines, partn_fn
    return name, n_lines, partn_fn
  else:
    return 'none',0,{}

def show_partn_fn(tag):
  """Make a plot of the partition function as a function of temperature
  for the molecule specified by tag."""
  metadata = get_mol_metadata(tag)
  Q_data = metadata[2]
  Q_temps = Q_data.keys()
  Q_temps.sort()
  Q_values = []
  for temp in Q_temps:
    Q_values.append(M.pow(10,float(Q_data[temp])))
  loglog(Q_temps,Q_values)
  axis([9,300,30,10000])
  title(metadata[0])
  xlabel('Temperature (K)',fontsize=12)
  ylabel('Partition Function')
  grid()
  show()

def einstein_a(intens,f,part_fn_300,g_upper,e_lower):
  """Returns Einstein A given the JPL Catalog intensity, the frequency in
  MHz, the 300 K partition function from the JPL Catalog, the degeneracy
  of the upper state, and the energy of the lower state in 1/cm."""
  # convert 1/cm to ergs
  e_lower *= c*h
  # calculate from e_lower and f, in ergs
  e_upper = e_lower + f*1000000.0*h
  # Eqn.6a in Catalogue explanatory supplement
  exp_term = M.exp(-e_lower/(k*300)) - M.exp(-e_upper/(k*300))
  return 10**intens * f**2 * (part_fn_300/g_upper) * 2.7964E-16/exp_term


def convert_JPL_QN(level):
  """Convert quantum numbers from first edition catalog to second edition
  version."""
  j = int(level[0])
  k = int(level[1])
  new_j = str(j)
  ka = str(k)
  kc = str(j - k)
  v = level[2]
  if v == '1':
    # A+ state
    if k%2 == 1:
      kc = str(j - k + 1)
    new_v = "+0"
  elif v == '2':
    # A- state
    if k%2 == 0:
      kc = str(j - k + 1)
    new_v = "-0"
  elif v == '3':
    # E1
    ka = "+"+str(k)
    new_v = "0"
  elif v == '4':
    # E2
    if k > 0:
      ka = str(-k)
    kc = str(j - k +1)
    new_v = "0"
  else:
    print "Unexpected value of v:",v,"in",level
  return [new_j,ka,kc,new_v]

def make_JPLv1_E_table():
  """Get a unique list sorted by increasing energy from the old catalog
  version and cross-reference it to the new version quantum labels.  Each
  item in the output list contains:
  sequence number, r1 quantum number list, energy, r2 quantum number list"""
  mol_fd = open("/home/kuiper/DOS/physics/radiat\'n/jpl_catl/c032003.cat","r")
  mol_data = mol_fd.readlines()
  mol_fd.close()
  old_transitions = []
  for line in mol_data:
    line_data = parse_catalog_line(line)
    old_transitions.append(line_data)
  old_levels = Radiation.Lines.Molec.find_levels(old_transitions,'q lower')
  old_sorted = sorted(old_levels,key=operator.itemgetter(1))
  # Compute new catalog quantum numbers
  index = 0
  for level in old_sorted:
    level.insert(0,index)
    new_QN = Radiation.Lines.Molec.jpl.convert_JPL_QN(level[1])
    level.append(new_QN)
    index += 1
  return old_sorted
