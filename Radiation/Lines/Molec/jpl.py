# -*- coding: latin-1 -*-
"""
module jpl - Routines to access the JPL Line Catalog

Reference
=========
http://spec.jpl.nasa.gov/
http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/catintro.pdf
"""
import logging
import urllib.request, urllib.error, urllib.parse
from Physics import *
import Math as M
import operator
from support import text
import re
from os.path import exists

logging.basicConfig(level=logging.WARNING)
module_logger = logging.getLogger(__name__)

diag = False

jpl_url = 'http://spec.jpl.nasa.gov/ftp/pub/catalog/'
local_path = '/usr/local/line_cat/jpl/'

# This keeps track of where to find the molecular data
version = {}

def catalog_versions():
  """
  Get the catalog versions available and save as a module global
  """
  global version
  # First find the available versions
  fd = urllib.request.urlopen(jpl_url)
  lines = fd.readlines()
  fd.close()
  versions = []
  for line in lines:
    txt = text.remove_html_tags(line.strip())
    if txt != '' and re.search('v',txt):
      parsed = txt.strip() # split('/')
      module_logger.debug("catalog_versions: %s found in '%s'", parsed, line)
      if parsed[1].isdigit():
        versions.append(parsed[1])
  # Now check each version's subdirectory, last first
  versions.sort(reverse=True) # should not be needed
  module_logger.info("catalog_versions: %s", versions)
  return versions

def get_codes(first_chars=None):
  """
  Get the codes for the available molecules
  
  This returns a dict of the form::
     {'OD': '18001', 'OH': '17001', 'PH3': '34003', 'PN': '45013', 
      'PO': '47006', 'PO+': '47005', 'PO2': '63008', 'PS': '63007',
  If first_chars is given, then only molecules whose names begin with those
  characters are returned.  first_chars may be a regular expression
  """
  doc_file = urllib.request.urlopen(jpl_url+'catdir.cat','r')
  doc_text = doc_file.readlines()
  doc_file.close()
  molecules = {}
  for line in doc_text:
    mol_data = line.strip().split()
    if first_chars:
      if re.match(first_chars,mol_data[1]):
        molecules[mol_data[1]] = mol_data[0]
    else:
      molecules[mol_data[1]] = mol_data[0]
  return molecules
  
def latest_version(code):
  """
  Returns the latest version for a given molecule
  """
  code_string = ('%06d' % code)
  filename = 'c'+code_string+'.cat'
  for v in catalog_versions():
    module_logger.debug("latest_version: Trying %s%s/%s",jpl_url, v, filename)
    try:
      mol_fd = urllib.request.urlopen(jpl_url+v+'/'+filename)
    except urllib.error.HTTPError:
      continue
    if code not in version:
      version[code] = v
    mol_fd.close()
    return v
  # Now try original catalog
  module_logger.debug("latest_version: trying %s%s",jpl_url,filename)
  try:
    mol_fd = urllib.request.urlopen(jpl_url+filename)
  except urllib.error.HTTPError:
    return None
  if code not in version:
    version[code] = v
  mol_fd.close()
  return v
  
def open_molecule(code, v=None):
  """
  Returns a file descriptor for the data file for the molecule
  with the specified code.

  Catalog Version
  ===============
  It looks for the latest version in the main catalog.  If a specific version
  is wanted, it looks in the sub-directory for that version.

  @type code : int
  @param code : molecule identifier
  """
  code_string = ('%06d' % code)
  filename = 'c'+code_string+'.cat'
  #v = latest_version(code)
  if v:
    catpath = jpl_url+'v'+v+'/'+filename
  else:
    catpath = jpl_url+filename
  module_logger.info("open_molecule: opening %s", catpath)
  try:
    mol_fd = urllib.request.urlopen(catpath)
    return mol_fd
  except urllib.error.HTTPError as details:
    module_logger.error("open_molecule: could not open %s", catpath)
    module_logger.error("open_molecule: %s", details)
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
  data['freq'] = float(line[0:13].strip())
  data['freq err'] = float(line[13:21].strip())
  data['int'] = float(line[21:29].strip())
  data['deg free'] = int(line[29:31].strip())
  data['E lower'] = float(line[31:41].strip())
  data['g upper'] = int(line[41:44].strip())
  data['tag'] = int(line[45:51])
  data['qn format'] = int(line[51:55])
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

def extract_range(moltag, f_low, f_up, min_int=-5):
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
  mol_fd = open_molecule(moltag)
  mol_data = mol_fd.readlines()
  mol_fd.close()
  result = []
  for line in mol_data:
    line_data = parse_catalog_line(line)
    module_logger.debug('extract_range: got %s', line_data)
    if ( ((line_data['freq']) > f_low)   & \
         ((line_data['freq']) < f_up)    & \
         ((line_data['int'] ) > min_int)   ):
            result.append(line_data)
    elif line_data['freq'] > f_up:
      break
  return(result)

def quantum_label(tag, q_nums, q_fmt, deg_free):
  """
  Decode the quantum numbers in the JPL catalog.

  References
  ==========
  Refer to section V of the JPL Submm Catalog Supplement
  and Townes and Schawlow, Dover edition 1975, p.176.
  This version follows the JPL convention in use on 2008 Jun 15

  Quantum number format
  =====================
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
  def decode_K(string):
    """
    Decode negative quantum numbers greater than -9
    """
    try:
      return int(string)
    except ValueError:
      prefix = ord(string[0])
      if prefix > 64 and prefix < 91:
        # upper case for 10, 11, 12, etc.
        return 10*(prefix - 55) + int(string[1])
      elif prefix > 96 and prefix < 123:
        # lower case for -1, -2, -3, etc.
        return 10*(96 - prefix) - int(string[1])
      else:
        module_logger.error("quantum_label.decode_K: invalid code %s",
                            string[0])
        raise ValueError
    
  module_logger.debug(" quantum_label: species tag = %d", tag)
  module_logger.debug(" quantum_label: QNs: %s", q_nums)
  quantum_code = str(q_fmt)
  mol_type = q_fmt/100
  module_logger.debug(" quantum_label: molecule type = %d", mol_type)
  # mol_type MOD 5 is the number of primary quantum numbers
  num_primary_qn = mol_type % 5
  module_logger.debug(" quantum_label: number of primary QNs = %d",
                      num_primary_qn)
  # half integer code
  half_int_code = (q_fmt-mol_type*100)/10
  module_logger.debug(" quantum_label: half integer code = %d", half_int_code)
  # number of quantum numbers for each state.
  num_qn = int(quantum_code[-1])
  module_logger.debug(" quantum_label: number of QNs = %d", num_qn)
  
  qn_id = ["","","","","",""]
  qn_str = ["","","","","",""]
  if mol_type == 0:
    # atoms:  (J),(F),· · ·
    module_logger.debug(" quantum_label: Processing atom")
    #    test half integer code
    # defaults
    qn_id[0] = "J"; qn_id[1] = "F"
    if half_int_code & 1: # if half_int_code[2] == '1':
      qn_id[0] = "J+0.5"
      qn_str[0] = str(2*q_nums[0]-1)+"/2"
    #    test first quantum number
    if q_nums[0] == 0:
      qn_id[1] = ""
    elif q_nums[0] == 1:
      qn_id[1] = "F+.5"
      qn_str[1] = "(F="+str(2*decode_K(q_nums[1])-1)+"/2)"
    elif q_nums[0] == 2:
      qn_id[1] = "F"
      qn_str[1] = "(F="+str(decode_K(q_nums[1]))+")"
  elif mol_type == 1 or mol_type == 2 or mol_type == 3 or mol_type == 8:
    module_logger.debug(" quantum_label: linear or diatomic molecule")
    # linear and diatomic molecules
    #    test molecular type code
    if mol_type == 1: # N, (J), (F1), (F2)(F)
      # sigma state
      # defaults
      qn_id[0] = "N"
      qn_id[1] = "J"
      qn_id[2] = "F1"
      qn_id[3] = "F2"
      qn_id[4] = "F"
      if num_qn == 1:
        qn_id[0] = "J"
        qn_str="(J="+str(q_nums[0])+")"
      elif num_qn == 2:
        qn_id[0] = "J"
        qn_id[1] = "K"
        qn_str="(J,K="+str(q_nums[0])+","+str(q_nums[1])+")"
    elif mol_type == 2:
      # linear case b
      if deg_free == 2: #  N, Λ, (F1), (F2), (F)
        qn_id[0] = "N"
        qn_id[1] = "\Lambda"
        qn_id[2] = "F1"
        qn_id[3] = "F2"
        qn_id[4] = "F"
      elif deg_free == 3: # N, K, (J), (F1), (F2), (F)
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
      if deg_free == 2: #  J, Ω, Λ, (F1), (F2), (F)
        # case a (2S+1 odd)
        qn_id[0] = "J"
        qn_id[1] = "\Omega"
        qn_id[2] = "\Lambda"
        qn_id[3] = "F1"
        qn_id[4] = "F2"
        qn_id[5] = "F"
      elif deg_free == 3: # N, K−1, K+1, (J), (F1), (F)
        # asymmetric rotor
        qn_id[0] = "N"
        qn_id[1] = "K_{-1}"
        qn_id[2] = "K_{+1}"
        qn_id[3] = "J"
        qn_id[4] = "F1"
        qn_id[5] = "F"
        if num_qn == 3:
          qn_id[0] = "J"
    elif mol_type == 8: # J+1/2, Ω+1/2, Λ, (F1), (F2), (F)
      # linear case a (2S+1 even)
      qn_id[0] = "J+\frac{1}{2}"
      qn_id[1] = "\Omega+\frac{1}{2}"
      qn_id[2] = "\Lambda"
      qn_id[3] = "F1"
      qn_id[4] = "F2"
      qn_id[5] = "F"
  elif mol_type == 13: # N, K, v, (J), (F1), (F)
    # symmetric rotor with vibration
    if tag == 32003: # methanol
      # The quantum number format is J, ±K, parity (±), vt. A-state transitions
      # have a ’+’ or ’-’ in the parity column while E states do not. E1 and E2
      # are denoted by the positive (unsigned) and negative sign of K, 
      # respectively, in the absence of a parity entry.
      J = int(q_nums[0])
      K = decode_K(q_nums[1])
      p = q_nums[2]
      if p == "-": # K < 0:
        state = "A^-"
      elif p == "+": # > 0:
        state = "A^+"
      else:
        if K < 0:
          state = "E1"
        else:
          state = "E2"
      qn_str = str(J)+"_{"+str(abs(K))+"} "+state
      qn_id[0] = "J"
      qn_id[1] = "K_a"
      qn_id[2] = "K_c"
      qn_id[3] = "v_t"
    else: 
      module_logger.debug(" quantum_label: symmetric rotor with vibration")
      qn_id[0] = "N"
      qn_id[1] = "K"
      qn_id[2] = "v"
      qn_id[3] = "J"
      qn_id[4] = "F1"
      qn_id[5] ="F"
      if num_qn == 3:
        qn_id[0] = "J"
      K = abs(decode_K(q_nums[1]))
      qn_str = "J_K = "+q_nums[0]+"_{"+str(K)+"} v="+q_nums[2]
  elif mol_type == 73: # 
    # ammonia with hyperfine splitting
    qn_id[0] = "J"
    qn_str[0] = str(q_nums[0])
    qn_id[1] = "K"
    qn_str[1] = str(decode_K(q_nums[1]))
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
    if half_int_code & 4: # if half_int_codes[0] == '1':
      qn_id[3] = "F_1+0.5"
      qn_str[3] = str(2*q_nums[3]-1)+"/2"
    else:
      qn_id[3] = "F_1"
      qn_str[3] = str(q_nums[3])
    if half_int_code & 2: # if half_int_code & 2:
      qn_id[4] = "I_{tot}+0.5"
      qn_str[4] = str(2*int(q_nums[4])-1)+"/2"
    else:
      qn_id[4] = "I_{tot}"
      qn_str[4] = str(q_nums[4])
    if half_int_code & 1: # if half_int_code & 1:
      qn_id[5] = "F+0.5"
      qn_str[5] = str(2*int(q_nums[5])-1)+"/2"
    else:
      qn_id[5] = "F"
      qn_str[5] = q_nums[5]
    qn_id.append("F2")
    if decode_K(q_nums[1]) % 3 == 0:
      q_nums.append(str(0))
    elif decode_K(q_nums[1]) % 3 == 1:
      q_nums.append(str(4))
    else:
      q_nums.append(str(2))
  elif mol_type ==14: # N, K−1, K+1, v, (J), (F)
    # asymmetric rotor with vibration
    module_logger.debug(" quantum_label: asymmetric rotor with vibration")
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
  
  Get the general data for the molecule specified by tag,
  returning a sequence consisting of the name, the number of
  lines in the data file, and the partition function.
  """
  code_string = ('%06d' % tag)
  filename = 'd'+code_string+'.cat'
  # Use local version if available
  if exists(local_path+filename):
    fd = open(local_path+filename,'r')
  else:
    #if not version.has_key(tag):
    #  v = latest_version(tag)+'/'
    #elif version[tag] != '':
    #  v = version[tag]+''
    #else:
    #  v = ''
    #doc_path = jpl_url+v+'doc/'+filename
    doc_path = jpl_url+'doc/'+filename
    #module_logger.debug("get_mol_metadata: looking for %s", doc_path)
    # Now process the doc file
    try:
      fd = urllib.request.urlopen(doc_path)
    except urllib.error.HTTPError as details:
      module_logger.error("get_mol_metadata: failed %s", details)
      return None
  lines = fd.readlines()
  fd.close()
  part_fn = {}
  for line in lines[:14]:
    data = line.lstrip('\\').strip().split('&')
    #module_logger.debug('get_mol_metadata: processing: %s', data)
    if len(data) > 1 or data[0]:
      # Ignore blank lines
      if re.search('Species',data[0]):
        name = data[-1].strip()
      if re.search('Lines Listed',data[0]):
        n_lines = int(data[1])
      if re.search('Q\(',data[2]):
        temp = float(data[2][3:8])
        Q = float(data[3])
        part_fn[temp] = M.log10(Q)
  return name, n_lines, part_fn
    
def get_mol_metadata2(tag):
  """
  Get the metadata for a molecule from the directory
  
  Get the general data for the molecule specified by tag,
  returning a sequence consisting of the name, the number of
  lines in the data file, and the partition function, same
  as get_mol_metadata()
  """
  q_temps = [300, 225, 150, 75, 37.5, 18.25, 9.375]
  doc_file = urllib.request.urlopen(jpl_url+'catdir.cat','r')
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

def einstein_a(intens, f, part_fn_300, g_upper, e_lower):
  """
  Returns Einstein A given the JPL Catalog intensity
  
  @param f : the frequency in MHz
  
  @param part_fn_300 : log10 of the 300 K partition function
  
  @param g_upper : the degeneracy of the upper state
  
  @param e_lower : the energy of the lower state in 1/cm.
  """
  # convert 1/cm to ergs
  e_lower *= c*h
  # calculate from e_lower and f, in ergs
  e_upper = e_lower + f*1000000.0*h
  # Eqn.6a in Catalogue explanatory supplement
  exp_term = M.exp(-e_lower/(k*300)) - M.exp(-e_upper/(k*300))
  return 10**intens * f**2 * (part_fn_300/g_upper) * 2.7964E-16/exp_term

def convert_JPL_QN(level):
  """
  Convert quantum numbers from first edition catalog to second edition version.
  """
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
    print("Unexpected value of v:",v,"in",level)
  return [new_j,ka,kc,new_v]

def make_JPLv1_E_table(species):
  """
  
  Make an energy level table
  
  Get a unique list sorted by increasing energy from the old catalog
  version and cross-reference it to the new version quantum labels.  Each
  item in the output list contains:
  sequence number, r1 quantum number list, energy, r2 quantum number list
  """
  filename = "c%06d.cat" % species
  mol_fd = open("/home/kuiper/DOS/physics/radiat\'n/jpl_catl/"+filename,"r")
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
