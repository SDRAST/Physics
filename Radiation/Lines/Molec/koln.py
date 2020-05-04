"""
koln_cat.py - Routines for the Cologne Database for Molecular Spectroscopy.

The database has one main format format, described in
http://www.astro.uni-koeln.de/site/vorhersagen/catalog/main_catalog.shtml.
Most datafiles follow the JPL Spectral Line Catalog format.  However, there
are exceptions. For example, the light hydrides have frequencies cm^-1
instead of MHz.  For methanol, the intensity entry is
log_10  (S mu^2 ) in D^2  and the quantum numbers have special meanings.
"""

import urllib.request, urllib.parse, urllib.error
import operator

cat_path = '/usr/local/line_cat/koln/'
koln_url = 'http://www.astro.uni-koeln.de/site/vorhersagen/catalog/'

        
def parse_catalog_line(line):
    """
    This parses a line of data extracted from a molecular data file.
    It returns a dictionary with fairly obvious keys.
    For the encoding details see
    http://www.astro.uni-koeln.de/site/vorhersagen/catalog/main_catalog.shtml\
      #description
    """
    line_data = line.split()
    data = {}
    data['freq'] = float(line_data[0])          # MHz
    data['freq err'] = float(line_data[1])      # MHz
    data['str'] = pow(10,float(line_data[2]))   # S mu^2, the line strength
    data['deg free'] = int(line_data[3])        # Degrees of freedom in the
                                                # rotational partition function
                                                #(0 for atoms,
                                                # 2 for linear molecules, and
                                                # 3 for nonlinear molecules).
    data['E lower'] = float(line_data[4])       # cm^{-1}
    data['g upper'] = int(line_data[5])         # for upper state
    data['tag'] = int(line_data[6])             # molecule code
    data['qn format'] = int(line_data[7])
    q_low = []
    q_up = []
    if data['qn format'] == 999:
      # in this case, the string range [55:] has eight fields of three
      for i in range(4):
        ptr_up = 55+3*i
        ptr_low = 55+12+3*i
        q_up.append(line[ptr_up:ptr_up+3].strip())
        q_low.append(line[ptr_low:ptr_low+3].strip())
      data['q upper'] = q_up
      data['q lower'] = q_low
    else:
      # JPL catalog format; the range [55:] has twelve fields of two
      for i in range(6):
        ptr_up = 55+2*i
        ptr_low = 55+12+2*i
        q_up.append(line[ptr_up:ptr_up+2].strip())
        q_low.append(line[ptr_low:ptr_low+2].strip())
      data['q upper'] = q_up
      data['q lower'] = q_low
    return(data)

def get_mol_metadata(tag):
    """
    Get the general data for the molecule specified by tag,
    returning a sequence consisting of the name, the number of lines in the
    data file, and a dictionary with partition function data keyed to
    temperatures
    """
    q_temps = [300, 225, 150, 75, 37.5, 18.25, 9.375]
    # Do we have a local copy?
    try:
      doc_file = open(cat_path+'partition_function.dat','r')
      doc_text = doc_file.readlines()[2:]
    except:
      doc_file = urllib.request.urlopen(koln_url+'partition_function.html','r')
      doc_text = doc_file.readlines()[14:]
    doc_file.close()
    for line in doc_text:
      print(line)
      id = int(line[:6])
      name = line[7:30].strip()
      mol_data = line[31:].strip().split()
      print(mol_data)
      if id == tag:
        n_lines = mol_data[0]
        partn_fn = {}
        for i in range(7):
          partn_fn[q_temps[i]] = mol_data[3+i]
        break
    if name != '':
      print(name, n_lines, partn_fn)
      return name, n_lines, partn_fn
    else:
      return 'none',0,{}

def einstein_a(freq,line_str,g_up):
    """
    Einstein A from line_str (S_g mu_g^2), freq (MHz) and g_up, as
    per eq. 6 in
    www.astro.uni-koeln.de/site/vorhersagen/catalog/main_catalog.shtml
    """
    return 1.16395e-20 * freq**3 * line_str / g_up


def make_Koln_E_table():
  """
  Get a unique list sorted by increasing energy from the Koln catalog.
  Each item in the output list contains:
  sequence number, quantum number list, energy
  """
  transitions = extract_lines(32504)
  levels = find_levels(transitions)
  sorted_levels = sorted(levels,key=operator.itemgetter(1))
  index = 0
  for level in sorted_levels:
    level.insert(0,index)
    index += 1
  return sorted_levels
