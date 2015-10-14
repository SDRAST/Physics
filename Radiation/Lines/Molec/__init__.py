"""
Molec - Package to support molecular line spectroscopy
This provides support for the JPL Catalog at
http://spec/jpl.nasa.gov
and the Cologne Database for Molecular Spectroscopy at
http://www.astro.uni-koeln.de/site/vorhersagen/
Molecule tags with a 0 in the hundreds position are in the JPL catalog.
Tags with a 5 in that position are in the Cologne catalog.

Generally, both catalogs use the same way to encode line data.  The
exception is for lines with the quantum number format code of 999.
In the Cologne directory, the tags for these molecules are followed
immediately by *. The Cologne catalog also appends the molecule name
on each line.

References::
  H. M. Pickett, R. L. Poynter, E. A. Cohen, M. L. Delitsky, J. C. Pearson,
    and H. S. P. Mueller,
    J. Quant. Spectrosc. Radiat. Transfer  60 (1998) 883 - 890.
  H. S. P. Mueller, F. Schloeder, J. Stutzki, and G. Winnewisser,
    J. Mol. Struct. 742, 215-227 (2005)
  H. S. P. Mueller, S. Thorwirth, D. A. Roth, and G. Winnewisser,
    Astronomy and Astrophysics 370, L49-L52 (2001)
"""

import urllib
import Physics as P
import math as M
from Physics.Radiation.Lines.Molec import koln
from Physics.Radiation.Lines.Molec import jpl

MHz2invCm = (P.pq(1e6,'Hz')/P.pq(1,'c')).inUnitsOf('1/cm').value

cat_path = '/usr/local/line_cat/'
jpl_url  = 'http://spec.jpl.nasa.gov/ftp/pub/catalog/'
koln_url = 'http://www.astro.uni-koeln.de/site/vorhersagen/catalog/'

def mol_file_basename(tag):
    """Convert a tag into a molecular data file name.
    This could be a tag string, if it is quoted and starts with a 0,
    e.g. '032504', or a molecule number (int), e.g. 32504"""
    if str(tag)[0] == '0':
      return 'c'+tag+'.'
    else:
      return 'c%06d.' % tag

def open_molecule(code):
    """Returns a file descriptor for the data file for the molecule
    with the specified code."""
    filename = mol_file_basename(code)+'cat'
    try:
      # Has the file been stored locally?
      if filename[4:5] == '0':
        mol_path = cat_path+'jpl/'+filename
      elif filename[4:5] == '5':
        mol_path = cat_path+'koln/'+filename
      else:
        return None
      mol_fd = open(mol_path)
    except:
      # See which catalog it is in
      if filename[4:5] == '0':
        mol_path = jpl_url+filename
      elif filename[4:5] == '5':
        mol_path = koln_url+filename
      else:
        return None
      print mol_path
      mol_fd = urllib.urlopen(mol_path)
    return mol_fd

def extract_lines(molecule,f_low=0,f_up=1000000):
    """Extract all the lines in the specified frequency range. 'molecule'
    can be either a list of lines extracted from a catalog or a tag for the
    designated molecule in the appropriate line catalog. Frequencies are in
    MHz.  The default is to extract all lines."""
    if isinstance(molecule,list):
      mol_data = molecule
    else:
      mol_fd = open_molecule(molecule)
      mol_data = mol_fd.readlines()
      mol_fd.close()
    result = []
    if len(mol_data) == 0:
      # return an empty list
      return result
    # check to see which format to be parsed
    if mol_data[0].split()[7] == '999':
      parser = koln.parse_catalog_line
    else:
      parser = jpl.parse_catalog_line
    for line in mol_data:
        line_data = parser(line)
        if ((float(line_data['freq']) > float(f_low)) & \
            (float(line_data['freq']) < float(f_up) ) ):
                result.append(line)
    return(result)

def find_levels(transitions,key='both'):
  """This extracts energy levels from a list of transitions according to a key
  such as 'q lower' or 'q upper' or some other key.  Redundancy is eliminated.
  The list items contain the quantum numbers and the energy if the key is
  'q upper' or 'q lower'.  Otherwise, the list just contains the selected key
  values."""
  levels = []
  for transition in transitions:
    found = False
    for level in levels:
      # Check the QN to see if the level was already found.
      if key == 'both':
        if (level[0] == transition['q lower']) or \
           (level[0] == transition['q upper']):
          found == True
          break
      elif level[0] == transition[key]:
        found = True
        break
    if found == False:
      # Add the state to the list of states
      Elo = float(transition['E lower'])
      Ehi = Elo + float(transition['freq'])*MHz2invCm
      if key == 'both':
        levels.append([transition['q lower'],Elo])
        levels.append([transition['q upper'],Ehi])
      elif key == 'q lower':
        levels.append([transition[key],Elo])
      elif key == 'q upper':
        levels.append([transition[key],Ehi])
      else:
        levels.append(transition[key])
  return levels
