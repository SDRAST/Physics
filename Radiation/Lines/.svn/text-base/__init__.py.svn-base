"""Molec -
Package to support molecular line spectroscopy

This uses two compilations, those at JPL in
http://spec.jpl.nasa.gov
and those at the University of Cologne:
http://www.astro.uni-koeln.de/site/vorhersagen/catalog/main_catalog.shtml

In the JPL catalog, the hundreds digit of the molecule tag is 0.  In the
Cologne catalog, it is 5."""

def mol_file_basename(tag):
    """Covert a tag into a molecular data file name"""
    return 'c%06d.' % tag
