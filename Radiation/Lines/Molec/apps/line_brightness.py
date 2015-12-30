"""
Plots the LTE intensity of spectral lines of a molecule
"""

import logging
from pylab import *

from Astronomy import pc
from Physics.Radiation import LTE_pop_ratio, n_upper_LTE
from Physics.Radiation import extinction_coefficient_MKS, optical_depth
from Physics.Radiation import source_function_K, brightness_temperature
from Physics.Radiation.Lines.Molec.jpl import extract_range, parse_catalog_line
from Physics.Radiation.Lines.Molec.jpl import einstein_a, quantum_label
from Physics.Radiation.Lines.Molec.jpl import get_mol_metadata

logging.basicConfig(level=logging.WARNING)
mylogger = logging.getLogger()
mylogger.setLevel(logging.DEBUG)

#species = 17002;  T_ks = [10, 20, 30]; n_species = 3e-4      # NH3
#species = 99002;   T_ks = [8, 10, 15];  n_species = 1e-4       # HC7N
#species = 147001; T_ks = [8, 10, 15]; n_species = 1e-4       # HC11N
species = 32003;  T_ks = [20, 40, 80, 160]; n_species = 3e-3 # CH3OH

 # K
colors = ['r','g','b','y']
L_path = 1 # pc
fmin = 17000
fmax = 27000

moldata = get_mol_metadata(species)
name = moldata[0]
if name == 'CH3OH':
  name = 'CH$_3$OH'
Q_300 = 10**moldata[2][300.0]
mylogger.info(" Q(300) = %f", Q_300)

# required arguments are moltag, f_low, f_up
transitions = extract_range(species, fmin, fmax, min_int=-7)
mylogger.info(" %d transitions found", len(transitions))

TBmin = 100
TBmax = 0
color_index = 0
for T_k in T_ks:
  freqs = []
  TBs = []
  for transition in transitions:
    label = quantum_label(species, transition['q upper'],
                                   transition['qn format'],
                                   transition['deg free'])[1]
    label += ' - '
    label += quantum_label(species, transition['q lower'],
                                    transition['qn format'],
                                    transition['deg free'])[1]
    mylogger.info(" transition is %s", label)

    # required arguments are transition,temp
    R = LTE_pop_ratio(transition, T_k, degenerate=False)
    mylogger.info(" n_u/n_l = %f", R)

    # required arguments are species_density, transition, kinetic temp
    n_u = n_upper_LTE(n_species, transition, T_k)
    mylogger.info(" n_u = %e", n_u)

    # required arguments: log10(intens), freq (MHz), part_fn_300, g_upper, e_lower
    freq = transition['freq']
    freqs.append(freq)
    A = einstein_a(transition['int'],
                   freq,
                   Q_300, 
                   transition['g upper'],
                   transition['E lower'])
    mylogger.info(" A_ul = %e", A)

    # required arguments: n_u (cm^-3), pop_ratio, freq (Hz), Aul (1/s), dv (m/s)
    # MKS because the speed of light is in m/s
    K = extinction_coefficient_MKS(n_u*1e6, R, freq*1e6, A, 1e3)/100
    mylogger.info(" extinction coefficient = %e", K)

    # required arguments are extinc_coef (1/cm), path_length (cm)
    # CGS extinction is (1/100)^3 of MKS extinction
    tau = optical_depth(K, L_path*pc*100)
    mylogger.info(" tau = %f", tau)

    # required arguments are freq (Hz), population_ratio
    S = source_function_K(freq*1e6, R)
    mylogger.info(" source function = %f K", S)

    T = brightness_temperature(S, tau)
    TBmin = min(TBmin, T)
    TBmax = max(TBmax, T)
    TBs.append(T)
    mylogger.info(" brightness temperature = %f", T)
    print "%s %s T_B = %5.2f" % (name, label, T)

  scatter(freqs, TBs, label=str(T_k)+" K", c=colors[color_index])
  color_index += 1
TBmin = max(TBmin,0.01)
xlim(fmin,fmax)
xlabel('Frequency (MHz)')
yscale('log')
ylim(TBmin*0.95,TBmax*1.05)
ylabel('T$_B$ (K)')
legend(loc='lower left', numpoints=1)
grid(True)
title(name+(' (%7.1e $cm^{-2}$)'% (n_species*L_path*pc)))
show()

