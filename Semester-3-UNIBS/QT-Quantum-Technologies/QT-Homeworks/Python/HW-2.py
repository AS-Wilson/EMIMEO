#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.special import jv, kn, jve, kve
from sympy import Point, Line, Segment
from mpl_toolkits.mplot3d import Axes3D

from engineering_notation import EngNumber
import siunits as u

#####################################################################################################################
##### matplotlib graph settings #####
# Produce all graphs as PGFs so that they can be imported natively to latex
# and still be edited if needed
plt.rcParams["pgf.texsystem"] = "pdflatex"
# Setting fonts for all the graphs etc.
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "ptserif"
plt.rcParams["font.size"] = 9

# Match the colouring if my Latex docs
COLOUR = '#2E2E2E'
plt.rcParams['text.color'] = COLOUR
plt.rcParams['axes.labelcolor'] = COLOUR
plt.rcParams['xtick.color'] = COLOUR
plt.rcParams['ytick.color'] = COLOUR
plt.rcParams['pgf.preamble'] = r"\usepackage[T1]{fontenc} \usepackage{mathpazo}"

# Set style to bmh and place the ticks for the axes on the outside to be readable
plt.style.use('bmh')
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'


##### Universal Constsnt Definitions #####
planck_const_J_s = 6.626e-34
planck_const_eV_s = 4.136e-15
celerity = 2.998e8


##### Question Three #####

power = 1e-3
radius = 1.1e-3/2
wavelength_green = 532e-9

SA = np.pi * radius**2
print('Surface Area is : ', EngNumber(SA),  u.m**2)

flux = power / SA
print('Flux is : ', EngNumber(flux), u.j / (u.m**2 * u.s), '\n')

flux_eV = flux / 1.602e-19

green_energy_eV = (planck_const_eV_s * celerity) / wavelength_green
print('Photon Energy is : ', EngNumber(green_energy_eV), 'electronVolts')

green_energy_J = (planck_const_J_s * celerity) / wavelength_green
print('Photon Energy is : ', green_energy_J, u.j, '\n')

q_flux = flux_eV / green_energy_eV
print('Quantum Flux is : ', q_flux, ' Quanta', u.m**-2 * u.s**-1)

q_flux = flux / green_energy_J
print('Quantum Flux is : ', q_flux, ' Quanta', u.m**-2 * u.s**-1, '\n')


q_flux_10 = q_flux * 10 * SA

print('Quantum Flux in 10 sec is : ', q_flux_10, ' Quanta')


half_q_flux_10 = q_flux_10 / 2
print('Quantum Flux in 10 sec with half covered hole is : ', half_q_flux_10, ' Quanta')
