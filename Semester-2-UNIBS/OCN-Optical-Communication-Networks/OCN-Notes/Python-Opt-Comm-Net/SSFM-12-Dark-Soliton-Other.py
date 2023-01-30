#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
from OCNToolbox import Split_Step_Fibre_Propagation, gaussian, make_3d_graph

#####################################################################################################################
##### IGNORE THIS IF YOU WISH, THIS IS SIMPLY TO MAKE PRETTY GRAPHS #####

##### matplotlib graph settings #####
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


#####################################################################################################################
no_of_samples = 1000000
signal_duration = 10

# Create an array of time values with X number of equally-spaced values, a larger array size for a given
#   signal duration makes for a better frequency domain transform
time = np.linspace(-signal_duration, signal_duration, no_of_samples)

eta = 1
t_shift = 0

##########
B1 = 0.5
##### Method laid out in Nonlinear Fibre Optics, Agrawl, 2013 (Fifth Edition). See Page 152, Equation 5.3.2:
## This method gives the correct shape graph but, the phase is shifted lower than it should be and the intensity is
## shifted to the right in the time domain.
# sigma1 = eta * B1 * (time - t_shift - (eta * B1 * np.sqrt(1 - (B1**2))))
# B1_input = eta * (B1 * np.tanh(sigma1) - complex(0, np.sqrt(1 - B1**2)))
# B1_intensity = abs(B1_input)**2
# B1_phase = np.angle(B1_input)

##### Method laid out in Nonlinear Fibre Optics, Agrawal, 2001 (Third Edition). See Page 160, Equation 5.3.2:
## The intensity graph can be obtained by considering equation V(tau) must be squared, cannot get phase graphed at all
B1_intensity = eta * (1 - (B1**2 * (1 - ( (np.tanh(eta * B1 * (time - t_shift)))*(np.tanh(eta * B1 * (time - t_shift))) ) )))
B1_phase = np.exp(((eta * np.sqrt(1 - B1**2) * time) + np.arctan(((B1 * np.tanh(eta * B1 * time)) / (np.sqrt(1 - B1**2))))).astype(complex))


##########
B2 = 0.8

# sigma2 = eta * B2 * (time - (eta * B2 * np.sqrt(1 - (B2**2))))
# B2_input = eta * (B2 * np.tanh(sigma2) - complex(0, np.sqrt(1 - B2**2))) * np.exp(complex(0, eta**2 * 0))
# B2_intensity = abs(B2_input)**2
# B2_phase = np.arctan2(B2_input.imag, B2_input.real)

# B2_intensity = eta * (1 - (B2**2 * (1 - ( (np.tanh(eta * B2 * time))*(np.tanh(eta * B2 * time)) ) )))
# B2_phase = (eta * np.sqrt(1 - B2**2) * time) + np.arctan((B2 * np.tanh(eta * B2 * time)) / (np.sqrt(1 - B2**2)))


##########
B3 = 1

# sigma3 = eta * B3 * (time - (eta * B3 * np.sqrt(1 - (B3**2))))
# B3_input = eta * (B3 * np.tanh(sigma3) - complex(0, np.sqrt(1 - B3**2)))
# B3_intensity = abs(B3_input)**2
# B3_phase = np.angle(B3_input)

# B3_intensity = eta * (1 - (B3**2 * (1 - ( (np.tanh(eta * B3 * time))*(np.tanh(eta * B3 * time)) ) )))
# B3_phase = (eta * np.sqrt(1 - B3**2) * time) + np.arctan((B3 * np.tanh(eta * B3 * time)) / (np.sqrt(1 - B3**2)))






#####################################################################################################################
########## GRAPHING ##########

## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the intensity envelope vs time
plt.plot(time, B1_intensity, color='blue', label='B = 0.5')
# plt.plot(time, B2_intensity, color='red', label='B = 0.8')
# plt.plot(time, B3_intensity, color='green', label='B = 1')


# Titles and labels
plt.title(r'Input Envelope Intensity, Grey and Dark Solitons - Time Domain')
plt.ylabel(r'Intensity, $|F(t)|^2$')
plt.xlabel(r'Time, t')
plt.legend(loc="upper left")

# Set limits
plt.xlim(-4, 4)
plt.ylim(-0.01, 1.2)

# Aspect Ratio
axes = plt.gca()
# axes.set_aspect(6)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# Save graph in this location
# plt.savefig('../Graphs/split-step-2d-fund-bright-soliton-time.pdf', bbox_inches='tight', pad_inches=0)





## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the intensity envelope vs time
plt.plot(time, B1_phase, color='blue', label='B = 0.5')
# plt.plot(time, B2_phase, color='red', label='B = 0.8')
# plt.plot(time, B3_phase, color='green', label='B = 1')


# Titles and labels
plt.title(r'Input Envelope Phase, Grey and Dark Solitons - Time Domain')
plt.ylabel(r'Phase, $\phi$')
plt.xlabel(r'Time, t')
plt.legend(loc="upper left")

# Set limits
plt.xlim(-4, 4)
plt.ylim(-2, 2)

# Aspect Ratio
axes = plt.gca()
# axes.set_aspect(6)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# Save graph in this location
# plt.savefig('../Graphs/split-step-2d-fund-bright-soliton-time.pdf', bbox_inches='tight', pad_inches=0)




plt.show()


