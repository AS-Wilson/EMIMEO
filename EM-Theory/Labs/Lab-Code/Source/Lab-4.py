import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scsp
from mpl_toolkits import mplot3d as plt3d
import os as sys
from sympy.functions.special.bessel import jn

# from matplotlib import rc

########################################
##### PROPAGATION CONSTANT CALCULATION #####

mag_permea = 1
elec_permit = 1
freq = 11.5e+09
c = 299792458

ang_freq = 2 * np.pi * freq
print("Angular Frequency is: ", ang_freq, "Rad/Sec")

k_0 = ang_freq / c
print("k_0 is: ", k_0)

k_c = np.sqrt( ((np.pi/(22.86e-03))**2) )
print("k_c is: ", k_c)

prop_const = np.sqrt( ((k_0**2) - (k_c**2)) )

print("The Propagation Constant is: ", prop_const, "Rad/m")


##### CUTOFF FREQUENCY CALCULATION #####
cutoff = (c / (2 * np.pi * np.sqrt(mag_permea / elec_permit))) * k_c
print("The Cutoff Frequency is: ", cutoff, " Hz")
print("The Cutoff Frequency is: ", np.round(cutoff/(10e+08), decimals=2), " GHz")

#def main():
    # step_index_graph()
    # normalised_spatial_freq_calcs()


#if __name__ == "__main__":
#    main()
    
    