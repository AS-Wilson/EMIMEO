import numpy as np
import matplotlib.pyplot as plt
import os as sys

# from matplotlib import rc

########################################
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




######################################

def gaussian_demo():
    time = np.linspace( -3 , 8 , 500 ) # Create an array of time values between -1 and 5 seconds, with 500 equally-spaced values
    intensity = 2
    init_waist = 1
    
    freq = 2
    ang_freq = 2 * np.pi * freq

    carrier = intensity * np.cos(ang_freq * (time))
    envelope = intensity * np.exp( -( ((time-3)**2) / (2*(init_waist**2)) ))

    gaussian_pulse = intensity * np.exp( -( ((time-3)**2) / (2*(init_waist**2)) )) * np.cos(ang_freq * (time))


    ########## Graphing ##########
    
    fig_gauss_decon, (sub_carrier, sub_envelope) = plt.subplots(1, 2, sharey='row') # Create a figure with three sub-plots side by side (sub_carrier, sub_envelope)

    
    ## Left hand plot
    sub_carrier.plot(time, carrier, color='black')
    # Set limits and aspect ratio
    sub_carrier.set_xlim(-1.5, 6.5)
    sub_carrier.set_aspect(1)
    # Titles and labels
    sub_carrier.set_title(r'Carrier Signal')
    sub_carrier.set_xlabel(r'Time (Seconds)')
    sub_carrier.set_ylabel(r'Intensity (Unitless)')
    
    ## Right hand plot
    sub_envelope.plot(time, envelope)
    sub_envelope.set_xlim(-1.5, 6.5) # Set limits
    sub_envelope.set_aspect(1) # Set aspect ratio
    # Titles and labels
    sub_envelope.set_title(r'|Envelope Signal|')
    sub_envelope.set_xlabel(r'Time (Seconds)')    
    
    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-decon.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    

    plt.figure()

    ## Right hand plot
    plt.plot(time, gaussian_pulse, color='green')
    plt.plot(time, envelope, color='orange', linestyle='dotted')
    plt.plot(time, -envelope, color='orange', linestyle='dotted')

    plt.annotate( 'Envelope signal', xy=(3.75,1.5), xytext=(5, 1.25), arrowprops=dict( arrowstyle="->" ))
    plt.annotate( 'Envelope signal', xy=(3.75,-1.5), xytext=(5, 1.25), arrowprops=dict( arrowstyle="->" ))
    
    plt.xlim(-1.5, 6.5) # Set limits
    # Titles and labels
    plt.title(r'Pulsed Gaussian Signal')
    plt.xlabel(r'Time (Seconds)')
    plt.ylabel(r'Intensity (Unitless)')

    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-demo.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location

    
def main():
    gaussian_demo()
    plt.show()


if __name__ == "__main__":
    main()
    
