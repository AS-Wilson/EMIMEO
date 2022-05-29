# To use this module just type include at start of your file:
# from OCNToolbox import *

from numpy import *
from sympy import *


class OCNToolbox:
    """Something - TODO Later"""

    @staticmethod
    def gaussian(time, FWHM, t_0, intensity=1, t_shift=0):
        """Public Method for producing a Gaussian Signal after providing various parameters."""

        try:
            t_0 = FWHM / (2 * sqrt(2 * log(2)))
            gaussian = intensity * exp(- ((time - t_shift)**2 / 2 * t_0**2))
            return gaussian

        except NameError:
            gaussian = intensity * exp(- ((time - t_shift) ** 2 / 2 * t_0 ** 2))
            return gaussian


    @staticmethod
    def split_step_fourier_method(input_spectrum, ang_freq, z_step_size, beta_2=0, gamma=0):
        """Something - TODO Later"""

        # Since gamma is 0 in this code, this will be an array with no value, but the code is still useful
        non_linear_term = gamma * z_step_size

        ## Produce all the dispersive terms first as they are frequency dependent
        for i in range(len(ang_freq)):
            try:
                dispersive_terms = np.concatenate(dispersive_terms, np.exp(complex(0, ((ang_freq[i]**2 / 2) * beta_2 * z_step_size))).T)

            except NameError:
                dispersive_terms = np.array(np.exp(complex(0, ((ang_freq[i]**2 / 2) * beta_2 * z_step_size))))


        output_spectrum = np.r_[complex(input_spectrum[0], 0)]
        output_time = np.r_[ifft(ifftshift(input_spectrum[0])) * ang_freq_interval / np.sqrt(2.0*np.pi) * len(input_spectrum)]

            envelope_complex = [complex(envelope[0], 0)]

            for i in range(1, no_of_samples):
                envelope_complex.append(complex(envelope[i], 0))

            pulse_output_time = np.r_[envelope_complex]
            pulse_output_spectrum = np.r_[envelope_raw_spectrum]

            # Create vectors to manipulate when carrying out our steps
            current_step_spectrum = envelope_raw_spectrum

            # The actual Split-Step Fourier Method Calculations
            for x in range(1, no_of_samples):

                for j in range(no_of_samples):
                    # Dispersive (Linear) Step
                    current_step_spectrum[j] = current_step_spectrum[j] * dispersive_terms[j]

                # Prep for Non-Linear Step
                current_step_time = ifft(ifftshift(current_step_spectrum)) * ang_freq_interval / np.sqrt(2.0 * np.pi) * no_of_samples

                for j in range(no_of_samples):
                    # Non-Linear Step
                    current_step_time[j] = current_step_time[j] * np.exp(complex(0, (non_linear_term *
                                                                                     (abs(current_step_time[j]) ** 2))))

                pulse_output_time = np.row_stack((pulse_output_time, current_step_time))
                pulse_output_spectrum = np.row_stack((pulse_output_spectrum, current_step_spectrum))

