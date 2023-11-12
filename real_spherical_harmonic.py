from scipy.special import sph_harm
import numpy as np

def real_spherical_harmonic(l, m, theta, phi):
    if m == 0:
        return sph_harm(m, l, theta, phi).real
    elif m > 0:
        return (sph_harm(m, l, theta, phi) + (-1)**m * sph_harm(-m, l, theta, phi)).real / np.sqrt(2)
    else:
        return (sph_harm(-m, l, theta, phi) - (-1)**m * sph_harm(m, l, theta, phi)).imag / np.sqrt(2)
    
def real_spherical_harmonic_array(l, m, theta, phi):
    if m == 0:
        return sph_harm(m, l, theta, phi).real
    elif m > 0:
        return (sph_harm(m, l, theta, phi) + (-1)**m * sph_harm(-m, l, theta, phi)).real / np.sqrt(2)
    else:
        return (sph_harm(-m, l, theta, phi) - (-1)**m * sph_harm(m, l, theta, phi)).imag / np.sqrt(2)
    