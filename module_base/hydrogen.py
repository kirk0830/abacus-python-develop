from scipy.special import assoc_laguerre as laguerre
import numpy as np

def quasi_hydrogen_orbitals_radial(n: int, l: int, r: np.array, Z: float) -> np.array:
    """
    Generate the radial part of the quasi-hydrogen orbitals.
    """
    rho = 2 * Z * r / n
    rho0 = 2 * Z / n
    rho1 = 2 * Z / (n - 1)
    rho2 = 2 * Z / (n - 2)

    if l == 0:
        return np.sqrt(rho0**3) * np.exp(-rho / 2) * laguerre(n - 1, 0, rho)
    elif l == 1:
        return np.sqrt(rho1**3) * np.exp(-rho / 2) * laguerre(n - 2, 1, rho)
    elif l == 2:
        return np.sqrt(rho2**3) * np.exp(-rho / 2) * laguerre(n - 3, 2, rho)
    else:
        raise ValueError('l must be 0, 1, or 2')
