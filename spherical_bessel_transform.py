from scipy.integrate import simps
from scipy.special import spherical_jn as jn
import numpy as np

def spherical_bessel_transform(f_in: np.array, r: np.array, l: int, k: float, dr: float) -> np.array:
    """
    Transform a function f(r) to f(k) using spherical bessel functions.
    """
    f_out = np.zeros(len(k), dtype=np.complex128)
    for i in range(len(k)):
        f_out[i] = simps(f_in * jn(l, k[i] * r) * r**2, r, dx=dr)
    return f_out

