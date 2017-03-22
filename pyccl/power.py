
import ccllib as lib
from pyutils import _vectorize_fn, _vectorize_fn2

def linear_matter_power(cosmo, a, k):
    """The linear matter power spectrum; Mpc^-3.

    Args:
        cosmo (:obj:`ccl.cosmology`): Cosmological parameters.
        a (float or array_like): Scale factor.
        k (float): Wavenumber; Mpc^-1.

    Returns:
        linear_matter_power (float or array_like): Linear matter power spectrum; Mpc^-3.

    """
    return _vectorize_fn2(lib.linear_matter_power, 
                          lib.linear_matter_power_vec, cosmo, k, a)

def nonlin_matter_power(cosmo, a, k):
    """The nonlinear matter power spectrum; Mpc^-3.

    Args:
        cosmo (:obj:`ccl.cosmology`): Cosmological parameters.
        a (float or array_like): Scale factor.
        k (float): Wavenumber; Mpc^-1.

    Returns:
        nonlin_matter_power (float or array_like): Nonlinear matter power spectrum; Mpc^-3.

    """
    return _vectorize_fn2(lib.nonlin_matter_power, 
                          lib.nonlin_matter_power_vec, cosmo, k, a)

def sigmaR(cosmo, R):
    """RMS variance in a top-hat sphere of radius R.

    Args:
        cosmo (:obj:`ccl.cosmology`): Cosmological parameters.
        R (float or array_like): Radius; Mpc.

    Returns:
        sigmaR (float or array_like): RMS variance in top-hat sphere.

    """
    return _vectorize_fn(lib.sigmaR, 
                         lib.sigmaR_vec, cosmo, R)

def sigma8(cosmo):
    """RMS variance in a top-hat sphere of radius 8 Mpc.

    Args:
        cosmo (:obj:`ccl.cosmology`): Cosmological parameters.

    Returns:
        sigma8 (float): RMS variance in top-hat sphere of radius 8 Mpc.

    """
    print cosmo['h']
    return sigmaR(cosmo,8./cosmo['h'])

