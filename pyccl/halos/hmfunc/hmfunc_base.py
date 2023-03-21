from ... import ccllib as lib
from ...core import check
from ...background import omega_x
from ...parameters import physical_constants
from ...base import CCLHalosObject, warn_api, deprecated, deprecate_attr
import numpy as np
import functools
from abc import abstractmethod


__all__ = ("MassFunc", "mass_function_from_name",)


class MassFunc(CCLHalosObject):
    """ This class enables the calculation of halo mass functions.
    We currently assume that all mass functions can be written as

    .. math::
        \\frac{dn}{d\\log_{10}M} = f(\\sigma_M)\\,\\frac{\\rho_M}{M}\\,
        \\frac{d\\log \\sigma_M}{d\\log_{10} M}

    where :math:`\\sigma_M^2` is the overdensity variance on spheres with a
    radius given by the Lagrangian radius for mass M.

    * Subclasses implementing analytical mass function parametrizations
      can be created by overriding the ``_get_fsigma`` method.

    * Subclasses may have particular implementations of
      ``_check_mass_def_strict`` to ensure consistency of the halo mass
      definition.

    Args:
        mass_def (:class:`~pyccl.halos.massdef.MassDef`):
            a mass definition object that fixes
            the mass definition used by this mass function
            parametrization.
        mass_def_strict (bool): if False, consistency of the mass
            definition will be ignored.
    """
    __repr_attrs__ = ("mass_def", "mass_def_strict",)
    __getattr__ = deprecate_attr(pairs=[('mdef', 'mass_def')]
                                 )(super.__getattribute__)

    @warn_api
    def __init__(self, *, mass_def=None, mass_def_strict=True):
        self.mass_def_strict = mass_def_strict
        # Check if mass definition was provided and check that it's sensible.
        if mass_def is not None:
            if self._check_mass_def(mass_def):
                raise ValueError("Mass function " + self.name +
                                 " is not compatible with mass definition" +
                                 " Delta = %s, " % (mass_def.Delta) +
                                 " rho = " + mass_def.rho_type)
            self.mass_def = mass_def
        self._setup()

    def _setup(self):
        """ Use this function to initialize any internal attributes
        of this object. This function is called at the very end of the
        constructor call.
        """

    def _check_mass_def_strict(self, mass_def):
        return False

    def _check_mass_def(self, mass_def):
        """ Return False if the input mass definition agrees with
        the definitions for which this mass function parametrization
        works. True otherwise. This function gets called at the
        start of the constructor call.

        Args:
            mass_def (:class:`~pyccl.halos.massdef.MassDef`):
                a mass definition object.

        Returns:
            bool: True if the mass definition is not compatible with \
                this mass function parametrization. False otherwise.
        """
        if self.mass_def_strict:
            return self._check_mass_def_strict(mass_def)
        return False

    def _get_consistent_mass(self, cosmo, M, a, mass_def_other):
        """ Transform a halo mass with a given mass definition into
        the corresponding mass definition that was used to initialize
        this object.

        Args:
            cosmo (:class:`~pyccl.core.Cosmology`): A Cosmology object.
            M (float or array_like): halo mass in units of M_sun.
            a (float): scale factor.
            mass_def_other (:class:`~pyccl.halos.massdef.MassDef`):
                a mass definition object.

        Returns:
            float or array_like: mass according to this object's \
                mass definition.
        """
        if mass_def_other is not None:
            M_use = mass_def_other.translate_mass(
                cosmo, M, a,
                mass_def_other=self.mass_def)
        else:
            M_use = M
        return np.log10(M_use)

    def _get_Delta_m(self, cosmo, a):
        """ For SO-based mass definitions, this returns the corresponding
        value of Delta for a rho_matter-based definition. This is useful
        mostly for the Tinker mass functions, which are defined for any
        SO mass in general, but explicitly only for Delta_matter.
        """
        delta = self.mass_def.get_Delta(cosmo, a)
        if self.mass_def.rho_type == 'matter':
            return delta
        else:
            om_this = omega_x(cosmo, a, self.mass_def.rho_type)
            om_matt = omega_x(cosmo, a, 'matter')
            return delta * om_this / om_matt

    @warn_api(pairs=[("mdef_other", "mass_def_other")])
    def get_mass_function(self, cosmo, M, a, *, mass_def_other=None):
        """ Returns the mass function for input parameters.

        Args:
            cosmo (:class:`~pyccl.core.Cosmology`): A Cosmology object.
            M (float or array_like): halo mass in units of M_sun.
            a (float): scale factor.
            mass_def_other (:class:`~pyccl.halos.massdef.MassDef`):
                the mass definition object that defines M.

        Returns:
            float or array_like: mass function \
                :math:`dn/d\\log_{10}M` in units of Mpc^-3 (comoving).
        """
        # Initialize sigma(M) splines if needed
        cosmo.compute_sigma()

        M_use = np.atleast_1d(M)
        logM = self._get_consistent_mass(cosmo, M_use,
                                         a, mass_def_other)

        # sigma(M)
        status = 0
        sigM, status = lib.sigM_vec(cosmo.cosmo, a, logM,
                                    len(logM), status)
        check(status, cosmo=cosmo)
        # dlogsigma(M)/dlog10(M)
        dlns_dlogM, status = lib.dlnsigM_dlogM_vec(cosmo.cosmo, a, logM,
                                                   len(logM), status)
        check(status, cosmo=cosmo)

        rho = (physical_constants.RHO_CRITICAL *
               cosmo['Omega_m'] * cosmo['h']**2)
        f = self._get_fsigma(cosmo, sigM, a, 2.302585092994046 * logM)
        mf = f * rho * dlns_dlogM / M_use

        if np.ndim(M) == 0:
            mf = mf[0]
        return mf

    @abstractmethod
    def _get_fsigma(self, cosmo, sigM, a, lnM):
        """ Get the :math:`f(\\sigma_M)` function for this mass function
        object (see description of this class for details).

        Args:
            cosmo (:class:`~pyccl.core.Cosmology`): A Cosmology object.
            sigM (float or array_like): standard deviation in the
                overdensity field on the scale of this halo.
            a (float): scale factor.
            lnM (float or array_like): natural logarithm of the
                halo mass in units of M_sun (provided in addition
                to sigM for convenience in some mass function
                parametrizations).

        Returns:
            float or array_like: :math:`f(\\sigma_M)` function.
        """

    @classmethod
    def from_name(cls, name):
        """ Returns mass function subclass from name string

        Args:
            name (string): a mass function name

        Returns:
            MassFunc subclass corresponding to the input name.
        """
        mass_functions = {c.name: c for c in cls.__subclasses__()}
        if name in mass_functions:
            return mass_functions[name]
        else:
            raise ValueError(f"Mass function {name} not implemented.")


@functools.wraps(MassFunc.from_name)
@deprecated(new_function=MassFunc.from_name)
def mass_function_from_name(name):
    return MassFunc.from_name(name)