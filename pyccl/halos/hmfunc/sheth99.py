from ... import ccllib as lib
from ...base import warn_api
from ...pyutils import check
from ..massdef import MassDef
from ..halo_model_base import MassFunc
import numpy as np


__all__ = ("MassFuncSheth99",)


class MassFuncSheth99(MassFunc):
    """ Implements mass function described in arXiv:astro-ph/9901122
    This parametrization is only valid for 'fof' masses.

    Args:
        mass_def (:class:`~pyccl.halos.massdef.MassDef`):
            a mass definition object.
            this parametrization accepts FoF masses only.
            The default is 'fof'.
        mass_def_strict (bool): if False, consistency of the mass
            definition will be ignored.
        use_delta_c_fit (bool): if True, use delta_crit given by
            the fit of Nakamura & Suto 1997. Otherwise use
            delta_crit = 1.68647.
    """
    __repr_attrs__ = ("mass_def", "mass_def_strict", "use_delta_c_fit",)
    name = 'Sheth99'

    @warn_api
    def __init__(self, *,
                 mass_def=MassDef('fof', 'matter'),
                 mass_def_strict=True,
                 use_delta_c_fit=False):
        self.use_delta_c_fit = use_delta_c_fit
        super().__init__(mass_def=mass_def, mass_def_strict=mass_def_strict)

    def _setup(self):
        self.A = 0.21615998645
        self.p = 0.3
        self.a = 0.707

    def _check_mass_def_strict(self, mass_def):
        if mass_def.Delta != 'fof':
            return True

    def _get_fsigma(self, cosmo, sigM, a, lnM):
        if self.use_delta_c_fit:
            status = 0
            delta_c, status = lib.dc_NakamuraSuto(cosmo.cosmo, a, status)
            check(status, cosmo=cosmo)
        else:
            delta_c = 1.68647

        nu = delta_c / sigM
        return nu * self.A * (1. + (self.a * nu**2)**(-self.p)) * \
            np.exp(-self.a * nu**2/2.)
