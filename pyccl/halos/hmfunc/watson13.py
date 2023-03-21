from ...base import warn_api
from ..massdef import MassDef200m
from .hmfunc_base import MassFunc
import numpy as np


__all__ = ("MassFuncWatson13",)


class MassFuncWatson13(MassFunc):
    """ Implements mass function described in arXiv:1212.0095.

    Args:
        mass_def (:class:`~pyccl.halos.massdef.MassDef`):
            a mass definition object.
            this parametrization accepts fof and any SO masses.
            The default is '200m'.
            If `None`, Delta = 200 (matter) will be used.
        mass_def_strict (bool): if False, consistency of the mass
            definition will be ignored.
    """
    name = 'Watson13'

    @warn_api
    def __init__(self, *,
                 mass_def=MassDef200m(),
                 mass_def_strict=True):
        super().__init__(mass_def=mass_def, mass_def_strict=mass_def_strict)

    def _setup(self):
        self.is_fof = self.mass_def.Delta == 'fof'

    def _check_mass_def_strict(self, mass_def):
        if mass_def.Delta == 'vir':
            return True
        return False

    def _get_fsigma(self, cosmo, sigM, a, lnM):
        if self.is_fof:
            pA = 0.282
            pa = 2.163
            pb = 1.406
            pc = 1.210
            return pA * ((pb / sigM)**pa + 1.) * np.exp(-pc / sigM**2)
        else:
            om = cosmo.omega_x(a, "matter")
            Delta_178 = self.mass_def.Delta / 178.0

            if a == 1.0:
                pA = 0.194
                pa = 1.805
                pb = 2.267
                pc = 1.287
            elif a < 0.14285714285714285:  # z>6
                pA = 0.563
                pa = 3.810
                pb = 0.874
                pc = 1.453
            else:
                pA = om * (1.097 * a**3.216 + 0.074)
                pa = om * (5.907 * a**3.058 + 2.349)
                pb = om * (3.136 * a**3.599 + 2.344)
                pc = 1.318

            f_178 = pA * ((pb / sigM)**pa + 1.) * np.exp(-pc / sigM**2)
            C = np.exp(0.023 * (Delta_178 - 1.0))
            d = -0.456 * om - 0.139
            Gamma = (C * Delta_178**d *
                     np.exp(0.072 * (1.0 - Delta_178) / sigM**2.130))
            return f_178 * Gamma