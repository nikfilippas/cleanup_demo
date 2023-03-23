from ..base import UnlockInstance, warn_api
from ..pk2d import Pk2D
from .profiles_2pt import Profile2pt
import numpy as np


__all__ = ("halomod_power_spectrum", "halomod_Pk2D",)


@warn_api(pairs=[("supress_1h", "suppress_1h")],
          reorder=["prof_2pt", "prof2", "p_of_k_a", "normprof1", "normprof2"])
def halomod_power_spectrum(cosmo, hmc, k, a, prof, *,
                           prof2=None, prof_2pt=None,
                           normprof1=None, normprof2=None,
                           p_of_k_a=None,
                           get_1h=True, get_2h=True,
                           smooth_transition=None, suppress_1h=None):
    """ Computes the halo model power spectrum for two
    quantities defined by their respective halo profiles.
    The halo model power spectrum for two profiles
    :math:`u` and :math:`v` is:

    .. math::
        P_{u,v}(k,a) = I^0_2(k,a|u,v) +
        I^1_1(k,a|u)\\,I^1_1(k,a|v)\\,P_{\\rm lin}(k,a)

    where :math:`P_{\\rm lin}(k,a)` is the linear matter
    power spectrum, :math:`I^1_1` is defined in the documentation
    of :meth:`~HMCalculator.I_1_1`, and :math:`I^0_2` is defined
    in the documentation of :meth:`~HMCalculator.I_0_2`.

    Args:
        cosmo (:class:`~pyccl.core.Cosmology`): a Cosmology object.
        hmc (:class:`HMCalculator`): a halo model calculator.
        k (float or array_like): comoving wavenumber in Mpc^-1.
        a (float or array_like): scale factor.
        prof (:class:`~pyccl.halos.profiles.HaloProfile`): halo
            profile.
        prof2 (:class:`~pyccl.halos.profiles.HaloProfile`): a
            second halo profile. If `None`, `prof` will be used as
            `prof2`.
        normprof1 (bool): (Deprecated - do not use)
            if `True`, this integral will be
            normalized by :math:`I^0_1(k\\rightarrow 0,a|u)`
            (see :meth:`~HMCalculator.I_0_1`), where
            :math:`u` is the profile represented by `prof`.
        normprof2 (bool): (Deprecated - do not use)
            if `True`, this integral will be
            normalized by :math:`I^0_1(k\\rightarrow 0,a|v)`
            (see :meth:`~HMCalculator.I_0_1`), where
            :math:`v` is the profile represented by `prof2`.
        prof_2pt (:class:`~pyccl.halos.profiles_2pt.Profile2pt`):
            a profile covariance object
            returning the the two-point moment of the two profiles
            being correlated. If `None`, the default second moment
            will be used, corresponding to the products of the means
            of both profiles.
        p_of_k_a (:class:`~pyccl.pk2d.Pk2D`): a `Pk2D` object to
            be used as the linear matter power spectrum. If `None`,
            the power spectrum stored within `cosmo` will be used.
        get_1h (bool): if `False`, the 1-halo term (i.e. the first
            term in the first equation above) won't be computed.
        get_2h (bool): if `False`, the 2-halo term (i.e. the second
            term in the first equation above) won't be computed.
        smooth_transition (function or None):
            Modify the halo model 1-halo/2-halo transition region
            via a time-dependent function :math:`\\alpha(a)`,
            defined as in HMCODE-2020 (``arXiv:2009.01858``): :math:`P(k,a)=
            (P_{1h}^{\\alpha(a)}(k)+P_{2h}^{\\alpha(a)}(k))^{1/\\alpha}`.
            If `None` the extra factor is just 1.
        suppress_1h (function or None):
            Suppress the 1-halo large scale contribution by a
            time- and scale-dependent function :math:`k_*(a)`,
            defined as in HMCODE-2020 (``arXiv:2009.01858``):
            :math:`\\frac{(k/k_*(a))^4}{1+(k/k_*(a))^4}`.
            If `None` the standard 1-halo term is returned with no damping.

    Returns:
        float or array_like: integral values evaluated at each
        combination of `k` and `a`. The shape of the output will
        be `(N_a, N_k)` where `N_k` and `N_a` are the sizes of
        `k` and `a` respectively. If `k` or `a` are scalars, the
        corresponding dimension will be squeezed out on output.
    """
    a_use = np.atleast_1d(a).astype(float)
    k_use = np.atleast_1d(k).astype(float)

    # Check inputs
    if smooth_transition is not None:
        if not (get_1h and get_2h):
            raise ValueError("transition region can only be modified "
                             "when both 1-halo and 2-halo terms are queried")
    if suppress_1h is not None:
        if not get_1h:
            raise ValueError("can't suppress the 1-halo term "
                             "when get_1h is False")

    if prof2 is None:
        prof2 = prof
    if prof_2pt is None:
        prof_2pt = Profile2pt()

    # TODO: Remove for CCLv3.
    if normprof1 is not None:
        with UnlockInstance(prof):
            prof.normprof = normprof1
    if normprof2 is not None:
        with UnlockInstance(prof2):
            prof2.normprof = normprof2

    # Power spectrum
    if isinstance(p_of_k_a, Pk2D):
        def pkf(sf):
            return p_of_k_a.eval(k_use, sf, cosmo)
    elif (p_of_k_a is None) or (str(p_of_k_a) == 'linear'):
        def pkf(sf):
            return cosmo.linear_matter_power(k_use, sf)
    elif str(p_of_k_a) == 'nonlinear':
        def pkf(sf):
            return cosmo.nonlin_matter_power(k_use, sf)
    else:
        raise TypeError("p_of_k_a must be `None`, \'linear\', "
                        "\'nonlinear\' or a `Pk2D` object")

    na = len(a_use)
    nk = len(k_use)
    out = np.zeros([na, nk])
    for ia, aa in enumerate(a_use):
        # Compute first profile normalization
        if prof.normprof:
            norm1 = hmc.profile_norm(cosmo, aa, prof)
        else:
            norm1 = 1
        # Compute second profile normalization
        if prof2 == prof:
            norm2 = norm1
        else:
            if prof2.normprof:
                norm2 = hmc.profile_norm(cosmo, aa, prof2)
            else:
                norm2 = 1
        norm = norm1 * norm2

        if get_2h:
            # Compute first bias factor
            i11_1 = hmc.I_1_1(cosmo, k_use, aa, prof)

            # Compute second bias factor
            if prof2 == prof:
                i11_2 = i11_1
            else:
                i11_2 = hmc.I_1_1(cosmo, k_use, aa, prof2)

            # Compute 2-halo power spectrum
            pk_2h = pkf(aa) * i11_1 * i11_2
        else:
            pk_2h = 0.

        if get_1h:
            pk_1h = hmc.I_0_2(cosmo, k_use, aa, prof,
                              prof2=prof2, prof_2pt=prof_2pt)
            if suppress_1h is not None:
                ks = suppress_1h(aa)
                pk_1h *= (k_use / ks)**4 / (1 + (k_use / ks)**4)
        else:
            pk_1h = 0.

        # Transition region
        if smooth_transition is None:
            out[ia, :] = (pk_1h + pk_2h) * norm
        else:
            alpha = smooth_transition(aa)
            out[ia, :] = (pk_1h**alpha + pk_2h**alpha)**(1/alpha) * norm

    if np.ndim(a) == 0:
        out = np.squeeze(out, axis=0)
    if np.ndim(k) == 0:
        out = np.squeeze(out, axis=-1)
    return out


@warn_api(pairs=[("supress_1h", "suppress_1h")],
          reorder=["prof_2pt", "prof2", "p_of_k_a", "normprof1", "normprof2"])
def halomod_Pk2D(cosmo, hmc, prof, *,
                 prof2=None, prof_2pt=None,
                 normprof1=None, normprof2=None,
                 p_of_k_a=None,
                 get_1h=True, get_2h=True,
                 lk_arr=None, a_arr=None,
                 extrap_order_lok=1, extrap_order_hik=2,
                 smooth_transition=None, suppress_1h=None):
    """ Returns a :class:`~pyccl.pk2d.Pk2D` object containing
    the halo-model power spectrum for two quantities defined by
    their respective halo profiles. See :meth:`halomod_power_spectrum`
    for more details about the actual calculation.

    Args:
        cosmo (:class:`~pyccl.core.Cosmology`): a Cosmology object.
        hmc (:class:`HMCalculator`): a halo model calculator.
        prof (:class:`~pyccl.halos.profiles.HaloProfile`): halo
            profile.
        prof2 (:class:`~pyccl.halos.profiles.HaloProfile`): a
            second halo profile. If `None`, `prof` will be used as
            `prof2`.
        prof_2pt (:class:`~pyccl.halos.profiles_2pt.Profile2pt`):
            a profile covariance object
            returning the the two-point moment of the two profiles
            being correlated. If `None`, the default second moment
            will be used, corresponding to the products of the means
            of both profiles.
        normprof1 (bool): if `True`, this integral will be
            normalized by :math:`I^0_1(k\\rightarrow 0,a|u)`
            (see :meth:`~HMCalculator.I_0_1`), where
            :math:`u` is the profile represented by `prof`.
        normprof2 (bool): if `True`, this integral will be
            normalized by :math:`I^0_1(k\\rightarrow 0,a|v)`
            (see :meth:`~HMCalculator.I_0_1`), where
            :math:`v` is the profile represented by `prof2`.
        p_of_k_a (:class:`~pyccl.pk2d.Pk2D`): a `Pk2D` object to
            be used as the linear matter power spectrum. If `None`,
            the power spectrum stored within `cosmo` will be used.
        get_1h (bool): if `False`, the 1-halo term (i.e. the first
            term in the first equation above) won't be computed.
        get_2h (bool): if `False`, the 2-halo term (i.e. the second
            term in the first equation above) won't be computed.
        a_arr (array): an array holding values of the scale factor
            at which the halo model power spectrum should be
            calculated for interpolation. If `None`, the internal
            values used by `cosmo` will be used.
        lk_arr (array): an array holding values of the natural
            logarithm of the wavenumber (in units of Mpc^-1) at
            which the halo model power spectrum should be calculated
            for interpolation. If `None`, the internal values used
            by `cosmo` will be used.
        extrap_order_lok (int): extrapolation order to be used on
            k-values below the minimum of the splines. See
            :class:`~pyccl.pk2d.Pk2D`.
        extrap_order_hik (int): extrapolation order to be used on
            k-values above the maximum of the splines. See
            :class:`~pyccl.pk2d.Pk2D`.
        smooth_transition (function or None):
            Modify the halo model 1-halo/2-halo transition region
            via a time-dependent function :math:`\\alpha(a)`,
            defined as in HMCODE-2020 (``arXiv:2009.01858``): :math:`P(k,a)=
            (P_{1h}^{\\alpha(a)}(k)+P_{2h}^{\\alpha(a)}(k))^{1/\\alpha}`.
            If `None` the extra factor is just 1.
        suppress_1h (function or None):
            Suppress the 1-halo large scale contribution by a
            time- and scale-dependent function :math:`k_*(a)`,
            defined as in HMCODE-2020 (``arXiv:2009.01858``):
            :math:`\\frac{(k/k_*(a))^4}{1+(k/k_*(a))^4}`.
            If `None` the standard 1-halo term is returned with no damping.

    Returns:
        :class:`~pyccl.pk2d.Pk2D`: halo model power spectrum.
    """
    if lk_arr is None:
        lk_arr = cosmo.get_pk_spline_lk()
    if a_arr is None:
        a_arr = cosmo.get_pk_spline_a()

    pk_arr = halomod_power_spectrum(
        cosmo, hmc, np.exp(lk_arr), a_arr,
        prof, prof2=prof2, prof_2pt=prof_2pt, p_of_k_a=p_of_k_a,
        normprof1=normprof1, normprof2=normprof2,  # TODO: remove for CCLv3
        get_1h=get_1h, get_2h=get_2h,
        smooth_transition=smooth_transition, suppress_1h=suppress_1h)

    return Pk2D(a_arr=a_arr, lk_arr=lk_arr, pk_arr=pk_arr,
                extrap_order_lok=extrap_order_lok,
                extrap_order_hik=extrap_order_hik,
                is_logp=False)
