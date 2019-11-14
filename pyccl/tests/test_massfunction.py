import numpy as np
import pytest
import pyccl as ccl


COSMO = ccl.Cosmology(
    Omega_c=0.27, Omega_b=0.045, h=0.67, sigma8=0.8, n_s=0.96,
    transfer_function='bbks', matter_power_spectrum='halofit')
MS = [1E13, [1E12, 1E15], np.array([1E12, 1E15])]
MF_EQUIV = {'tinker10': 'Tinker10',
            'tinker': 'Tinker08',
            'watson': 'Watson13',
            'shethtormen': 'Sheth99',
            'angulo': 'Angulo12'}
MF_TYPES = sorted(list(MF_EQUIV.keys()))


@pytest.mark.parametrize('mf_type', MF_TYPES)
def test_massfunc_smoke(mf_type):
    cosmo = ccl.Cosmology(
        Omega_c=0.27, Omega_b=0.045, h=0.67, sigma8=0.8, n_s=0.96,
        transfer_function='bbks', matter_power_spectrum='linear',
        mass_function=mf_type)
    hmf_cls = ccl.halos.mass_function_from_name(MF_EQUIV[mf_type])
    hmf = hmf_cls(cosmo)
    for m in MS:
        nm_old = ccl.massfunc(cosmo, m, 1.)
        nm_new = hmf.get_mass_function(cosmo, m, 1.)
        assert np.all(np.isfinite(nm_old))
        assert np.shape(nm_old) == np.shape(m)
        assert np.all(np.array(nm_old) ==
                      np.array(nm_new))


def test_halo_bias_smoke():
    cosmo = ccl.Cosmology(
        Omega_c=0.27, Omega_b=0.045, h=0.67, sigma8=0.8, n_s=0.96,
        transfer_function='bbks', matter_power_spectrum='linear',
        mass_function='tinker10')
    hmf = ccl.halos.HaloBiasTinker10(cosmo)
    for m in MS:
        bm_old = ccl.halo_bias(cosmo, m, 1.)
        bm_new = hmf.get_halo_bias(cosmo, m, 1.)
        assert np.all(np.isfinite(bm_old))
        assert np.shape(bm_old) == np.shape(m)
        assert np.all(np.array(bm_old) ==
                      np.array(bm_new))


def test_m2r_smoke():
    for m in MS:
        r_old = ccl.massfunc_m2r(COSMO, m)
        r_new = ccl.halos.mass2radius_lagrangian(COSMO, m)
        assert np.all(np.isfinite(r_old))
        assert np.shape(r_old) == np.shape(m)
        assert np.all(np.array(r_old) ==
                      np.array(r_new))


@pytest.mark.parametrize('m', [
    1e14,
    int(1e14),
    [1e14, 1e15],
    np.array([1e14, 1e15])])
def test_massfunc_smoke(m):
    a = 0.8
    mf = ccl.massfunc(COSMO, m, a)
    assert np.all(np.isfinite(mf))
    assert np.shape(mf) == np.shape(m)


@pytest.mark.parametrize('m', [
    1e14,
    int(1e14),
    [1e14, 1e15],
    np.array([1e14, 1e15])])
def test_massfunc_m2r_smoke(m):
    r = ccl.massfunc_m2r(COSMO, m)
    assert np.all(np.isfinite(r))
    assert np.shape(r) == np.shape(m)


@pytest.mark.parametrize('m', [
    1e14,
    int(1e14),
    [1e14, 1e15],
    np.array([1e14, 1e15])])
def test_sigmaM_smoke(m):
    a = 0.8
    s = ccl.sigmaM(COSMO, m, a)
    assert np.all(np.isfinite(s))
    assert np.shape(s) == np.shape(m)


@pytest.mark.parametrize('m', [
    1e14,
    int(1e14),
    [1e14, 1e15],
    np.array([1e14, 1e15])])
def test_halo_bias_smoke(m):
    a = 0.8
    b = ccl.halo_bias(COSMO, m, a)
    assert np.all(np.isfinite(b))
    assert np.shape(b) == np.shape(m)
