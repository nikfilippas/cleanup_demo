/** @file */
#ifndef __CCL_CLASS_H_INCLUDED__
#define __CCL_CLASS_H_INCLUDED__

CCL_BEGIN_DECLS

/*
 * Spline the linear power spectrum for mu-Sigma MG cosmologies.
 * @param cosmo Cosmological parameters
 ^ @param psp The linear power spectrum to spline.
 * @param status, integer indicating the status
 */
void ccl_rescale_musigma_s8(ccl_cosmology* cosmo, ccl_f2d_t *psp,
                            int mg_rescale, int* status);

CCL_END_DECLS

#endif
