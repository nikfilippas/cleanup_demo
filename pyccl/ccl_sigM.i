%module ccl_massfunc

%{
/* put additional #include here */
%}

%include "../include/ccl_massfunc.h"

// Enable vectorised arguments for arrays
%apply (double* IN_ARRAY1, int DIM1) {(double* logM, int nM)};
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int nout, double* output)};

/* The python code here will be executed before all of the functions that
   follow this directive. */
%feature("pythonprepend") %{
    if numpy.shape(logM) != (nout,):
        raise CCLError("Input shape for `halo_mass` must match `(nout,)`!")
%}

%inline %{

void sigM_vec(ccl_cosmology * cosmo, double a,
	      double *logM, int nM,
	      int nout, double* output, int *status)
{
  // Check if sigma has already been calculated
  if (!cosmo->computed_sigma) {
    *status = CCL_ERROR_SIGMA_INIT;
    ccl_cosmology_set_status_message(
      cosmo,
      "ccl_massfunc.c: ccl_sigmaM(): linear power spctrum has not been computed!");
    return;
  }

  for(int i=0; i<nM; i++)
    output[i] = ccl_sigmaM(cosmo, logM[i], a, status);
}

void dlnsigM_dlogM_vec(ccl_cosmology * cosmo,
		       double *logM, int nM,
		       int nout, double* output, int *status)
{
  // Check if sigma has already been calculated
  if (!cosmo->computed_sigma) {
    *status = CCL_ERROR_SIGMA_INIT;
    ccl_cosmology_set_status_message(
      cosmo,
      "ccl_massfunc.c: ccl_sigmaM(): linear power spctrum has not been computed!");
    return;
  }

  for(int i=0; i<nM; i++)
    output[i] = ccl_dlnsigM_dlogM(cosmo, logM[i], status);
}

%}

/* The directive gets carried between files, so we reset it at the end. */
%feature("pythonprepend") %{ %}
