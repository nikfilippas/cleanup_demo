%module ccl_core

%{
#define SWIG_FILE_WITH_INIT
#include "../include/ccl_core.h"
%}

// Automatically document arguments and output types of all functions
%feature("autodoc", "1");

// Strip the ccl_ prefix from function names
%rename("%(strip:[ccl_])s") "";

%include "../include/ccl_core.h"

// Enable vectorised arguments for arrays
%apply (double* IN_ARRAY1, int DIM1) {
            (double* zarr, int nz),
            (double* dfarr, int nf),
            (double* M_nu, int n_m)
};

%inline %{
ccl_parameters parameters_create_vec(
                        double Omega_c, double Omega_b, double Omega_k, 
                        double N_nu_rel, double w0, double wa, double h, 
                        double norm_pk, double n_s, 
                        double* zarr, int nz,
                        double* dfarr, int nf, double* M_nu, int n_m, int* status)
{

    assert(nz == nf);
    if (nz == 0){ nz = -1; }
    return ccl_parameters_create(Omega_c, Omega_b, Omega_k, N_nu_rel, n_m, M_nu, 
                                 w0, wa, h, norm_pk, n_s, 
                                 nz, zarr, dfarr, status);
    
    
}

ccl_parameters parameters_create_nuvec(
                        double Omega_c, double Omega_b, double Omega_k, 
                        double N_nu_rel, double w0, double wa, double h, 
                        double norm_pk, double n_s, int n_mg, double *z_mg, double* df_mg, double* M_nu, int n_m, int* status)
{


    return ccl_parameters_create(Omega_c, Omega_b, Omega_k, N_nu_rel, n_m, M_nu, 
                                 w0, wa, h, norm_pk, n_s, 
                                 n_mg, z_mg, df_mg, status);
    
    
}
%}
