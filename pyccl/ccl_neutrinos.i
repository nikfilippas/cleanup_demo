%module ccl_neutrinos

%{
#define SWIG_FILE_WITH_INIT
#include "../include/ccl_neutrinos.h"
%}

// Strip the ccl_ prefix from function names
%rename("%(strip:[ccl_])s") "";

// Automatically document arguments and output types of all functions
%feature("autodoc", "1");

%include "../include/ccl_neutrinos.h"

// Enable vectorised arguments for arrays
%apply (double* IN_ARRAY1, int DIM1) {(double* a, int na), (double* mnu, int nm)};
%apply (double* ARGOUT_ARRAY1, int DIM1){(double* output, int nout)};

%inline %{

void Omeganuh2_vec(double Neff, double TCMB,
                   double* a, int na, 
                   double* mnu, int nm, 
                   double* output, int nout,
                   int* status)
{
    assert(nout == na);
    assert(nm == 3);
    for(int i=0; i < na; i++){
        output[i] = ccl_Omeganuh2(a[i], Neff, mnu, TCMB, NULL, status);
    }   
}

void nu_masses_vec(double OmNuh2, int label, double TCMB,
                          double* output, int nout,
                          int* status)
{
    double* mnu;
    mnu = ccl_nu_masses(OmNuh2, label, TCMB, NULL, status);
    for(int i=0; i < nout; i++){
        output[i] = *(mnu+i);
    }
}

%}
