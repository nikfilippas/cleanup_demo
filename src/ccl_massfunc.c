#include "ccl.h"
#include "ccl_core.h"
#include "ccl_utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "ccl_power.h"
#include "ccl_massfunc.h"
#include "ccl_error.h"

void ccl_cosmology_compute_hmfparams(ccl_cosmology * cosmo, int *status)
{
  if(cosmo->computed_hmfparams)
    return;

  // declare parameter splines on case-by-case basis
  switch(cosmo->config.mass_function_method){
  case ccl_tinker10:
    printf("");
    double delta[9] = {200.0, 300.0, 400.0, 600.0, 800.0, 1200.0, 1600.0, 2400.0, 3200.0};
    double lgdelta[9];
    double alpha[9] = {0.368, 0.363, 0.385, 0.389, 0.393, 0.365, 0.379, 0.355, 0.327};
    double beta[9] = {0.589, 0.585, 0.544, 0.543, 0.564, 0.623, 0.637, 0.673, 0.702};
    double gamma[9] ={0.864, 0.922, 0.987, 1.09, 1.20, 1.34, 1.50, 1.68, 1.81};
    double phi[9] = {-0.729, -0.789, -0.910, -1.05, -1.20, -1.26, -1.45, -1.50, -1.49};
    double eta[9] = {-0.243, -0.261, -0.261, -0.273, -0.278, -0.301, -0.301, -0.319, -0.336};
    int nd = 9;
    int i;

    for(i=0; i<nd; i++){
      lgdelta[i] = log10(delta[i]);
    }

    gsl_spline * alphahmf = gsl_spline_alloc(D_SPLINE_TYPE, nd);
    *status = gsl_spline_init(alphahmf, lgdelta, alpha, nd);
    if (*status){
      gsl_spline_free(alphahmf);
      *status = CCL_ERROR_SPLINE ;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_hmfparams(): Error creating alpha(D) spline\n");
      return;
    }

    gsl_spline * betahmf  = gsl_spline_alloc(D_SPLINE_TYPE, nd);
    *status = gsl_spline_init(betahmf, lgdelta, beta, nd);
    if (*status){
      gsl_spline_free(alphahmf);
      gsl_spline_free(betahmf);
      *status = CCL_ERROR_SPLINE ;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_hmfparams(): Error creating beta(D) spline\n");
      return;
    }

    gsl_spline * gammahmf = gsl_spline_alloc(D_SPLINE_TYPE, nd);
    *status = gsl_spline_init(gammahmf, lgdelta, gamma, nd);
    if (*status){
      gsl_spline_free(alphahmf);
      gsl_spline_free(betahmf);
      gsl_spline_free(gammahmf);
      *status = CCL_ERROR_SPLINE ;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_hmfparams(): Error creating gamma(D) spline\n");
      return;
    }

    gsl_spline * phihmf   = gsl_spline_alloc(D_SPLINE_TYPE, nd);
    *status = gsl_spline_init(phihmf, lgdelta, phi, nd);
    if (*status){
      gsl_spline_free(alphahmf);
      gsl_spline_free(betahmf);
      gsl_spline_free(gammahmf);
      gsl_spline_free(phihmf);
      *status = CCL_ERROR_SPLINE ;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_hmfparams(): Error creating phi(D) spline\n");
      return;
    }

    gsl_spline * etahmf   = gsl_spline_alloc(D_SPLINE_TYPE, nd);
    *status = gsl_spline_init(etahmf, lgdelta, eta, nd);
    if (*status){
      gsl_spline_free(alphahmf);
      gsl_spline_free(betahmf);
      gsl_spline_free(gammahmf);
      gsl_spline_free(phihmf);
      gsl_spline_free(etahmf);
      *status = CCL_ERROR_SPLINE ;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_hmfparams(): Error creating eta(D) spline\n");
      return;
    }
    if(cosmo->data.accelerator_d==NULL)
      cosmo->data.accelerator_d=gsl_interp_accel_alloc();
    cosmo->data.alphahmf = alphahmf;
    cosmo->data.betahmf = betahmf;
    cosmo->data.gammahmf = gammahmf;
    cosmo->data.phihmf = phihmf;
    cosmo->data.etahmf = etahmf;
    cosmo->computed_hmfparams = true;

    break;
  default:
  //TODO: Error message goes here.
    break;
  }
}

/*----- ROUTINE: ccl_massfunc_f -----
INPUT: cosmology+parameters, a halo mass, and scale factor
TASK: Outputs fitting function for use in halo mass function calculation;
  currently only supports:
    ccl_tinker (arxiv 0803.2706 )
    ccl_tinker10 (arxiv 1001.3162 )
    ccl_angulo (arxiv 1203.3216 ) 
    ccl_watson (arxiv 1212.0095 )
*/

static double massfunc_f(ccl_cosmology *cosmo, double halomass, double a, double odelta, int *status)
{
  double fit_A, fit_a, fit_b, fit_c, fit_d, overdensity_delta;
  double Omega_m_a;
  double delta_c_Tinker, nu;

  double sigma=ccl_sigmaM(cosmo, halomass, a, status);

  switch(cosmo->config.mass_function_method){
  case ccl_tinker:
    if (odelta == 200.0){
      fit_A = 0.186;
      fit_a = 1.47;
      fit_b = 2.57; 
      fit_c = 1.19;
    }
    else if (odelta == 300.0){
      fit_A = 0.200;
      fit_a = 1.52;
      fit_b = 2.25;
      fit_c = 1.27;
    }
    else if (odelta == 400.0){
      fit_A = 0.212;
      fit_a = 1.56;
      fit_b = 2.05;
      fit_c = 1.34;
    }
    else if (odelta == 600.0){
      fit_A = 0.218;
      fit_a = 1.61;
      fit_b = 1.87;
      fit_c = 1.45;
    }
    else if (odelta == 800.0){
      fit_A = 0.248;
      fit_a = 1.87;
      fit_b = 1.59;
      fit_c = 1.58;
    }
    else if (odelta == 1200.0){
      fit_A = 0.255;
      fit_a = 2.13;
      fit_b = 1.51;
      fit_c = 1.80;
    }
    else if (odelta == 1600.0){
      fit_A = 0.260;
      fit_a = 2.30;
      fit_b = 1.46;
      fit_c = 1.97;
    }
    else if (odelta == 2400.0){
      fit_A = 0.260;
      fit_a = 2.53;
      fit_b = 1.44;
      fit_c = 2.24;
    }
    else if (odelta == 3200.0){
      fit_A = 0.260;
      fit_a = 2.66;
      fit_b = 1.41;
      fit_c = 2.44;
    }
    else{
      *status = CCL_ERROR_HMF_INTERP;
      strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_massfunc_f(): Tinker 2008 only supported for Delta = {200, 300, 400, 600, 800, 1200, 1600, 2400, and 3200. Calculation continues assuming Delta =200.\n");
      odelta = 200;
      fit_A = 0.186;
      fit_a = 1.47;
      fit_b = 2.57; 
      fit_c = 1.19;
    }

    fit_a = fit_a*pow(a, 0.06);
    fit_d = pow(10, -1.0*pow(0.75 / log10(odelta / 75.0), 1.2));
    fit_b = fit_b*pow(a, fit_d);

    return fit_A*(pow(sigma/fit_b,-fit_a)+1.0)*exp(-fit_c/sigma/sigma);
    break;
    //this version uses f(nu) parameterization from Eq. 8 in Tinker et al. 2010
    // use this for consistency with Tinker et al. 2010 fitting function for halo bias
  case ccl_tinker10:
    if (!cosmo->computed_hmfparams){
        ccl_cosmology_compute_hmfparams(cosmo, status);
        ccl_check_status(cosmo, status);
    }
    //critical collapse overdensity assumed in this model
    delta_c_Tinker = 1.686;
    nu = delta_c_Tinker/(sigma);

    fit_A = gsl_spline_eval(cosmo->data.alphahmf, log10(odelta), cosmo->data.accelerator_d); //alpha in Eq. 8
    fit_a = gsl_spline_eval(cosmo->data.etahmf, log10(odelta), cosmo->data.accelerator_d)*pow(a, -0.27); //eta in Eq. 8
    fit_b = gsl_spline_eval(cosmo->data.betahmf, log10(odelta), cosmo->data.accelerator_d)*pow(a, -0.20); //beta in Eq. 8
    fit_c = gsl_spline_eval(cosmo->data.gammahmf, log10(odelta), cosmo->data.accelerator_d)*pow(a, 0.01); //gamma in Eq. 8
    fit_d = gsl_spline_eval(cosmo->data.phihmf, log10(odelta), cosmo->data.accelerator_d)*pow(a, 0.08); //phi in Eq. 8;

    printf("%le %le %le %le %le\n", fit_A, fit_a, fit_b, fit_c, fit_d);

    return nu*fit_A*(1.+pow(fit_b*nu,-2.*fit_d))*pow(nu, 2.*fit_a)*exp(-0.5*fit_c*nu*nu);
    break;

  case ccl_watson:
    *status = CCL_ERROR_HMF_INTERP;
    strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_massfunc_f(): Watson HMF only supported for Delta = 200. Calculation continues assuming Delta = 200.\n");
    Omega_m_a = ccl_omega_x(cosmo, a, ccl_omega_m_label);
    fit_A = Omega_m_a*(0.990*pow(a,3.216)+0.074);
    fit_a = Omega_m_a*(5.907*pow(a,3.599)+2.344);
    fit_b = Omega_m_a*(3.136*pow(a,3.058)+2.349);
    fit_c = 1.318;

    return fit_A*(pow(sigma/fit_b,-fit_a)+1.0)*exp(-fit_c/sigma/sigma);

  case ccl_angulo:
    *status = CCL_ERROR_HMF_INTERP;
    strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_massfunc_f(): Angulo HMF only supported for Delta = 200. Calculation continues assuming Delta = 200.\n");
    fit_A = 0.201;
    fit_a = 2.08;
    fit_b = 1.7;
    fit_c = 1.172;

    return fit_A*pow( (fit_a/sigma)+1.0, fit_b)*exp(-fit_c/sigma/sigma);

  default:
    *status = CCL_ERROR_MF;
    sprintf(cosmo->status_message ,
	    "ccl_massfunc.c: ccl_massfunc(): Unknown or non-implemented mass function method: %d \n",
	    cosmo->config.mass_function_method);
    return 0;
  }
}
static double ccl_halo_b1(ccl_cosmology *cosmo, double halomass, double a, double odelta, int * status)
{
  double fit_A, fit_B, fit_C, fit_a, fit_b, fit_c, overdensity_delta, y;
  double delta_c_Tinker, nu;
  double sigma=ccl_sigmaM(cosmo,halomass,a, status);
  switch(cosmo->config.mass_function_method){

    //this version uses b(nu) parameterization, Eq. 6 in Tinker et al. 2010
    // use this for consistency with Tinker et al. 2010 fitting function for halo bias
  case ccl_tinker10:
    
    y = log10(odelta);
    //critical collapse overdensity assumed in this model
    delta_c_Tinker = 1.686;
    //peak height - note that this factorization is incorrect for e.g. massive neutrino cosmologies
    nu = delta_c_Tinker/(sigma);
    // Table 2 in https://arxiv.org/pdf/1001.3162.pdf
    fit_A = 1.0 + 0.24*y*exp(-pow(4./y,4.)); 
    fit_a = 0.44*y-0.88; 
    fit_B = 0.183; 
    fit_b = 1.5; 
    fit_C = 0.019+0.107*y+0.19*exp(-pow(4./y,4.)); 
    fit_c = 2.4; 

    return 1.-fit_A*pow(nu,fit_a)/(pow(nu,fit_a)+pow(delta_c_Tinker,fit_a))+fit_B*pow(nu,fit_b)+fit_C*pow(nu,fit_c);
    break;

  default:
    *status = CCL_ERROR_MF;
    cosmo->status = 11;
    sprintf(cosmo->status_message ,
      "ccl_massfunc.c: ccl_halo_b1(): No b(M) fitting function implemented for mass_function_method: %d \n",
      cosmo->config.mass_function_method);
    return 0;
  }
}

void ccl_cosmology_compute_sigma(ccl_cosmology * cosmo, int *status)
{
    if(cosmo->computed_sigma)
        return;

    // create linearly-spaced values of the mass.
    int nm=LOGM_SPLINE_NM;
    double * m = ccl_linear_spacing(LOGM_SPLINE_MIN, LOGM_SPLINE_MAX, nm);
    if (m==NULL ||
        (fabs(m[0]-LOGM_SPLINE_MIN)>1e-5) ||
        (fabs(m[nm-1]-LOGM_SPLINE_MAX)>1e-5) ||
        (m[nm-1]>10E17)
        ) {
       *status =CCL_ERROR_LINSPACE;
       strcpy(cosmo->status_message,"ccl_cosmology_compute_sigmas(): Error creating linear spacing in m\n");
       return;
    }
    
    // allocate space for y, to be filled with sigma and dlnsigma_dlogm
    double *y = malloc(sizeof(double)*nm);
    double smooth_radius; 
   
   // fill in sigma
   for (int i=0; i<nm; i++){
     smooth_radius = ccl_massfunc_m2r(cosmo, pow(10,m[i]), status);
     y[i] = log10(ccl_sigmaR(cosmo, smooth_radius));
   }
   gsl_spline * logsigma = gsl_spline_alloc(M_SPLINE_TYPE, nm);
   *status = gsl_spline_init(logsigma, m, y, nm);
   if (*status){
     free(m);
     free(y);
     gsl_spline_free(logsigma);
     *status = CCL_ERROR_SPLINE ;
     strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_sigma(): Error creating sigma(M) spline\n");
     return;
   }
   for (int i=0; i<nm; i++){
     if(i==0){
       y[i] = log(pow(10, gsl_spline_eval(logsigma, m[i], NULL)))-log(pow(10,gsl_spline_eval(logsigma, m[i]+LOGM_SPLINE_DELTA/2., NULL)));
       y[i] = 2.*y[i] / LOGM_SPLINE_DELTA;
     }
     else if (i==nm-1){
       y[i] = log(pow(10, gsl_spline_eval(logsigma, m[i]-LOGM_SPLINE_DELTA/2., NULL)))-log(pow(10,gsl_spline_eval(logsigma, m[i], NULL)));
       y[i] = 2.*y[i] / LOGM_SPLINE_DELTA;
     }
     else{
       y[i] = (log(pow(10,gsl_spline_eval(logsigma, m[i]-LOGM_SPLINE_DELTA/2., NULL)))-log(pow(10,gsl_spline_eval(logsigma, m[i]+LOGM_SPLINE_DELTA/2., NULL))));
       y[i] = y[i] / LOGM_SPLINE_DELTA;
     }
   }

   gsl_spline * dlnsigma_dlogm = gsl_spline_alloc(M_SPLINE_TYPE, nm);
   *status = gsl_spline_init(dlnsigma_dlogm, m, y, nm);
   if (*status){
     free(m);
     free(y);
     gsl_spline_free(logsigma);
     *status = CCL_ERROR_SPLINE ;
     strcpy(cosmo->status_message, "ccl_massfunc.c: ccl_cosmology_compute_sigma(): Error creating dlnsigma/dlogM spline\n");
     return;
   }

   if(cosmo->data.accelerator_m==NULL)
     cosmo->data.accelerator_m=gsl_interp_accel_alloc();
   cosmo->data.logsigma = logsigma;
   cosmo->data.dlnsigma_dlogm = dlnsigma_dlogm;
   cosmo->computed_sigma = true;

   free(m);
   free(y);
}

/*----- ROUTINE: ccl_massfunc -----
INPUT: ccl_cosmology * cosmo, double halo mass in units of Msun, double scale factor
TASK: returns halo mass function as dn / dlog10 m
*/

double ccl_massfunc(ccl_cosmology *cosmo, double halomass, double a, double odelta, int * status)
{
  if (!cosmo->computed_sigma){
    ccl_cosmology_compute_sigma(cosmo, status);
    ccl_check_status(cosmo, status);
  }

  double f,deriv,rho_m,logmass;
  
  logmass = log10(halomass);
  rho_m = RHO_CRITICAL*cosmo->params.Omega_m*cosmo->params.h*cosmo->params.h;
  f=massfunc_f(cosmo,halomass,a,odelta,status);
  ccl_check_status(cosmo, status);
  deriv = gsl_spline_eval(cosmo->data.dlnsigma_dlogm, logmass, cosmo->data.accelerator_m);
  return f*rho_m*deriv/halomass;
}

/*----- ROUTINE: ccl_halob1 -----
INPUT: ccl_cosmology * cosmo, double halo mass in units of Msun, double scale factor
TASK: returns linear halo bias
*/

double ccl_halo_bias(ccl_cosmology *cosmo, double halomass, double a, double odelta, int * status)
{
  if (!cosmo->computed_sigma){
    ccl_cosmology_compute_sigma(cosmo, status);
    ccl_check_status(cosmo, status);
  }

  double f;
  f = ccl_halo_b1(cosmo,halomass,a,odelta, status);  
  ccl_check_status(cosmo, status);  
  return f;
}
/*---- ROUTINE: ccl_massfunc_m2r -----
INPUT: ccl_cosmology * cosmo, halomass in units of Msun
TASK: takes halo mass and converts to halo radius
  in units of Mpc.
*/
double ccl_massfunc_m2r(ccl_cosmology * cosmo, double halomass, int * status)
{
    double rho_m, smooth_radius;

    //TODO: make this neater
    rho_m = RHO_CRITICAL*cosmo->params.Omega_m*cosmo->params.h*cosmo->params.h;

    smooth_radius = pow((3.0*halomass) / (4*M_PI*rho_m), (1.0/3.0));

    return smooth_radius;
}

/*----- ROUTINE: ccl_sigma_M -----
INPUT: ccl_cosmology * cosmo, double halo mass in units of Msun, double scale factor
TASK: returns sigma from the sigmaM interpolation. Also computes the sigma interpolation if
necessary.
*/

double ccl_sigmaM(ccl_cosmology * cosmo, double halomass, double a, int * status)
{
    double sigmaM;

    if (!cosmo->computed_sigma){
        ccl_cosmology_compute_sigma(cosmo, status);
        ccl_check_status(cosmo, status);
    }

    sigmaM = pow(10,gsl_spline_eval(cosmo->data.logsigma, log10(halomass), cosmo->data.accelerator_m));
    sigmaM = sigmaM*ccl_growth_factor(cosmo, a, status);

    return sigmaM;
}
