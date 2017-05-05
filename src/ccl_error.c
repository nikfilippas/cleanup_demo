#include "ccl_core.h"
#include "ccl_error.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

// RH 
void ccl_check_status(ccl_cosmology *cosmo, int * status){
	switch (*status){
		case 0: //all good, nothing to do
			return;
		case CCL_ERROR_LINSPACE:	// spacing allocation error, always terminate		
			fprintf(stderr,"%s",cosmo->status_message);
			exit(1);	
		case CCL_ERROR_SPLINE:	// spline allocation error, always terminate	
			fprintf(stderr,"%s",cosmo->status_message);
			exit(1);
		case CCL_ERROR_COMPUTECHI:	// compute_chi error //RH
			fprintf(stderr,"%s",cosmo->status_message);
			exit(1);
        case CCL_ERROR_HMF_INTERP: // terminate if hmf definition not supported
            fprintf(stderr,"%s",cosmo->status_message);
            exit(1);
        case CCL_ERROR_NU_INT: // error in getting the neutrino integral spline: exit. No status_message in cosmo because can't pass cosmology to the function. //DL
			fprintf(stderr, "%s", "Error, in ccl_neutrinos.c. ccl_calculate_nu_phasespace_spline(): Error in setting neutrino phasespace spline.");
			exit(1);


		// implement softer error handling, e.g. for integral convergence here
			
		default:		
			fprintf(stderr,"%s",cosmo->status_message);
			exit(1);
	}
}


/* ------- ROUTINE: ccl_check_status_nocosmo ------
   INPUTS: pointer to a status integer
   TASK: Perform a check on status for the case where it is not possible to have a cosmology object where the status check is required.
*/
void ccl_check_status_nocosmo(int * status){
	switch (*status){
		case 0: //all good, nothing to do
			return;
		case CCL_ERROR_LINSPACE:	// spacing allocation error, always terminate		
			fprintf(stderr,"%s", "CCL_ERROR_LINSPACE");
			exit(1);	
		case CCL_ERROR_SPLINE:	// spline allocation error, always terminate	
			fprintf(stderr,"%s", "CCL_ERROR_SPLINE");
			exit(1);
		case CCL_ERROR_COMPUTECHI:	// compute_chi error //RH
			fprintf(stderr,"%s", "CCL_ERROR_COMPUTECHI");
			exit(1);
        case CCL_ERROR_HMF_INTERP: // terminate if hmf definition not supported
            fprintf(stderr,"%s", "CCL_ERROR_HMF_INTERP");
            exit(1);
        case CCL_ERROR_NU_INT: // error in getting the neutrino integral spline: exit. No status_message in cosmo because can't pass cosmology to the function. //DL
			fprintf(stderr, "%s", "Error, in ccl_neutrinos.c. ccl_calculate_nu_phasespace_spline(): Error in setting neutrino phasespace spline.");
			exit(1);
			
		default:		
			fprintf(stderr,"%s", "OTHER ERROR");
			exit(1);
	}
}
