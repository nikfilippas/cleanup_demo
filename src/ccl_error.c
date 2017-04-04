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
                case CCL_ERROR_HMF_INTERP: // continue computation w/ Delta=200
                        fprintf(stderr,"%s",cosmo->status_message);
                        return;
        case CCL_ERROR_NU_INT: // error in getting the neutrino integral spline: exit. No status_message in cosmo because can't pass cosmology to the function. //DL
			fprintf(stderr, "%s", "Error, in ccl_neutrinos.c. ccl_calculate_nu_phasespace_spline(): Error in setting neutrino phasespace spline.");
			exit(1);

		// implement softer error handling, e.g. for integral convergence here
			
		default:		
			fprintf(stderr,"%s",cosmo->status_message);
			exit(1);
	}
}
