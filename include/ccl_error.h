#pragma once
#include "ccl_core.h"

// Whether to do bounds checks on interpolated quantities
#define CCL_BOUNDS_CHECK_INTERP

#define CCL_ERROR_MEMORY 1025
#define CCL_ERROR_LINSPACE 1026
#define CCL_ERROR_INCONSISTENT 1027
#define CCL_ERROR_SPLINE 1028
#define CCL_ERROR_SPLINE_EV 1029
#define CCL_ERROR_INTEG 1030
#define CCL_ERROR_ROOT 1031
#define CCL_ERROR_CLASS 1032
#define CCL_ERROR_COMPUTECHI 1033
#define CCL_ERROR_MF 1034
#define CCL_ERROR_HMF_INTERP 1035
#define CCL_ERROR_PARAMETERS 1036
#define CCL_ERROR_NU_INT 1037

typedef enum {
    CCL_ERROR_POLICY_EXIT = 0,
    CCL_ERROR_POLICY_CONTINUE = 1,
} CCLErrorPolicy;

void ccl_raise_exception(int err, char* msg);
void ccl_set_error_policy(CCLErrorPolicy error_policy);
void ccl_check_status(ccl_cosmology *cosmo, int* status);
void ccl_check_status_nocosmo(int* status);
