#include "ccl.h"
#include "ctest.h"
#include <stdio.h>
#include <math.h>

#define EH_TOLERANCE 1.0E-4

CTEST_DATA(eh) {
  double Omega_c;
  double Omega_b;
  double h;
  double A_s;
  double n_s;
  double sigma_8;
  double Omega_n;
  double Omega_v[1];
  double Omega_k[1];
  double w_0[1];
  double w_a[1];
};

CTEST_SETUP(eh) {
  data->Omega_c = 0.25;
  data->Omega_b = 0.05;
  data->h = 0.7;
  data->A_s = 2.1e-9;
  data->sigma_8=0.8;
  data->n_s = 0.96;

  double Omega_v[1]={0.7};
  double w_0[1] = {-1.0};
  double w_a[1] = {0.0};

  for(int i=0;i<1;i++) {
    data->Omega_v[i] = Omega_v[i];
    data->w_0[i] = w_0[i];
    data->w_a[i] = w_a[i];
    data->Omega_n = 0.0;
    data->Omega_k[i] = 1.0 - data->Omega_c - data->Omega_b - data->Omega_n - data->Omega_v[i];
  }
}

static int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

static void compare_eh(int i_model,struct eh_data * data)
{
  int nk,i,j;
  char fname[256],str[1024];
  FILE *f;
  ccl_configuration config = default_config;
  config.transfer_function_method = ccl_eisenstein_hu;
  ccl_parameters params = ccl_parameters_create(data->Omega_c,data->Omega_b,
						data->Omega_k[i_model-1],data->Omega_n,
						data->w_0[i_model-1],data->w_a[i_model-1],
						data->h,data->A_s,data->n_s,-1,NULL,NULL);
  params.sigma_8=data->sigma_8;
  params.Omega_g=0;
  ccl_cosmology * cosmo = ccl_cosmology_create(params, config);
  ASSERT_NOT_NULL(cosmo);
  
  sprintf(fname,"./tests/benchmark/model%d_pk_eh_david.txt",i_model);
  f=fopen(fname,"r");
  if(f==NULL) {
    fprintf(stderr,"Error opening file %s\n",fname);
    exit(1);
  }
  nk=linecount(f)-1; rewind(f);
  
  fgets(str, 1024, f);
  for(i=0;i<nk;i++) {
    double k_h,k;
    int stat;
    stat=fscanf(f,"%lf",&k_h);
    if(stat!=1) {
      fprintf(stderr,"Error reading file %s, line %d\n",fname,i+2);
      exit(1);
    }
    k=k_h*data->h;
    for(j=0;j<1;j++) {
      double pk_h,pk_bench,pk_ccl,err;
      double z=2*j+0.;
      int status=0;
      stat=fscanf(f,"%lf",&pk_h);
      if(stat!=1) {
	fprintf(stderr,"Error reading file %s, line %d\n",fname,i+2);
	exit(1);
      }
      pk_bench=pk_h/pow(data->h,3);
      pk_ccl=ccl_linear_matter_power(cosmo,1./(1+z),k,&status);
      if (status) printf("%s\n",cosmo->status_message);
      err=fabs(pk_ccl/pk_bench-1);
      ASSERT_DBL_NEAR_TOL(err,0.,EH_TOLERANCE);
    }
  }
  fclose(f);

  ccl_cosmology_free(cosmo);
}

CTEST2(eh,model_1) {
  int model=1;
  compare_eh(model,data);
}

/*
CTEST2(eh,model_2) {
  int model=2;
  compare_eh(model,data);
}

CTEST2(eh,model_3) {
  int model=3;
  compare_eh(model,data);
}
*/
