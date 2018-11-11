#include "ccl.h"
#include "../include/ccl_params.h"
#include "ctest.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SZ_VAL 0.4 //This will cancel the magnification contribution
#define CLS_TOLERANCE 1E-3
#define ELS_TOLERANCE 0.1
#define NELLS 3001

CTEST_DATA(cls) {
  double Omega_c;
  double Omega_b;
  double h;
  double A_s;
  double n_s;
  double sigma8;
};

CTEST_SETUP(cls) {
  data->Omega_c = 0.30;
  data->Omega_b = 0.00;
  data->h = 0.7;
  data->A_s = 2.1e-9;
  data->sigma8=0.8;
  data->n_s = 0.96;
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

static void compare_cls(char *compare_type,struct cls_data * data)
{
  int status=0;
  double zlss=1100.;


  ccl_configuration config = default_config;
  config.transfer_function_method = ccl_bbks;
  config.matter_power_spectrum_method = ccl_linear;
  ccl_parameters params = ccl_parameters_create_flat_lcdm(data->Omega_c,data->Omega_b,data->h,
							  data->A_s,data->n_s, &status);
  params.Omega_n_rel=0;
  params.Omega_l = 0.7;
  params.sigma8=data->sigma8;
  ccl_cosmology * cosmo = ccl_cosmology_create(params, config);
  ASSERT_NOT_NULL(cosmo);

  int nz;
  double *zarr_1,*pzarr_1,*zarr_2,*pzarr_2,*bzarr;
  if(!strcmp(compare_type,"analytic")) {
    //Create arrays for N(z)
    double zmean_1=1.0,sigz_1=0.15;
    double zmean_2=1.5,sigz_2=0.15;
    nz=512;
    zarr_1=malloc(nz*sizeof(double));
    pzarr_1=malloc(nz*sizeof(double));
    zarr_2=malloc(nz*sizeof(double));
    pzarr_2=malloc(nz*sizeof(double));
    bzarr=malloc(nz*sizeof(double));
    for(int ii=0;ii<nz;ii++) {
      double z1=zmean_1-5*sigz_1+10*sigz_1*(ii+0.5)/nz;
      double z2=zmean_2-5*sigz_2+10*sigz_2*(ii+0.5)/nz;
      double pz1=exp(-0.5*((z1-zmean_1)*(z1-zmean_1)/(sigz_1*sigz_1)));
      double pz2=exp(-0.5*((z2-zmean_2)*(z2-zmean_2)/(sigz_2*sigz_2)));
      zarr_1[ii]=z1;
      zarr_2[ii]=z2;
      pzarr_1[ii]=pz1;
      pzarr_2[ii]=pz2;
      bzarr[ii]=1.;
    }
  }
  else {
    char str[1024];
    char* rtn;
    int stat;
    FILE *fnz1=fopen("./tests/benchmark/codecomp_step2_outputs/bin1_histo.txt","r");
    ASSERT_NOT_NULL(fnz1);
    FILE *fnz2=fopen("./tests/benchmark/codecomp_step2_outputs/bin2_histo.txt","r");
    ASSERT_NOT_NULL(fnz2);
    nz=linecount(fnz1)-1; rewind(fnz1);
    zarr_1=malloc(nz*sizeof(double));
    pzarr_1=malloc(nz*sizeof(double));
    zarr_2=malloc(nz*sizeof(double));
    pzarr_2=malloc(nz*sizeof(double));
    bzarr=malloc(nz*sizeof(double));
    rtn = fgets(str,1024,fnz1);
    rtn = fgets(str,1024,fnz2);
    for(int ii=0;ii<nz;ii++) {
      double z1,z2,nz1,nz2;
      stat = fscanf(fnz1,"%lf %lf",&z1,&nz1);
      stat = fscanf(fnz2,"%lf %lf",&z2,&nz2);
      zarr_1[ii]=z1; zarr_2[ii]=z2;
      pzarr_1[ii]=nz1; pzarr_2[ii]=nz2;
      bzarr[ii]=1.;
    }
  }

  char fname[256];
  FILE *fi_dd_11,*fi_dd_12,*fi_dd_22;
  FILE *fi_dl_12;
  FILE *fi_dc_1,*fi_dc_2;
  FILE *fi_ll_11,*fi_ll_12,*fi_ll_22;
  FILE *fi_lc_1,*fi_lc_2;
  FILE *fi_cc;
  CCL_ClTracer *tr_nc_1=ccl_cl_tracer_number_counts_simple(cosmo,nz,zarr_1,pzarr_1,nz,zarr_1,bzarr,&status);
  ASSERT_NOT_NULL(tr_nc_1);
  CCL_ClTracer *tr_nc_2=ccl_cl_tracer_number_counts_simple(cosmo,nz,zarr_2,pzarr_2,nz,zarr_2,bzarr,&status);
  ASSERT_NOT_NULL(tr_nc_2);
  CCL_ClTracer *tr_wl_1=ccl_cl_tracer_lensing_simple(cosmo,nz,zarr_1,pzarr_1,&status);
  ASSERT_NOT_NULL(tr_wl_1);
  CCL_ClTracer *tr_wl_2=ccl_cl_tracer_lensing_simple(cosmo,nz,zarr_2,pzarr_2,&status);
  ASSERT_NOT_NULL(tr_wl_2);
  CCL_ClTracer *tr_cl=ccl_cl_tracer_cmblens(cosmo,zlss,&status);
  ASSERT_NOT_NULL(tr_cl);

  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b1%s_log_cl_dd.txt",compare_type);
  fi_dd_11=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dd_11);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b2%s_log_cl_dd.txt",compare_type);
  fi_dd_12=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dd_12);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b2b2%s_log_cl_dd.txt",compare_type);
  fi_dd_22=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dd_22);

  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b1%s_log_cl_dc.txt",compare_type);
  fi_dc_1=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dc_1);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b2b2%s_log_cl_dc.txt",compare_type);
  fi_dc_2=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dc_2);

  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b1%s_log_cl_ll.txt",compare_type);
  fi_ll_11=fopen(fname,"r"); ASSERT_NOT_NULL(fi_ll_11);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b2%s_log_cl_ll.txt",compare_type);
  fi_ll_12=fopen(fname,"r"); ASSERT_NOT_NULL(fi_ll_12);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b2b2%s_log_cl_ll.txt",compare_type);
  fi_ll_22=fopen(fname,"r"); ASSERT_NOT_NULL(fi_ll_22);
  sprintf(fname,"tests/benchmark/run_b1b2%s_log_cl_dl.txt",compare_type);
  fi_dl_12=fopen(fname,"r"); ASSERT_NOT_NULL(fi_dl_12);

  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b1b1%s_log_cl_lc.txt",compare_type);
  fi_lc_1=fopen(fname,"r"); ASSERT_NOT_NULL(fi_lc_1);
  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_b2b2%s_log_cl_lc.txt",compare_type);
  fi_lc_2=fopen(fname,"r"); ASSERT_NOT_NULL(fi_lc_2);

  sprintf(fname,"tests/benchmark/codecomp_step2_outputs/run_log_cl_cc.txt");
  fi_cc=fopen(fname,"r"); ASSERT_NOT_NULL(fi_cc);
  
  int *ells=malloc(NELLS*sizeof(int));
  double *cls_dd_11_b=malloc(NELLS*sizeof(double));
  double *cls_dd_12_b=malloc(NELLS*sizeof(double));
  double *cls_dd_22_b=malloc(NELLS*sizeof(double));
  double *cls_dl_11_b=malloc(NELLS*sizeof(double));
  double *cls_dl_12_b=malloc(NELLS*sizeof(double));
  double *cls_dl_22_b=malloc(NELLS*sizeof(double));
  double *cls_dc_1_b=malloc(NELLS*sizeof(double));
  double *cls_dc_2_b=malloc(NELLS*sizeof(double));
  double *cls_ll_11_b=malloc(NELLS*sizeof(double));
  double *cls_ll_12_b=malloc(NELLS*sizeof(double));
  double *cls_ll_22_b=malloc(NELLS*sizeof(double));
  double *cls_lc_1_b=malloc(NELLS*sizeof(double));
  double *cls_lc_2_b=malloc(NELLS*sizeof(double));
  double *cls_cc_b=malloc(NELLS*sizeof(double));
  double *cls_dd_11_h=malloc(NELLS*sizeof(double));
  double *cls_dd_12_h=malloc(NELLS*sizeof(double));
  double *cls_dd_22_h=malloc(NELLS*sizeof(double));
  double *cls_dl_11_h=malloc(NELLS*sizeof(double));
  double *cls_dl_12_h=malloc(NELLS*sizeof(double));
  double *cls_dl_22_h=malloc(NELLS*sizeof(double));
  double *cls_dc_1_h=malloc(NELLS*sizeof(double));
  double *cls_dc_2_h=malloc(NELLS*sizeof(double));
  double *cls_ll_11_h=malloc(NELLS*sizeof(double));
  double *cls_ll_12_h=malloc(NELLS*sizeof(double));
  double *cls_ll_22_h=malloc(NELLS*sizeof(double));
  double *cls_lc_1_h=malloc(NELLS*sizeof(double));
  double *cls_lc_2_h=malloc(NELLS*sizeof(double));
  double *cls_cc_h=malloc(NELLS*sizeof(double));

  for(int ii=0;ii<NELLS;ii++) {
    int l,stat;
    stat=fscanf(fi_dd_11,"%d %lf",&l,&(cls_dd_11_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_dd_12,"%d %lf",&l,&(cls_dd_12_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_dd_22,"%d %lf",&l,&(cls_dd_22_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_dl_12,"%d %lf",&l,&(cls_dl_12_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_dc_1,"%d %lf",&l,&(cls_dc_1_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_dc_2,"%d %lf",&l,&(cls_dc_2_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_ll_11,"%d %lf",&l,&(cls_ll_11_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_ll_12,"%d %lf",&l,&(cls_ll_12_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_ll_22,"%d %lf",&l,&(cls_ll_22_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_lc_1,"%d %lf",&l,&(cls_lc_1_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_lc_2,"%d %lf",&l,&(cls_lc_2_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    stat=fscanf(fi_cc,"%d %lf",&l,&(cls_cc_b[ii]));
    if(stat!=2) {
      fprintf(stderr,"Error reading benchmark file");
      exit(1);
    }
    ells[ii]=l;
  }

  fclose(fi_dd_11);
  fclose(fi_dd_12);
  fclose(fi_dd_22);
  fclose(fi_dl_12);
  fclose(fi_dc_1);
  fclose(fi_dc_2);
  fclose(fi_ll_11);
  fclose(fi_ll_12);
  fclose(fi_ll_22);
  fclose(fi_lc_1);
  fclose(fi_lc_2);
  fclose(fi_cc);

  double l_logstep = 1.05;
  double l_linstep = 20.;
  CCL_ClWorkspace *w=ccl_cl_workspace_new_limber(NELLS,l_logstep,l_linstep,&status);

  ccl_angular_cls(cosmo,w,tr_nc_1,tr_nc_1,NULL,NELLS,ells,cls_dd_11_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_1,tr_nc_2,NULL,NELLS,ells,cls_dd_12_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_2,tr_nc_2,NULL,NELLS,ells,cls_dd_22_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_1,tr_wl_1,NULL,NELLS,ells,cls_dl_11_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_1,tr_wl_2,NULL,NELLS,ells,cls_dl_12_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_1,tr_cl,NULL,NELLS,ells,cls_dc_1_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_2,tr_cl,NULL,NELLS,ells,cls_dc_2_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_nc_2,tr_wl_2,NULL,NELLS,ells,cls_dl_22_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_wl_1,tr_wl_1,NULL,NELLS,ells,cls_ll_11_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_wl_1,tr_wl_2,NULL,NELLS,ells,cls_ll_12_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_wl_2,tr_wl_2,NULL,NELLS,ells,cls_ll_22_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_wl_1,tr_cl,NULL,NELLS,ells,cls_lc_1_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_wl_2,tr_cl,NULL,NELLS,ells,cls_lc_2_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  ccl_angular_cls(cosmo,w,tr_cl,tr_cl,NULL,NELLS,ells,cls_cc_h,&status);
  if (status) printf("%s\n",cosmo->status_message);
  
  
  double fraction_failed=0;
  for(int ii=2;ii<w->n_ls-1;ii++) {
    int l=w->l_arr[ii];
    double ell_correct,ell_correct2;
    double el_dd_11,el_dd_12,el_dd_22;
    double el_dl_12;
    double el_dc_1,el_dc_2;
    double el_ll_11,el_ll_12,el_ll_22;
    double el_lc_1,el_lc_2;
    double el_cc;
    double cl_dd_11,cl_dd_12,cl_dd_22;
    double cl_dl_12;
    double cl_dc_1,cl_dc_2;
    double cl_ll_11,cl_ll_12,cl_ll_22;
    double cl_lc_1,cl_lc_2;
    double cl_cc;
    double cl_dd_11_h,cl_dd_12_h,cl_dd_22_h;
    double cl_dl_12_h;
    double cl_dc_1_h,cl_dc_2_h;
    double cl_ll_11_h,cl_ll_12_h,cl_ll_22_h;
    double cl_lc_1_h,cl_lc_2_h;
    double cl_cc_h;
    if(l<=1)
      ell_correct=1;
    else{
      ell_correct=l*(l+1.)/sqrt((l+2.)*(l+1.)*l*(l-1.));
      ell_correct2=(l+0.5)*(l+0.5)/sqrt((l+2.)*(l+1.)*l*(l-1.));
    }
    cl_dd_11  =cls_dd_11_b[l];
    cl_dd_12  =cls_dd_12_b[l];
    cl_dd_22  =cls_dd_22_b[l];
    cl_dl_12  =cls_dl_12_b[l];
    cl_dc_1  =cls_dc_1_b[l];
    cl_dc_2  =cls_dc_2_b[l];
    cl_ll_11  =cls_ll_11_b[l];
    cl_ll_12  =cls_ll_12_b[l];
    cl_ll_22  =cls_ll_22_b[l];
    cl_lc_1  =cls_lc_1_b[l];
    cl_lc_2  =cls_lc_2_b[l];
    cl_cc    =cls_cc_b[l];
    el_dd_11=fmax(ELS_TOLERANCE*sqrt((cl_dd_11*cl_dd_11+cl_dd_11*cl_dd_11)/(2*l+1.)),
		  cl_dd_11*CLS_TOLERANCE);
    el_dd_12=fmax(ELS_TOLERANCE*sqrt((cl_dd_11*cl_dd_22+cl_dd_12*cl_dd_12)/(2*l+1.)),
		  cl_dd_12*CLS_TOLERANCE);
    el_dd_22=fmax(ELS_TOLERANCE*sqrt((cl_dd_22*cl_dd_22+cl_dd_22*cl_dd_22)/(2*l+1.)),
		  cl_dd_22*CLS_TOLERANCE);
    el_dl_12=fmax(ELS_TOLERANCE*sqrt((cl_dd_11*cl_ll_22+cl_dl_12*cl_dl_12)/(2*l+1.)),
		  cl_dl_12*CLS_TOLERANCE);
    el_dc_1=fmax(ELS_TOLERANCE*sqrt((cl_dd_11*cl_cc+cl_dc_1*cl_dc_1)/(2*l+1.)),
		 cl_dc_1*CLS_TOLERANCE);
    el_dc_2=fmax(ELS_TOLERANCE*sqrt((cl_dd_22*cl_cc+cl_dc_2*cl_dc_2)/(2*l+1.)),
		 cl_dc_2*CLS_TOLERANCE);
    el_ll_11=fmax(ELS_TOLERANCE*sqrt((cl_ll_11*cl_ll_11+cl_ll_11*cl_ll_11)/(2*l+1.)),
		  cl_ll_11*CLS_TOLERANCE);
    el_ll_12=fmax(ELS_TOLERANCE*sqrt((cl_ll_11*cl_ll_22+cl_ll_12*cl_ll_12)/(2*l+1.)),
		  cl_ll_12*CLS_TOLERANCE);
    el_ll_22=fmax(ELS_TOLERANCE*sqrt((cl_ll_22*cl_ll_22+cl_ll_22*cl_ll_22)/(2*l+1.)),
		  cl_ll_22*CLS_TOLERANCE);
    el_lc_1=fmax(ELS_TOLERANCE*sqrt((cl_ll_11*cl_cc+cl_lc_1*cl_lc_1)/(2*l+1.)),
		 cl_lc_1*CLS_TOLERANCE);
    el_lc_2=fmax(ELS_TOLERANCE*sqrt((cl_ll_22*cl_cc+cl_lc_2*cl_lc_2)/(2*l+1.)),
		 cl_lc_2*CLS_TOLERANCE);
    el_cc=fmax(ELS_TOLERANCE*sqrt((cl_cc*cl_cc+cl_cc*cl_cc)/(2*l+1.)),
	       cl_cc*CLS_TOLERANCE);
    cl_dd_11_h=cls_dd_11_h[l];
    cl_dd_12_h=cls_dd_12_h[l];
    cl_dd_22_h=cls_dd_22_h[l];
    cl_dl_12_h=cls_dl_12_h[l]*ell_correct2;
    cl_dc_1_h=cls_dc_1_h[l];
    cl_dc_2_h=cls_dc_2_h[l];
    cl_ll_11_h=cls_ll_11_h[l]*ell_correct*ell_correct;
    cl_ll_12_h=cls_ll_12_h[l]*ell_correct*ell_correct;
    cl_ll_22_h=cls_ll_22_h[l]*ell_correct*ell_correct;
    cl_lc_1_h=cls_lc_1_h[l]*ell_correct;
    cl_lc_2_h=cls_lc_2_h[l]*ell_correct;
    cl_cc_h=cls_cc_h[l];

    ASSERT_TRUE(fabs(cl_dd_11_h-cl_dd_11)<el_dd_11);
    ASSERT_TRUE(fabs(cl_dd_12_h-cl_dd_12)<el_dd_12);
    ASSERT_TRUE(fabs(cl_dd_22_h-cl_dd_22)<el_dd_22);
    ASSERT_TRUE(fabs(cl_dl_12_h-cl_dl_12)<el_dl_12);
    ASSERT_TRUE(fabs(cl_dc_1_h-cl_dc_1)<el_dc_1);
    ASSERT_TRUE(fabs(cl_dc_2_h-cl_dc_2)<el_dc_2);
    ASSERT_TRUE(fabs(cl_ll_11_h-cl_ll_11)<el_ll_11);
    ASSERT_TRUE(fabs(cl_ll_12_h-cl_ll_12)<el_ll_12);
    ASSERT_TRUE(fabs(cl_ll_22_h-cl_ll_22)<el_ll_22);
    ASSERT_TRUE(fabs(cl_lc_1_h-cl_lc_1)<el_lc_1);
    ASSERT_TRUE(fabs(cl_lc_2_h-cl_lc_2)<el_lc_2);
    ASSERT_TRUE(fabs(cl_cc_h-cl_cc)<el_cc);
  }
  ccl_cl_workspace_free(w);
    
  free(ells);
  free(cls_dd_11_b); free(cls_dd_12_b); free(cls_dd_22_b); 
  free(cls_dl_12_b);
  free(cls_dc_1_b); free(cls_dc_2_b);
  free(cls_ll_11_b); free(cls_ll_12_b); free(cls_ll_22_b); 
  free(cls_lc_1_b); free(cls_lc_2_b);
  free(cls_cc_b);

  free(cls_dd_11_h); free(cls_dd_12_h); free(cls_dd_22_h); 
  free(cls_dl_12_h);
  free(cls_dc_1_h); free(cls_dc_2_h);
  free(cls_ll_11_h); free(cls_ll_12_h); free(cls_ll_22_h); 
  free(cls_lc_1_h); free(cls_lc_2_h);
  free(cls_cc_h);

  free(zarr_1);
  free(zarr_2);
  free(pzarr_1);
  free(pzarr_2);
  free(bzarr);
  ccl_cl_tracer_free(tr_nc_1);
  ccl_cl_tracer_free(tr_nc_2);
  ccl_cl_tracer_free(tr_wl_1);
  ccl_cl_tracer_free(tr_wl_2);
  ccl_cl_tracer_free(tr_cl);
  ccl_cosmology_free(cosmo);
}

CTEST2(cls,analytic) {
  compare_cls("analytic",data);
}

CTEST2(cls,histo) {
  compare_cls("histo",data);
}
