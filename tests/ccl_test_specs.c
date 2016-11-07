#include "ccl_lsst_specs.h"
#include <stdio.h>
#include <math.h>

int main(int argc,char **argv){

	double z_test, a_test;
	double dNdz_tomo;
	double dNdz_clust;
	double sigz_src;
	double sigz_clust;
	double clust_bias;	
	int z;
	double Omega_c = 0.25;
        double Omega_b = 0.05;
        double h = 0.7;
        double A_s = 2.1e-9;
        double n_s = 0.96;
	FILE * output;
	ccl_cosmology* cosmo_1;

	// Declare a struct to pass parameters to the function which calculates the unnormalised dNdz of sources
	// The parameter to set here is 'type' - this sets the choice of dNdz for lensing sources. These are all from Chang et al 2013. type = 1 is k=0.5, type = 2 is k=1, type = 3 is k=2.
	struct dNdz_sources_params * p_test, val_params;
	val_params.type_ = 1;  // Set type 
	p_test = &val_params;

	/* Call the dNdz for weak lensing sources, the clustering dNdz, 
           and the photometric redshift errors for src and clustering 
           galaxies at each z and save the output */
	// This is to test that the functions can be called and give something reasonable.

	// Open file for writing:
	output = fopen("./tests/specs_test.out", "w");	

	// Test also the function for the bias in the clustering sample 
	// (requires setting up a cosmology to get the growth rate)
	ccl_parameters params = ccl_parameters_create_flat_lcdm(Omega_c, Omega_b, h, A_s, n_s);
	cosmo_1 = ccl_cosmology_create(params, default_config);

	//Try splitting dNdz (lensing) into 5 redshift bins
	double tmp1,tmp2,tmp3,tmp4,tmp5;
	output = fopen("./tests/specs_test_tomo_lens.out", "w");     
	for (z=0; z<100; z=z+1){
		z_test = 0.035*z;
		dNdz_tomo = ccl_specs_dNdz_tomog(z_test, p_test, 0.,6.,dNdz_sources_unnormed,sigmaz_sources); 
		tmp1 = ccl_specs_dNdz_tomog(z_test, p_test, 0.,0.6,dNdz_sources_unnormed,sigmaz_sources); 
		tmp2 = ccl_specs_dNdz_tomog(z_test, p_test, 0.6,1.2,dNdz_sources_unnormed,sigmaz_sources); 
		tmp3 = ccl_specs_dNdz_tomog(z_test, p_test, 1.2,1.8,dNdz_sources_unnormed,sigmaz_sources); 
		tmp4 = ccl_specs_dNdz_tomog(z_test, p_test, 1.8,2.4,dNdz_sources_unnormed,sigmaz_sources); 
		tmp5 = ccl_specs_dNdz_tomog(z_test, p_test, 2.4,3.0,dNdz_sources_unnormed,sigmaz_sources); 
		fprintf(output, "%f %f %f %f %f %f %f\n", z_test,tmp1,tmp2,tmp3,tmp4,tmp5,dNdz_tomo);
	}

	fclose(output);

	//Try splitting dNdz (clustering) into 5 redshift bins
	output = fopen("./tests/specs_test_tomo_clu.out", "w");     
	for (z=0; z<100; z=z+1){
		z_test = 0.035*z;
		dNdz_tomo = ccl_specs_dNdz_tomog(z_test, NULL,0.,6.,dNdz_clustering,sigmaz_clustering); 
		tmp1 = ccl_specs_dNdz_tomog(z_test, NULL,0.,0.6,dNdz_clustering,sigmaz_clustering); 
		tmp2 = ccl_specs_dNdz_tomog(z_test, NULL,0.6,1.2,dNdz_clustering,sigmaz_clustering); 
		tmp3 = ccl_specs_dNdz_tomog(z_test, NULL,1.2,1.8,dNdz_clustering,sigmaz_clustering);  
		tmp4 = ccl_specs_dNdz_tomog(z_test, NULL,1.8,2.4,dNdz_clustering,sigmaz_clustering); 
		tmp5 = ccl_specs_dNdz_tomog(z_test, NULL,2.4,3.0,dNdz_clustering,sigmaz_clustering); 
		fprintf(output, "%f %f %f %f %f %f %f\n", z_test,tmp1,tmp2,tmp3,tmp4,tmp5,dNdz_tomo);
	}

	fclose(output);


}


