#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>
 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>

#include "../cosmolike_light/theory/baryons.h"
#include "../cosmolike_light/theory/basics.c"
#include "../cosmolike_light/theory/structs.c"
#include "../cosmolike_light/theory/parameters.c"
#include "../cosmolike_light/emu17/P_cb/emu.c"
#include "../cosmolike_light/theory/recompute.c"
#include "../cosmolike_light/theory/cosmo3D.c"
#include "../cosmolike_light/theory/redshift_spline.c"
#include "../cosmolike_light/theory/halo.c"
#include "../cosmolike_light/theory/HOD.c"
#include "../cosmolike_light/theory/cosmo2D_fourier.c"
#include "../cosmolike_light/theory/IA.c"
#include "../cosmolike_light/theory/BAO.c"
#include "../cosmolike_light/theory/external_prior.c"
#include "../cosmolike_light/theory/covariances_3D.c"
#include "../cosmolike_light/theory/covariances_fourier.c"
#include "../cosmolike_light/theory/covariances_real.c"
#include "../cosmolike_light/theory/run_covariances_real.c"
#include "../cosmolike_light/theory/init.c"


int main(int argc, char** argv)
{
  int hit=atoi(argv[1]);
  FILE *F1,*F2;
  int i,l,m,n,o,s,p,output;
  double ktmp;
  char OUTFILE[400],filename[400];

  printf("-----------------\n");
  printf("This is not production mode covariance code\n");
  printf("Please only use with care and at your own risk\n");
  printf("-----------------\n");
  
  Ntable.N_a=20;
  set_cov_parameters_to_("cov_Y1/cov_y1_mcal_revision.ini",1);
  //here: setting values internally

  // set this to zero to quickly run Gaussian-only covariances for testing
  if (covparams.ng==1){
    NG = 1;
  }
  else {
    NG = 0;
  }

  // set this to one to output details about inputs for diagnostics
  output = 0;
  FILE *F;
  printf("running multi_covariance_real with NG = %d\n",NG);
  
  set_cosmological_parameters_to_("cov_Y1/cov_y1_mcal_revision.ini",1);
  //here: setting values internally 
  // cosmology.Omega_m   = 0.286;
  // cosmology.Omega_v   = 1.0-cosmology.Omega_m;
  // cosmology.sigma_8   = 0.82;
  // cosmology.n_spec    = 0.96;
  
  // cosmology.w0=-1.;
  // cosmology.wa=0.;
  // cosmology.omb=0.04868;
  // cosmology.h0=0.673;
  // cosmology.coverH0= 2997.92458; 

  // cosmology.rho_crit = 7.4775e+21;
  // cosmology.f_NL = 0.0;
  // pdeltaparams.runmode="Halofit"

  set_survey_parameters_to_("cov_Y1/cov_y1_mcal_revision.ini",1);
  //here: setting values internally 
  //survey.area=1000.0;
  //survey.m_lim=
  //survey.name=
  //survey.Kcorrect_File=
  // redshift.shear_REDSHIFT_FILE=source.nz
  // redshift.clustering_REDSHIFT_FILE=lens.nz
  // survey.sourcephotoz="multihisto";
  // survey.lensphotoz="multihisto";
  // survey.galsample=
  // tomo.shear_Nbin=5;
  // redshift.clustering_histogram_zbins =
  // tomo.clustering_Nbin=5;
  // tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  // survey.sigma_e=0.24;
  // survey.ggl_overlap_cut=
  // tomo.n_source[0]=1.07535691109;
  // tomo.n_source[1]=1.1606931903;
  // tomo.n_source[2]=1.2275933389;
  // tomo.n_source[3]=1.23453951309;
  // tomo.n_source[4]=1.3387876584;
  
  // tomo.n_lens[0]=0.02143277;
  // tomo.n_lens[1]=0.05426833;
  // tomo.n_lens[2]=0.08070083;
  // tomo.n_lens[3]=0.05130111;
  // tomo.n_lens[4]=0.01515694;
  // for (i=0;i<tomo.clustering_Nbin:i++) survey.n_lens+=tomo.n_lens[i];
  // for (i=0;i<tomo.shear_Nbin:i++) survey.n_gal+=tomo.n_source[i];

  init_source_sample_();
  init_lens_sample_();

  //init_clusters();
  //init_IA("none", "GAMA");
  //printf("test values: %d, %d, %s",redshift.clustering_photoz,tomo.clustering_Nbin,redshift.clustering_REDSHIFT_FILE);
  // printf("end of setup in main\n");

  double theta_min= covparams.tmin;
  double theta_max= covparams.tmax;
  int Ntheta=covparams.ntheta; 

  double logdt=(log(theta_max)-log(theta_min))/Ntheta;
  double *theta,*thetamin,*thetamax, *dtheta;
  theta=create_double_vector(0,Ntheta-1);
  thetamin=create_double_vector(0,Ntheta);
  thetamax=create_double_vector(0,Ntheta-1);
  dtheta=create_double_vector(0,Ntheta-1);
  for(i=0; i<Ntheta ; i++){
    thetamin[i]=exp(log(theta_min)+(i+0.0)*logdt);
    thetamax[i]=exp(log(theta_min)+(i+1.0)*logdt);
    theta[i] = 2./3.*(pow(thetamax[i],3.)-pow(thetamin[i],3.))/(pow(thetamax[i],2.)-pow(thetamin[i],2.));
    dtheta[i]=thetamax[i]-thetamin[i];
  }
  thetamin[Ntheta] = thetamax[Ntheta-1];

  int k=1;
  if (strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_ssss_++_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);

    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,1,1,k);}
        }
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_--_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,0,0,k);}  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_+-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,1,0,k);}  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0)
  {
    sprintf(OUTFILE,"%s_llll_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){ 
      for (m=l;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_clustering_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,k);}  
        }        
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_lsls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=l;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,k);}  
        }        
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ls,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_lsss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_ggl_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,1,k);}  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_lsss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_ggl_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,0,k);}  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0 && strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_llss_+_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_clustering_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,1,k);}  
        } 
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_llss_-_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_clustering_shear_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,0,k);}  
        } 
        k=k+1;
      }
    }
  }
  if (strcmp(covparams.ll,"true")==0 && strcmp(covparams.ls,"true")==0)
  {
    sprintf(OUTFILE,"%s_llls_cov_Ntheta%d_Ntomo%d",covparams.filename,Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {run_cov_clustering_ggl_real_bin(OUTFILE,covparams.outdir,thetamin,Ntheta,l,m,k);}  
        }
        k=k+1;
      }
    }
  }

  if (hit==0)
  {
    sprintf(OUTFILE,"%s%s",covparams.outdir,covparams.filename);
    write_gglensing_zbins(OUTFILE);

    sprintf(OUTFILE,"%s%s.blocks",covparams.outdir,covparams.filename);
    F1 = fopen(OUTFILE,"w");
    fprintf(F1,"%d\n",k-1);
    fclose(F1);
  }

  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}





