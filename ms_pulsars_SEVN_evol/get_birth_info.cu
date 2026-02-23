#include <stdio.h>
#include <stdlib.h>
#include "macro.h"
#include "initialize.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>

void get_info_birth(void *params){

	struct func_params *part= (struct func_params *)params;
	FILE *birth_info=NULL;
	FILE *taud_file=NULL;
	FILE *beta_file=NULL;
	FILE *betainit_file=NULL;
	FILE *subms_info=NULL;
	long np=0;
	birth_info=fopen("/home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_pop_ms_all.txt","r");
	taud_file=fopen("/home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_taud_all.txt","r");
	beta_file=fopen("/home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_beta_all.txt","r");
	betainit_file=fopen("/home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_beta_init_all.txt","r");
	subms_info=fopen("/home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_subms_all.txt","r");
	while(fscanf(birth_info,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",&part->period[np],&part->Pdot[np],&part->Edot[np],&part->alpha[np],&part->B[np],&part->age_pulsar[np],&part->Pinit[np],&part->Binit[np],&part->alpha0[np],&part->t_acc[np],&part->spin_up_or_not[np],&part->M_NS_f[np],&part->M_NS_i[np],&part->ellipticity[np],&part->Mc_i[np],&part->Mc_f[np],&part->a_bin[np],&part->ecc_bin[np],&part->comp_type[np],&part->comp_rem_type[np]) == 20) {np++;}
	printf("Check recovery of the data: %e\n",part->period[np-1]);
	np=0;
	/*for(np=0;np<part->Npulsars;np++){
		part->age_pulsar[np]=((5e9-5e7)*gsl_rng_uniform(part->r)+5e7)*365*24*3600.0;
	}       //When considering AIC of WD if I did not change already the age that I recover from SEVN
        np=0;*/
	while(fscanf(taud_file,"%le\n",&part->taud_evol[np]) == 1) {np++;}
	np=0;
	while(fscanf(beta_file,"%le\n",&part->beta[np]) == 1) {np++;}
	np=0;
	while(fscanf(betainit_file,"%le\n",&part->beta0[np]) == 1) {np++;}
	np=0;
	while(fscanf(subms_info,"%le\n",&part->subms[np]) == 1) {np++;}
	fclose(subms_info);
	fclose(birth_info);
	fclose(taud_file);
	fclose(beta_file);
	fclose(betainit_file);
}
