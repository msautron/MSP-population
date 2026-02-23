#include"initialize.h"
#include"macro.h"
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
void initialize(void *params){
       struct func_params *part= (struct func_params *)params;

       FILE *lines;
       lines=fopen("lines.txt","r");
       system("wc -l /home/matteo.sautron/Documents/ms_pulsars/ms_pulsars_SEVN_evol/data_pop_ms_all.txt > lines.txt");
       fscanf(lines,"%ld",&part->Npulsars);
       fclose(lines);

       part->v_old              =	70.0;//265.0;
       part->R			=	12000;//m
       part->zexp               =	0.18; // kpc
       part->Rexp		=	4.5; //kpc
       part->sigma_v		=	70.0;//265.; //  km/s 
       part->v_young		=	70.0;//265.; //  km/s 
       part->vacuum		=	0;
       part->vacuum_evol	=	0;
       part->ff_evol            =       1;
       part->Bfield_var         =       1;


       part->Pinit = (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pinit == NULL) printf("Pinit: allocation failed"); // check if allocation succeeded 
       part->Binit = (double *)calloc(part->Npulsars,sizeof(double)); // initialize pointer (allocate) 
                if (part->Binit == NULL) printf("Binit: allocation failed"); // check if allocation succeeded 
       part->alpha= (double *)calloc(part->Npulsars,sizeof(double)); // initialize pointer (allocate) 
                if (part->alpha== NULL) printf("alpha: allocation failed"); // check if allocation succeeded 
       part->age_pulsar= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->age_pulsar == NULL) printf("t_pulsar: allocation failed"); // check if allocation succeeded 
       part->period= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->period == NULL) printf("period: allocation failed"); // check if allocation succeeded 
       part->Edot= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Edot== NULL) printf("Edot: allocation failed"); // check if allocation succeeded 
       part->dist= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->dist== NULL) printf("dist: allocation failed"); // check if allocation succeeded 
       part->x= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x== NULL) printf("x: allocation failed"); // check if allocation succeeded 
       part->y= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y== NULL) printf("y: allocation failed"); // check if allocation succeeded 
       part->z= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z== NULL) printf("z: allocation failed"); // check if allocation succeeded 
       part->Pdot= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pdot== NULL) printf("Pdot: allocation failed"); // check if allocation succeeded 
       part->Fr= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fr== NULL) printf("Fr: allocation failed"); // check if allocation succeeded 
       part->Fg= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fg== NULL) printf("Fg: allocation failed"); // check if allocation succeeded 

       part->x0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x0== NULL) printf("x0: allocation failed"); // check if allocation succeeded 
       part->y0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y0== NULL) printf("y0: allocation failed"); // check if allocation succeeded 
       part->z0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z0== NULL) printf("z0: allocation failed"); // check if allocation succeeded 
       part->xi= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->xi== NULL) printf("xi: allocation failed"); // check if allocation succeeded 
       part->rho= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->rho== NULL) printf("rho: allocation failed"); // check if allocation succeeded 
       part->w_r_fast= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r_fast== NULL) printf("w_r_fast: allocation failed"); // check if allocation succeeded 
       part->w_r_pmps= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r_pmps== NULL) printf("w_r_pmps: allocation failed"); // check if allocation succeeded
       part->Smin_fast= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_fast== NULL) printf("Smin_fast: allocation failed"); // check if allocation succeeded 
       part->alpha0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->alpha0== NULL) printf("alpha0: allocation failed"); // check if allocation succeeded 
       part->B= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->B== NULL) printf("B: allocation failed"); // check if allocation succeeded 
       part->gl= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gl== NULL) printf("gl: allocation failed"); // check if allocation succeeded 
       part->gb= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gb== NULL) printf("gb: allocation failed"); // check if allocation succeeded 

       part->vx0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->vx0== NULL) printf("vx0: allocation failed");	// check if allocation succeede

       part->vy0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy0== NULL) printf("vy0: allocation failed"); // check if allocation succeede

       part->vz0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz0== NULL) printf("vz0: allocation failed"); // check if allocation succeede
       part->vz= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz== NULL) printf("vz: allocation failed"); // check if allocation succeede

       part->vy= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy== NULL) printf("vy: allocation failed"); // check if allocation succeede
       part->vx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vx== NULL) printf("vx: allocation failed"); // check if allocation succeede

       part->err_rel_g=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->err_rel_g== NULL) printf("err_rel_g: allocation failed"); // check if allocation succeede
       part->detec= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec== NULL) printf("detec: allocation failed"); // check if allocation succeeded
       part->detec_rad= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rad== NULL) printf("detec_rad: allocation failed"); // check if allocation succeeded
       part->detec_gam= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_gam== NULL) printf("detec_gam: allocation failed"); // check if allocation succeeded
       part->detec_rg= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rg== NULL) printf("detec_rg: allocation failed"); // check if allocation succeeded
       part->n_omega_z= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_z== NULL) printf("n_omega_z: allocation failed"); // check if allocation succeed
       part->n_omega_y= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_y== NULL) printf("n_omega_y: allocation failed"); // check if allocation succeeded
       part->n_omega_x= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_x== NULL) printf("n_omega_x: allocation failed"); // check if allocation succeeded
       part->PA= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->PA== NULL) printf("PA: allocation failed");
       part->DM=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->DM== NULL) printf("DM: allocation failed");
       part->t_acc=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->t_acc== NULL) printf("t_acc: allocation failed");
       part->M_NS_f=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->M_NS_f== NULL) printf("M_NS_f: allocation failed");
       part->M_NS_i=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->M_NS_i== NULL) printf("M_NS_i: allocation failed");
       part->spin_up_or_not=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->spin_up_or_not== NULL) printf("spin_up_or_not: allocation failed");
       part->Mc_f=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Mc_f== NULL) printf("Mc_f: allocation failed");
       part->a_bin=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->a_bin== NULL) printf("a_bin: allocation failed");
       part->ecc_bin=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ecc_bin== NULL) printf("ecc_bin: allocation failed");
       part->comp_type=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->comp_type== NULL) printf("comp_type: allocation failed");
       part->comp_rem_type=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->comp_rem_type== NULL) printf("comp_rem_type: allocation failed");
       part->Nb_orb= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Nb_orb== NULL) printf("Nb_orb: allocation failed");
       part->delta= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->delta== NULL) printf("delta: allocation failed");
       part->Mc_i=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Mc_i== NULL) printf("Mc_i: allocation failed");
       part->ellipticity=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ellipticity== NULL) printf("ellipticity: allocation failed");
       part->taud_evol=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->taud_evol== NULL) printf("taud_evol: allocation failed");
       part->temp=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->temp== NULL) printf("temp: allocation failed");
       part->Smin_pmps= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_pmps== NULL) printf("Smin_pmps: allocation failed");
       part->Smin_fermi= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_fermi== NULL) printf("Smin_fermi: allocation failed");	
       part->beta= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->beta== NULL) printf("beta: allocation failed");
       part->beta0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->beta0== NULL) printf("beta0: allocation failed");
       part->fomega=(double **)calloc(91,sizeof(double *));
                if (part->fomega== NULL) printf("fomega: allocation failed");
       for(int i=0;i<91;i++) {part->fomega[i]=(double *)calloc(91,sizeof(double)); if (part->fomega[i]==NULL) printf("fomega[i]: allocation failed");}
       part->Smin_SKAlow= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->Smin_SKAlow== NULL) printf("Smin_SKAlow: allocation failed");
       part->Smin_SKAmid= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_SKAmid== NULL) printf("Smin_SKAmid: allocation failed");
       part->Lgamma= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Lgamma== NULL) printf("Lgamma: allocation failed");
       part->subms= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->subms== NULL) printf("subms: allocation failed");
       part->w_int= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_int== NULL) printf("w_int: allocation failed");
       part->w_r_htru= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r_htru== NULL) printf("w_r_htru: allocation failed");
       part->Smin_htru= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_htru== NULL) printf("Smin_htru: allocation failed");
}
