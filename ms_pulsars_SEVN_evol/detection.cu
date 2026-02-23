#include"macro.h"
#include <stdlib.h>
#include<string.h>
#include <stdio.h>
#include"initialize.h"
#include<gsl/gsl_rng.h>
#include <math.h>
#include "detection.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include<time.h>
#include <stdbool.h>
#include "cn.h"

double pdf_density_rho(double r){

        double A=37.6;double a=1.64;double b=4.01;double R_1=0.55;double R_sun=8.5; // Yusifov & Kucuk (2004)
        double rho=A*pow((r+R_1)/(R_sun+R_1),a)*exp(-b*((r-R_sun)/(R_sun+R_1)));
        return rho;

}

double theta_arms(double r,int choice){

        //Faucher-Giguère & Kaspi (2006)
        /*double k1=4.25;double r0_1=3.48;double theta0_1=1.57; //Norma arm
        double k2=4.25;double r0_2=3.48;double theta0_2=4.71; //Carina-Sagittarius
        double k3=4.89;double r0_3=4.90;double theta0_3=4.09; //Perseus
        double k4=4.89;double r0_4=4.90;double theta0_4=0.95; //Crux-Scutum */
        //Yao et al. (2017)
        double k1=4.95;double r0_1=3.35;double theta0_1=0.77; //Norma arm
        double k2=5.46;double r0_2=3.56;double theta0_2=3.82; //Carina-Sagittarius
        double k3=5.77;double r0_3=3.71;double theta0_3=2.09; //Perseus
        double k4=5.37;double r0_4=3.67;double theta0_4=5.76; //Crux-Scutum
        double theta;
        if (choice==1) {theta=k1*log(r/r0_1)+theta0_1;}
        else if (choice==2) {theta=k2*log(r/r0_2)+theta0_2;}
        else if (choice==3) {theta=k3*log(r/r0_3)+theta0_3;}
        else if (choice==4) {theta=k4*log(r/r0_4)+theta0_4;}
        return theta;

}

void distrib_init_2(void *params){
        struct func_params *part= (struct func_params*)params;
        long np=0;
        double p,rng_0_1_d;int rng_0_1;
        FILE *save_coord=NULL;
        int choice=1;
        bool sample=false;double *r;double pdf_val;double comp_val;double *theta;
        r=(double *)calloc(part->Npulsars,sizeof(double));
        theta=(double *)calloc(part->Npulsars,sizeof(double));
        double p_plus_or_minus_1;
        double omega_MW=(2*M_PI)/(250e6);
        save_coord=fopen("save_coord.txt","w+");

                for(np=0;np<part->Npulsars;np++){
                        part->z0[np]    =  gsl_ran_exponential(part->r,part->zexp); //kpc
                        rng_0_1=(rand() % 1000000001);rng_0_1_d=rng_0_1;p=rng_0_1_d/1000000000.0;
                        if(p<0.5) part->z0[np] = part->z0[np];
                        else if(p>=0.5)  part->z0[np]  = -part->z0[np];
                        part->x[np]=0;
                        part->y[np]=0;
			while(sample==false){
                   		r[np]=11.3+gsl_ran_gaussian_ziggurat(part->r,1.80);
                   		if(r[np]>=0){
                           		sample=true;
                   		}
           		}
                        /*while (sample==false){

                           r[np]=20*gsl_rng_uniform(part->r);
                           pdf_val=pdf_density_rho(r[np]);
                           comp_val=37.6*gsl_rng_uniform(part->r);
                           if (comp_val<=pdf_val) {sample=true;}

                        }*/
                        sample=false;
                        choice=(int)(4*gsl_rng_uniform(part->r)+1);
                        theta[np]=theta_arms(r[np],choice);
                        //Correction
                        double arg_exp=0.35*r[np];double sigma_r=0.07;
                        double theta_corr=2*M_PI*gsl_rng_uniform(part->r)*exp(-arg_exp);
                        double theta_corr_2=omega_MW*part->age_pulsar[np]/(365*24*3600);
                        theta[np]+=theta_corr+theta_corr_2;
                        double r_corr=1.0+gsl_ran_gaussian_ziggurat(part->r,sigma_r);
                        r[np]=r[np]*r_corr;
                        p_plus_or_minus_1=gsl_rng_uniform(part->r);
                        if (p_plus_or_minus_1<0.5) {part->y[np]=r[np]*tan(theta[np])/pow(1.0+pow(tan(theta[np]),2),0.5);part->x[np]=pow(pow(r[np],2)-pow(part->y[np],2),0.5);}
                        else {part->y[np]=-r[np]*tan(theta[np])/pow(1.0+pow(tan(theta[np]),2),0.5);part->x[np]=-pow(pow(r[np],2)-pow(part->y[np],2),0.5);}
                        fprintf(save_coord,"%e %e\n",part->x[np],part->y[np]);
			part->x0[np]= part->x[np];
                        part->y0[np]= part->y[np];
                        part->z[np]= part->z0[np];
                        }
		fclose(save_coord);
                free(r);free(theta);
                }

FILE *kick(void *params){

    	struct func_params *part= (struct func_params*)params;
    	long np;
	double two_pi=2*M_PI;
	double v=-1;
	double cos_theta,phi;
	double vx,vy,vz;
	double x_s,y_s,z_s;
	double dx,dy,dz;
	const double kpc2km=3.0856775807e16; //kpc to km
	//const double kpc2m=3.0856775807e19; //kpc to km
	const double yr_sec=365*24*3600; // year to sec
	double RAD_2 = 180/M_PI;
	double r;
	double glr,gbr; // galactic latitude, galactic longitude in radians
	double age_pulsar_yr;
	FILE *fp_gal=NULL;
	if ( ( fp_gal = fopen("galactic_coord.dat","w+")) == NULL){
		fprintf(stderr,"Couldn't open file galactic_coord.dat \n");
	exit(-1);
		}
			
           	for(np=0;np<part->Npulsars;np++){  
  	 		age_pulsar_yr    =  part->age_pulsar[np]/yr_sec; // /!\   should use a local variable
			if (age_pulsar_yr<3*1e6) part->sigma_v=part->v_young;
				else part->sigma_v=part->v_old;

			//Computation of the velocity if we ignore the fact that pulsars are going in the same direction as the rotation axis
			/*vx = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v); //km/s
			vy = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
			vz = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
        		v  = sqrt(vx*vx + vy*vy + vz*vz); // is  v really useful??*/

			//Computation of the velocity if we do not ignore the fact above
			cos_theta   =  2*gsl_rng_uniform(part->r)-1;
                        phi         =  two_pi*gsl_rng_uniform(part->r);
                        part->n_omega_x[np]=sqrt(1-sq(cos_theta))*cos(phi);
			part->n_omega_y[np]=sqrt(1-sq(cos_theta))*sin(phi);
                        part->n_omega_z[np]=cos_theta;
			while (v<0) {
                                v=sqrt(8.0/M_PI)*part->sigma_v+gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
                        }
                        //v=sqrt(8.0/M_PI)*part->sigma_v+gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
                        vx=v*part->n_omega_x[np];
			part->vx0[np]=vx;
			part->vx[np]=vx;
                        vy=v*part->n_omega_y[np];
			part->vy0[np]=vy;
			part->vy[np]=vy;
                        vz=v*part->n_omega_z[np];
			part->vz0[np]=vz;
			part->vz[np]=vz;

			dx = vx*part->age_pulsar[np]/kpc2km;
			dy = vy*part->age_pulsar[np]/kpc2km;
			dz = vz*part->age_pulsar[np]/kpc2km;

                        // coordinates in the Galactocentric system 
			part->x[np]=part->x0[np]+dx; 
			part->y[np]=part->y0[np]+dy;
			part->z[np]=part->z0[np]+dz;
		        	
                        // coordinates relative to the Sun in kpc 
			x_s            = part->x[np]; // shift center of the Galaxy to the Sun
			y_s            = part->y[np]-8.5;
			z_s            = part->z[np]-0.015; 
    			r	       = sqrt(x_s*x_s+y_s*y_s);

			part->dist[np]= sqrt(sq(x_s)+sq(y_s)+sq(z_s));

			// Galactic coordinates 

			if(part->dist[np]<1e-15){
			     gbr=0;
			 }else{
			     gbr=asin(z_s/part->dist[np]);
			 }
			part->gb[np]=gbr*RAD_2;
			    
			if(r<1e-15){
			      glr=0;
			}else{
			      if(x_s>=0){
				glr=acos(-y_s/r);
			      }else{
				glr=acos(y_s/r)+M_PI;
			      }
			 }

			part->gl[np]=glr*RAD_2;
			v=-1;

			part->dist[np]=5*gsl_rng_uniform(part->r);
			//part->gb[np] = RAD*asin(part->z[np]/part->dist[np]);
			//part->gl[np] = 180*acos(xx/part->dist[np]/M_PI);
			//part->gb[np] = 180*asin(part->z[np]/part->dist[np]/M_PI);
			//part->gb[np] = 180*asin(part->z[np]/part->dist[np])/M_PI;


			fprintf(fp_gal,"Gal %e %e %e %d \n",part->gl[np],part->gb[np],part->dist[np],2);
			//fprintf(file,"%e %e\n",part->x[np],part->y[np]);
	//		printf("Gal %e %e %e %d \n",part->gl[np],part->gb[np],part->dist[np],2);
			//printf("%e \n",part->gl[np]);
			//printf("%e  %e \n",part->x[np],part->z[np]);
        	} 
		//fclose(file);


return fp_gal;

/*
    FILE *fp=NULL;
       if ( ( fp= fopen("file.dat","w+")) == NULL){
                printf("Couldn't open file file.dat \n");
                exit(-1);
        }
       int i;
    for (i=0;i<10;i++){
	x[i]=i*i;	
 	//printf("%d \n", x[i]);
 	fprintf(fp,"%d \n",x[i]);
 //	fprintf(fp,"%s \n", "This is file.dat");
    }
    return fp;
*/
}



/*int geometry(void *params){ // calculated the geometry of the pulsar (i.e xi, and rho) 


    	struct func_params *part= (struct func_params*)params;
    	long np=0;
	//double two_pi=2*M_PI;   
	//double phi,cos_theta;
	double Npulsars=part->Npulsars;
	double n[3];
	//double n_omega[3];
	double norm;
	double xi,alpha;
	double rho;
	double ratio;
	//double k=0.118940963;
	
	
           	for (np=0;np<Npulsars;np++){ 

			part->np       = np;
                        if (part->alpha[np]<=M_PI/2) alpha=part->alpha[np]; // just to make it easier to read    
                        else if (part->alpha[np]>M_PI/2) alpha=M_PI-part->alpha[np];                    
			//rho            = 3*sqrt(M_PI*hem/(2*part->period[np]*SI_C)); // replace by a numeric value
			//rho            = k/sqrt(part->period[np]); // 3*sqrt((PI*hem)/2*P*c) Equation 2 of Johnston et al. (2020)
			//part->rho[np]  = rho;
	 		//phi	       =   two_pi*gsl_rng_uniform(part->r); // azimuth of the rotation axis 
			//if(np%2==0)  cos_theta   =  gsl_rng_uniform(part->r); // theta of the rotation axis
		       	  //    else   cos_theta   =  -gsl_rng_uniform(part->r); // theta is between -1 and 1 
			//cos_theta   =  gsl_rng_uniform(part->r); // theta of the rotation axis

			// line of sight vector 
			norm=sqrt(sq(part->x[np])+sq(part->y[np]-8.5)+sq(part->z[np]-0.015));
			n[0]=(part->x[np])/norm;  
			n[1]=(part->y[np]-8.5)/norm;
			n[2]=(part->z[np]-0.015)/norm;

			// unit vector of the rotation axis 
			//n_omega[0]=((1-sq(cos_theta))*cos(phi));part->n_omega_x[np]=n_omega[0];
			n_omega[1]=((1-sq(cos_theta))*sin(phi));part->n_omega_y[np]=n_omega[1];
			n_omega[2]=cos_theta;part->n_omega_z[np]=n_omega[2];

			// angle between vectors n and n_omega	
		        xi=acos(part->n_omega_x[np]*n[0]+part->n_omega_y[np]*n[1]+part->n_omega_z[np]*n[2]); 
			part->xi[np]=xi;
			if (part->xi[np]<=M_PI/2) xi=part->xi[np]; // just to make it easier to read
                        else if (part->xi[np]>M_PI/2) xi=M_PI-part->xi[np];

			// Computation of rho for ms pulsars  
                        //part->rho[np]=acos(cos(alpha)*cos(xi)+sin(alpha)*sin(xi)*cos((2*M_PI/part->period[np])*(part->w_r[np]/2.0)));
                        //part->rho[np]=(M_PI/180.0)*(1.5*pow(part->period[np],-0.5));
                        double r_kg=40*pow(1.4,-0.26)*pow(part->Pdot[np]/1e-15,0.07)*pow(part->period[np],0.3);
                        part->rho[np]=(M_PI/180.0)*1.24*sqrt(r_kg)*pow(part->period[np],-0.5);rho=part->rho[np];


			 //calculates the width of the radio beam w_r, from eq 22 of our paper 
				if(fabs(alpha-xi)<= rho && alpha >= rho && xi < M_PI/2){ 
					ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
					part->w_r[np]=2*acos(ratio);
				} else if(fabs(xi-(M_PI-alpha))<=rho && alpha >= rho && xi > M_PI/2){ 
					//alpha=M_PI-alpha;
					ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
					part->w_r[np]=2*acos(ratio);		
				}	
		}

    return(0);
}*/


void spinvel_angle(void *params){

	struct func_params *part= (struct func_params*)params;
	long np=0;
	double norm_v;double norm_omega;
	for(np=0;np<part->Npulsars;np++){
	   norm_v=sqrt(sq(part->vx[np]*1e3)+sq(part->vy[np]*1e3)+sq(part->vz[np]*1e3));
	   norm_omega=sqrt(sq(part->n_omega_x[np])+sq(part->n_omega_y[np])+sq(part->n_omega_z[np]));
	   part->PA[np]=(part->n_omega_x[np]*(part->vx[np]*1e3)+part->n_omega_y[np]*(part->vy[np]*1e3)+part->n_omega_z[np]*(part->vz[np]*1e3))/(norm_v*norm_omega);
           //if (np==0) printf("value computed for PA: %e",part->PA[np]);
	}
	printf("Spinvel angle computed for eveyrone\n");

}

void save_all(void *params){

	struct func_params *part= (struct func_params*)params;
	FILE *file_data=NULL;
	file_data=fopen("P_Pdot_positions.txt","w+");
	for(long np=0;np<part->Npulsars;np++){
		fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|1|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->alpha0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np],part->t_acc[np],part->spin_up_or_not[np],part->M_NS_i[np],part->M_NS_f[np],part->Mc_f[np],part->a_bin[np],part->ecc_bin[np],part->comp_type[np],part->comp_rem_type[np]);
	}
	fclose(file_data);
}

void gamma_ray_peak_sep(void *params){

        struct func_params *part= (struct func_params*)params;
        long np=0;
        for(np=0;np<part->Npulsars;np++){
                part->delta[np]=(1.0/M_PI)*acos(fabs((cos(part->xi[np])/sin(part->xi[np]))*(cos(part->alpha[np])/sin(part->alpha[np]))));
                if(part->delta[np]>0.5) part->delta[np]=1.0-part->delta[np];
        }
	printf("Gamma-ray peak separation computed for everyone\n");
}

void count_nb_msp_formed(void *params){

	struct func_params *part= (struct func_params*)params;
        long np=0;
	int count_nb_msp_formed=0;
	double R_msp,BS_MSP;
	double age_max_msp=13.9e9;
	double age_min_msp=5e7;
	for(np=0;np<part->Npulsars;np++){
		if (part->spin_up_or_not[np]==1 && part->t_acc[np]>=3.1536e12 && part->period[np]<=0.6) count_nb_msp_formed++;
	}
	R_msp=count_nb_msp_formed/(age_max_msp-age_min_msp);
	BS_MSP=1.0/R_msp;
	printf("Number of MSP formed : %d\nBirth rate of MSP : %e 1/yr\nBirth spacing of MSP : %e yr\n(According to Story et al. (2007), BR MSP between 3e-6 and 6.5e-6)\n",count_nb_msp_formed,R_msp,BS_MSP);
}

void sky_temp_Fmin_fermi(void *params){

	struct func_params *part= (struct func_params*)params;
	long np=0;
	FILE *lb_coord=NULL;
	FILE *temp_data=NULL;
	FILE *fermi_fmin=NULL;
	lb_coord=fopen("l_b_coord_sim.txt","w+");
	temp_data=fopen("temp.txt","r");
	fermi_fmin=fopen("fermi_fmin.txt","r");
	for(np=0;np<part->Npulsars;np++){
		fprintf(lb_coord,"%e %e\n",part->gl[np],part->gb[np]);
	}
	np=0;
	fclose(lb_coord);
	system("python3 get_temp.py");
	system("python3 sensitivity_3PC.py");
	while(fscanf(temp_data,"%le\n",&part->temp[np])==1) {np++;}
	np=0;
	while(fscanf(fermi_fmin,"%le\n",&part->Smin_fermi[np])==1) {np++;}
	np=0;
	/*for(np=0;np<part->Npulsars;np++){
                part->temp[np]=part->temp[np]*pow(408.0/1374.0,2.6); //Adapt the temperature to fast survey 
        }*/
        fclose(temp_data);	
	fclose(fermi_fmin);
}

void pos_all_MSP(void *params){

	struct func_params *part= (struct func_params*)params;
        long np=0;
	FILE *pos_all=NULL;
	pos_all=fopen("pos_all.txt","w+");
	for(np=0;np<part->Npulsars;np++){
		if (part->spin_up_or_not[np]==1 && part->t_acc[np]>=3.1536e12 && part->period[np]<=0.6){
			fprintf(pos_all,"%e %e %e %e %e\n",part->x[np],part->y[np],part->z[np],part->Lgamma[np],part->gl[np]);
		}
	}
	fclose(pos_all);
}

void subms_info(void *params){

	struct func_params *part= (struct func_params*)params;
        long np=0;
	FILE *subms_info=NULL;
	subms_info=fopen("subms_info_complete.txt","w+");
	for(np=0;np<part->Npulsars;np++){
		fprintf(subms_info,"%e %e %e %e %e %ld %ld %ld\n",part->period[np],part->Pdot[np],part->M_NS_f[np],part->M_NS_i[np],part->subms[np],part->detec_rad[np],part->detec_gam[np],part->detec_rg[np]);
	}
	fclose(subms_info);

}

void save_info_AIC(void *params){

	struct func_params *part= (struct func_params*)params;
        long np=0;
        FILE *AIC_info=NULL;
        AIC_info=fopen("AIC_info.txt","w+");
	double P_orb=0;
	double RSUN=696430e3;
        for(np=0;np<part->Npulsars;np++){
		if ((part->detec_rad[np]==1 || part->detec_gam[np]==1 || part->detec_rg[np]==1) && (part->period[np] <= 6e-1) && (part->Pdot[np] <=5e-17)){
			P_orb=sqrt(4*sq(M_PI)*pow(part->a_bin[np]*RSUN,3.0)/(G_grav*1e9*(part->Mc_f[np]+part->M_NS_f[np])*MSUN))/(24*60*60); //P_orb in days
			if (!isnan(P_orb) && !isnan(part->period[np])){
                		fprintf(AIC_info,"%e %e\n",part->period[np],P_orb);
				}
        	}
	}
        fclose(AIC_info);
}

void compute_wr_htru_pks(void *params){

	struct func_params *part= (struct func_params*)params;
        long np=0;
	double ratio;
	double xi;
	double rho;
	double alpha;
	double tau_sc;
	double e=1.6e-19;//electronic charge in C
        double m_e=9.1e-31; //electron mass in kg
	double tau_samp_htru=64e-6;
	double tau_dm_htru;
	double f_htru=1.352e9;
	double deltaf_htru=390.625e3;
	double DM_SI;
	for(np=0;np<part->Npulsars;np++){

		if (part->xi[np]<=M_PI/2.0) xi=part->xi[np];
                else if (part->xi[np]>M_PI/2.0) xi=M_PI-part->xi[np];
                rho=part->rho[np];
                if (part->alpha[np]<=M_PI/2.0) alpha=part->alpha[np];
                else if (part->alpha[np]>M_PI/2.0) alpha=M_PI-part->alpha[np];
		if(fabs(alpha-xi)<= rho){
			ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                        part->w_int[np]=2*acos(ratio);
                } else if(fabs(xi-(M_PI-alpha))<=rho){
                        ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                        part->w_int[np]=2*acos(ratio);
                        }
		DM_SI=part->DM[np]*3.086e22;
		tau_sc=3.6e-9*pow(part->DM[np],2.2)*(1+1.94e-3*sq(part->DM[np])); //in s
		tau_dm_htru=(sq(e)*deltaf_htru*DM_SI)/(4*M_PI*SI_eps0*M_PI*m_e*SI_C*cube(f_htru));
		part->w_r_htru[np]=sqrt(sq(part->w_int[np]*part->period[np]/(2*M_PI))+sq(tau_samp_htru)+sq(tau_dm_htru)+sq(tau_sc))*((2*M_PI)/(part->period[np]));
	}

}

int detection(void *params){ //check the flux of each pulsar and if the beam sweps the Earth

	printf("Started detection pipeline\n");
    	struct func_params *part= (struct func_params*)params;
    	long np=0;
	long count_radio_tot=0;
	long count_radio=0;
	long count_gamma=0;
	long count_radio_gamma=0;
	long Nsup_gamma=0;
	long Nsup2_gamma=0;
	long Nsup=0;
	long Nsup2=0;
	long Nsup_radio=0;
	long Nsup_radio_gamma=0;
	long Nsup2_radio=0;
	long Nsup2_radio_gamma=0;
	long Nr=0;
	long Ng=0;
	long count_rg_lowF=0;
	double ratio;
	double xi,rho,alpha;
	double Smin_gamma;
        double glat;
	long Ntot;
	long Nbeam=0;
        int detec;
	int count_possible_gdetec=0;
	long S_Nmin_fast=9.0; // signal to noise ratio min detectable FAST
	long S_Nmin_pmps=9.0; // signal to noise ratio min detectable pmps
	long S_Nmin_htru=9.0; // signal to noise ratio min detectable htru
	FILE *file_data=NULL;
	FILE *file_taud=NULL;
	FILE *file_wr=NULL;
	FILE *file_wint=NULL;
	FILE *beta_file=NULL;
	FILE *betainit_file=NULL;
	FILE *save_tempo=NULL;
	FILE *pos_init=NULL;
	double eta=0.15;double alpha_l=45*(M_PI/180);double T6=2;double b=40;
	double alpha_l2;double T6_2;double b2;
	double P_dot_line;
	int detec_fast=0;int detec_pmps=0;int detec_htru=0;
	int current_detec_fast,current_detec_pmps,current_detec_htru;
	int count_dead=0;
	int count_dm_sup=0;
	FILE *gamma_peak_sep=NULL;
        gamma_peak_sep=fopen("gamma_peak_sep.txt","w+");
	file_wr=fopen("wr.txt","w+");
	beta_file=fopen("data_beta_detec_all.txt","w+");
	betainit_file=fopen("data_beta_detec_init_all.txt","w+");
	file_data=fopen("P_Pdot_positions.txt","w+");
	file_taud=fopen("data_taud_detec.txt","w+");
	save_tempo=fopen("xi_rho_data.txt","w+");
	pos_init=fopen("pos_init.txt","w+");
	file_wint=fopen("wint.txt","w+");
			
           	for (np=0;np<part->Npulsars;np++){	
			//printf("Pulsar number : %ld\n",np);
			if (part->xi[np]<=M_PI/2.0) xi=part->xi[np];
			else if (part->xi[np]>M_PI/2.0) xi=M_PI-part->xi[np];
			rho=part->rho[np];
			if (part->alpha[np]<=M_PI/2.0) alpha=part->alpha[np];
                        else if (part->alpha[np]>M_PI/2.0) alpha=M_PI-part->alpha[np];
			P_dot_line=(3.16e-4*pow(T6,4)*sq(part->period[np])*1e-15)/(sq(eta)*b*sq(cos(alpha_l)));
			if (P_dot_line/part->Pdot[np] >= pow(10,-0.55) && P_dot_line/part->Pdot[np] <= pow(10,1.15)) {
                                alpha_l2=65*(M_PI/180)*gsl_rng_uniform(part->r);
                                T6_2=(2.8-1.9)*gsl_rng_uniform(part->r)+1.9;
                                b2=30*gsl_rng_uniform(part->r)+30;
                                P_dot_line=(3.16e-4*pow(T6_2,4)*sq(part->period[np])*1e-15)/(sq(eta)*b2*sq(cos(alpha_l2)));
                        }

			/* calculates the width of the radio beam w_r, from eq 22 of our paper */
      			if(fabs(alpha-xi)<= rho){
                		ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                		part->w_int[np]=2*acos(ratio);
      			} else if(fabs(xi-(M_PI-alpha))<=rho){
                		ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                		part->w_int[np]=2*acos(ratio);
                		}

			glat = (180*asin(part->z[np]/part->dist[np]))/M_PI;
			detec=0; //Initialization of detec param for detection, it will be one if detected in one survey 
			current_detec_fast=0;current_detec_pmps=0;

			//SKA-low detection
			//if (part->Fr[np]/part->Smin_SKAlow[np] > S_Nmin_fast && part->Smin_SKAlow[np]!=0 && part->gb[np]>-90 && part->gb[np]<30) {detec = 1;detec_fast+=1;current_detec_fast=1;}

			//SKA-mid detection
                        //if (part->Fr[np]/part->Smin_SKAmid[np] > S_Nmin_pmps && part->Smin_SKAmid[np]!=0 && part->gb[np]>-90 && part->gb[np]<30) {detec = 1;detec_pmps+=1;current_detec_pmps=1;}

			//HTRU-mid detection 
                        //if (part->Fr[np]/part->Smin_htru[np] > S_Nmin_htru && part->Smin_htru[np]!=0 && (part->gl[np]<=30 || part->gl[np]>=240) && abs(part->gb[np])<15) {detec = 1;detec_htru+=1;current_detec_htru=1;}
			//FAST detection
                        if (part->Fr[np]/part->Smin_fast[np] > S_Nmin_fast && part->Smin_fast[np]!=0 && abs(part->gb[np])<10) {detec = 1;detec_fast+=1;current_detec_fast=1;}
			//PMPS detection
			if (part->Fr[np]/part->Smin_pmps[np] > S_Nmin_pmps && part->Smin_pmps[np]!=0 && (part->gl[np]<=50 || part->gl[np]>=260) && abs(part->gb[np])<5) {detec = 1;detec_pmps+=1;current_detec_pmps=1;}				

			/* radio detection */
				if(fabs(xi-alpha)<= rho){ 
					Nbeam++;
					//printf("Id : %ld, S/N low : %e, S/N mid : %e, detec=%d\n",np,part->Fr[np]/part->Smin_SKAlow[np],part->Fr[np]/part->Smin_SKAmid[np],detec);
                			if (detec==1){
						//printf("Got in!\n");
						Nr=1;  
						part->detec[np]=1;
                                                part->detec_rad[np]=1;
						count_radio_tot++;
					}
				} else if (fabs(xi-(M_PI-alpha))<=rho){ 
					Nbeam++;
					//printf("Id : %ld, S/N low : %e, S/N mid : %e, detec=%d\n",np,part->Fr[np]/part->Smin_SKAlow[np],part->Fr[np]/part->Smin_SKAmid[np],detec);
                			if (detec==1){
						//printf("Got in!\n");
						count_radio_tot++;
						part->detec[np]=1;
                                                part->detec_rad[np]=1;
						Nr=1;  
						}
					}
			
			/* gamma detection */
			  	if(fabs(xi-M_PI/2.)<=alpha){ 
					//if(part->Fg[np]>1e-16) {Ng=1;part->detec[np]=1;part->detec_gam[np]=1;} //Predictions CTA
					if(part->Fg[np]>part->Smin_fermi[np]) {Ng=1;part->detec[np]=1;part->detec_gam[np]=1;}
					if((part->Smin_fermi[np]/part->Fg[np])<=5 && part->Pdot[np]<1e-19) {count_possible_gdetec+=1;}
					//if(Nr==1 && glat < 2.) {  //deep follow-up surveys at 1.4 GHz of gamma-ray sources
					//    Smin_gamma=4e-15; // in W.m^-2 
					//} else Smin_gamma=16e-15;
				//Smin_gamma=5e-12;
                			//if(part->Fg[np]>Smin_gamma) {Ng=1;part->detec[np]=1;part->detec_gam[np]=1; }  
				}  //else continue; 

			/* X-rays detection */
				//if (part->Fx[np]>Smin_x){Nx=1;part->detec_x[np]=1;part->detec[np]=1;}

				if(Ng==1 && Nr==1){  //both radio and gamma are detected
				        if (P_dot_line>part->Pdot[np]) {Ng=0;Nr=0;part->detec_rg[np]=0;part->detec[np]=0;part->detec_rad[np]=0;part->detec_gam[np]=0;count_dead+=1;}
					else{
					   count_radio_gamma++;part->detec_rg[np]=1;part->detec[np]=1;part->detec_rad[np]=0;part->detec_gam[np]=0;	
		
			        	       	   if(part->Edot[np]>1e28) { //28
							   Nsup_radio_gamma+=1;
						   }
			        		   if(part->Edot[np]>1e31) { //31
							   Nsup2_radio_gamma+=1;
						   }
						   if(part->Fr[np]*1e-3<30e-6){
							   count_rg_lowF+=1;
						   }
						   if(part->DM[np]>=2.5){
							   count_dm_sup+=1;
						   }
					} 
				}

				else if(Nr==1 && Ng==0){ //radio only
                       
					if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_rad[np]=0;count_dead+=1;}
					else{
					   count_radio++;part->detec[np]=1;part->detec_rad[np]=1;

			        		   if(part->Edot[np]>1e28) { 
							   Nsup_radio+=1;
						   }
			        		   if(part->Edot[np]>1e31) { 
							   Nsup2_radio+=1;
						   }
				}
				}
				else if(Ng==1 && Nr==0){ //gamma only
				
					//if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_gam[np]=0;}
					//else{
					   count_gamma++;part->detec[np]=1;part->detec_gam[np]=1;

			        		   if(part->Edot[np]>1e28) { 
							   Nsup_gamma+=1;
						   } 
			        		   if(part->Edot[np]>1e31) { 	
							   Nsup2_gamma+=1;
						   }
						   if(part->DM[np]>=2.5){
                                                           count_dm_sup+=1;
                                                   }
				//} 
				}

				/*else if(Ng==0 && Nr==0 && Nx==1){ //X-rays only

                                        if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_x[np]=0;}
                                        else{
                                           count_x++;part->detec[np]=1;part->detec_x[np]=1;

                                                   if(part->Edot[np]>1e28) {
                                                           Nsup_x+=1;
                                                   }
                                                   if(part->Edot[np]>1e31) {
                                                           Nsup2_x+=1;
                                                   }
                                                   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
                                }
				}

				else if(Ng==1 && Nr==0 && Nx==1){ //gamma and x-rays

                                        if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_gam[np]=0;part->detec_x[np]=0;part->detec_gam_x[np]=0;}
                                        else{
                                           count_gamma_x++;part->detec[np]=1;part->detec_gam[np]=0;part->detec_x[np]=0;part->detec_gam_x[np]=1;

                                                   if(part->Edot[np]>1e28) {
                                                           Nsup_gamma_x+=1;
                                                   }
                                                   if(part->Edot[np]>1e31) {
                                                           Nsup2_gamma_x+=1;
                                                   }
                                                   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
                                }
				}

				else if(Ng==0 && Nr==1 && Nx==1){ //radio and x-rays

                                        if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_rad[np]=0;part->detec_rad_x[np]=0;part->detec_x[np]=0;}
                                        else{
                                           count_rad_x++;part->detec[np]=1;part->detec_x[np]=0;part->detec_rad[np]=0;part->detec_rad_x[np]=1;

                                                   if(part->Edot[np]>1e28) {
                                                           Nsup_rad_x+=1;
                                                   }
                                                   if(part->Edot[np]>1e31) {
                                                           Nsup2_rad_x+=1;
                                                   }
                                                   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
                                }
				}

				else if(Ng==1 && Nr==1 && Nx==1){ //gamma radio and x-rays

                                        if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_gam[np]=0;part->detec_x[np]=0;part->detec_rad[np]=0;}
                                        else{
                                           count_rad_gam_x++;part->detec[np]=1;part->detec_gam[np]=0;part->detec_x[np]=0;part->detec_rad[np]=0;part->detec_rad_gam_x[np]=1;

                                                   if(part->Edot[np]>1e28) {
                                                           Nsup_rad_gam_x+=1;
                                                   }
                                                   if(part->Edot[np]>1e31) {
                                                           Nsup2_rad_gam_x+=1;
                                                   }
                                                   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
                                }
				}*/

				if((Ng==1 || Nr==1)) {
			        		if(part->Edot[np]>1e28) { 
							Nsup++;
						}
			        		if(part->Edot[np]>1e31) { 
							Nsup2++;
						}	
				}
			Nr=0;
			Ng=0;
			detec=0;
			       if(part->detec_rad[np]==1 && isnan(part->period[np])!=true && isnan(part->Pdot[np])!=true){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|1|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->alpha0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np],part->t_acc[np],part->spin_up_or_not[np],part->M_NS_i[np],part->M_NS_f[np],part->Mc_i[np],part->Mc_f[np],part->a_bin[np],part->ecc_bin[np],part->comp_type[np],part->comp_rem_type[np],part->ellipticity[np]);
				  fprintf(file_taud,"%e\n",part->taud_evol[np]);
				  if (current_detec_fast==1) fprintf(file_wr,"%e\n",part->w_r_fast[np]*180.0/M_PI);
				  else if (current_detec_pmps==1) fprintf(file_wr,"%e\n",part->w_r_pmps[np]*180.0/M_PI);
				  //else if (current_detec_htru==1) fprintf(file_wr,"%e\n",part->w_r_htru[np]*180.0/M_PI);
				  fprintf(gamma_peak_sep,"%e\n",1e55);
				  fprintf(beta_file,"%e\n",part->beta[np]);
				  fprintf(betainit_file,"%e\n",part->beta0[np]);
				  fprintf(save_tempo,"%e %e\n",xi,rho);
				  fprintf(pos_init,"%e %e %e\n",part->x0[np],part->y0[np],part->z0[np]);
				  fprintf(file_wint,"%e\n",part->w_int[np]*180.0/M_PI);
                        }

                               else if(part->detec_gam[np]==1 && isnan(part->period[np])!=true && isnan(part->Pdot[np])!=true){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|2|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->alpha0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np],part->t_acc[np],part->spin_up_or_not[np],part->M_NS_i[np],part->M_NS_f[np],part->Mc_i[np],part->Mc_f[np],part->a_bin[np],part->ecc_bin[np],part->comp_type[np],part->comp_rem_type[np],part->ellipticity[np]);
                                  fprintf(file_taud,"%e\n",part->taud_evol[np]);
				  fprintf(file_wr,"%e\n",part->w_r_fast[np]*180.0/M_PI);
				  fprintf(gamma_peak_sep,"%e\n",part->delta[np]);
				  fprintf(beta_file,"%e\n",part->beta[np]);
				  fprintf(betainit_file,"%e\n",part->beta0[np]);
				  fprintf(save_tempo,"%e %e\n",xi,rho);
				  fprintf(pos_init,"%e %e %e\n",part->x0[np],part->y0[np],part->z0[np]);
				  fprintf(file_wint,"%e\n",part->w_int[np]*180.0/M_PI);
			}

                               else if(part->detec_rg[np]==1 && isnan(part->period[np])!=true && isnan(part->Pdot[np])!=true){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|3|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->alpha0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np],part->t_acc[np],part->spin_up_or_not[np],part->M_NS_i[np],part->M_NS_f[np],part->Mc_i[np],part->Mc_f[np],part->a_bin[np],part->ecc_bin[np],part->comp_type[np],part->comp_rem_type[np],part->ellipticity[np]);
                                  fprintf(file_taud,"%e\n",part->taud_evol[np]);
				  fprintf(gamma_peak_sep,"%e\n",part->delta[np]);
				  if (current_detec_fast==1) fprintf(file_wr,"%e\n",part->w_r_fast[np]*180.0/M_PI);
                                  else if (current_detec_pmps==1) fprintf(file_wr,"%e\n",part->w_r_pmps[np]*180.0/M_PI);
				  //else if (current_detec_htru==1) fprintf(file_wr,"%e\n",part->w_r_htru[np]*180.0/M_PI);
				  fprintf(beta_file,"%e\n",part->beta[np]);
				  fprintf(betainit_file,"%e\n",part->beta0[np]);
				  fprintf(save_tempo,"%e %e\n",xi,rho);
				  fprintf(pos_init,"%e %e %e\n",part->x0[np],part->y0[np],part->z0[np]);
				  fprintf(file_wint,"%e\n",part->w_int[np]*180.0/M_PI);
			}

			       /*else if(part->detec_x[np]==1 && part->accretion[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|4|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }

			       else if(part->detec_gam_x[np]==1 && part->accretion[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|5|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }

			       else if(part->detec_rad_x[np]==1 && part->accretion[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|6|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }


			       else if(part->detec_rad_gam_x[np]==1 && part->accretion[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|7|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }*/

		}

 
			//Ntot = count_radio+count_gamma+count_radio_gamma+count_x+count_gamma_x+count_rad_x+count_rad_gam_x;
			Ntot = count_radio+count_gamma+count_radio_gamma;


			printf("# Number of pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup);
			printf("# Number of radio-only pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_radio);
			printf("# Number of gamma-only pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_gamma);
			printf("# Number of radio+gamma pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_radio_gamma);
			/*printf("# Number of x-rays only pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_x);
			printf("# Number of gamma-x-rays pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_gamma_x);
			printf("# Number of radio-x-rays pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_rad_x);
			printf("# Number of gamma, radio and x-rays pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_rad_gam_x);*/
			printf(" \n");
			printf("# Number of pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2);
			printf("# Number of radio-only pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_radio);
			printf("# Number of gamma-only pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_gamma);
			printf("# Number of radio+gamma pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_radio_gamma);
			/*printf("# Number of x-rays only pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_x);
                        printf("# Number of gamma-x-rays pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_gamma_x);
                        printf("# Number of radio-x-rays pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_rad_x);
                        printf("# Number of gamma, radio and x-rays pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_rad_gam_x);*/
			printf(" \n");
			printf("# Total number of detected pulsars  %ld \n",Ntot);
			printf("# Total number of radio-only pulsars detected %ld \n",count_radio);
			printf("# Total number of gamma-only pulsars detected %ld \n",count_gamma);
			printf("# Total number of radio+gamma pulsars detected %ld \n",count_radio_gamma);
			/*printf("# Total number of x-rays-only pulsars detected %ld \n",count_x);
			printf("# Total number of gamma, x-rays pulsars detected %ld \n",count_gamma_x);
			printf("# Total number of radio-x-rays pulsars detected %ld \n",count_rad_x);
			printf("# Total number of gamma,radio and x-rays pulsars detected %ld \n",count_rad_gam_x);*/
			printf("# Number of radio pulsars beaming to us %ld \n",Nbeam);
			printf("# Number of radio pulsars detected with PMPS: %d\n",detec_pmps);
			printf("# Number of radio pulsars detected with fast: %d\n",detec_fast);
			printf("# Number of radio pulsars detected with fast: %d\n",detec_htru);
			printf("# Number of pulsars geometrically detectable with flux maximum five times less than required %d\n",count_possible_gdetec);
			printf("# Number of pulsars detected both in radio and gamma, with Fr<30 µJy : %ld\n",count_rg_lowF);
			printf("# Number of pulsars who died in the valley : %d\n",count_dead); 
			printf("# Number of detected gamma pulsars with DM>2.5 pc.cm^-3 : %d\n",count_dm_sup);
			printf(" \n");
			fclose(file_data);
			fclose(file_taud);
			fclose(file_wr);
			fclose(beta_file);
			fclose(betainit_file);
			fclose(gamma_peak_sep);
			fclose(save_tempo);
			fclose(pos_init);
			fclose(file_wint);
			//fclose(check_val2);

return(0);
}

/* /!\  units should be changed to SI units */
int radio_flux(void *params){ // returns Smin and radio flux Fr in mJy 
    	struct func_params *part= (struct func_params*)params;
   	double Fj; //scatter term
	long np=0;
	double wtilde_fast,wtilde_pmps,wtilde_htru;
	double S_0=0.05; //Radiometer equation for long period pulsars, Johnston et al. (2020)
	double C_thres,T_rec,n_pol,Gain,B_mhz,int_t;
	double pedot=-1.0;
	int count_comp_lorimer25=0;
	n_pol=2;

           	for (np=0;np<part->Npulsars;np++){ 
		
			part->np=np;	
			wtilde_fast=part->w_r_fast[np]*part->period[np]/(2*M_PI); // converting the width from radians to s if predections, fast ->SKAlow 
			wtilde_pmps=part->w_r_pmps[np]*part->period[np]/(2*M_PI); // converting the width from radians to s if predections, pmps ->SKAmid
			wtilde_htru=part->w_r_htru[np]*part->period[np]/(2*M_PI);

			//SKA-Low
			C_thres=9;T_rec=46;Gain=23.0;int_t=600;B_mhz=100e6;
                        S_0=(C_thres*(T_rec+(part->temp[np]*pow(408.0/300.0,2.6)))/(Gain*sqrt(n_pol*B_mhz*int_t)))*1e3; //in mJy
                        if (part->period[np]-wtilde_fast>0){
                                part->Smin_SKAlow[np]=S_0*sqrt(wtilde_fast/((part->period[np])-wtilde_fast)); //in mJy
                        }
                        else if (part->period[np]-wtilde_fast<=0){
                                part->Smin_SKAlow[np]=1e55;
                        }

			//SKA-mid
                        C_thres=9;T_rec=30;Gain=15.0;int_t=2100;B_mhz=300e6;
                        S_0=(C_thres*(T_rec+(part->temp[np]*pow(408.0/1400.0,2.6)))/(Gain*sqrt(n_pol*B_mhz*int_t)))*1e3; //in mJy
                        if (part->period[np]-wtilde_pmps>0){
                                part->Smin_SKAmid[np]=S_0*sqrt(wtilde_pmps/((part->period[np])-wtilde_pmps)); //in mJy
                        }
                        else if (part->period[np]-wtilde_pmps<=0){
                                part->Smin_SKAmid[np]=1e55;
                        }

			//FAST
			C_thres=9;T_rec=25;Gain=16.;int_t=300;B_mhz=450e6;
			S_0=(C_thres*(T_rec+(part->temp[np]*pow(408.0/1374.0,2.6)))/(Gain*sqrt(n_pol*B_mhz*int_t)))*1e3; //in mJy 
			if (part->period[np]-wtilde_fast>0){
                                part->Smin_fast[np]=S_0*sqrt(wtilde_fast/((part->period[np])-wtilde_fast)); //in mJy
                        }
                        else if (part->period[np]-wtilde_fast<=0){
                                part->Smin_fast[np]=1e55;
                        }

			//PMPS survey
			C_thres=9;T_rec=21;Gain=0.735;int_t=2100;B_mhz=288e6;
			S_0=(C_thres*(T_rec+(part->temp[np]*pow(408.0/1374.0,2.6)))/(Gain*sqrt(n_pol*B_mhz*int_t)))*1e3; //in mJy

			if (part->period[np]-wtilde_pmps>0){
		        	part->Smin_pmps[np]=S_0*sqrt(wtilde_pmps/((part->period[np])-wtilde_pmps)); //in mJy	
			}
			else if (part->period[np]-wtilde_pmps<=0){
				part->Smin_pmps[np]=1e55;
			}
			//HTRU-mid
			C_thres=9;T_rec=23;Gain=0.735;int_t=540;B_mhz=340e6;
                        S_0=(C_thres*(T_rec+(part->temp[np]*pow(408.0/1352.0,2.6)))/(Gain*sqrt(n_pol*B_mhz*int_t)))*1e3; //in mJy
                        if (part->period[np]-wtilde_htru>0){
                                part->Smin_htru[np]=S_0*sqrt(wtilde_htru/((part->period[np])-wtilde_htru)); //in mJy
                        }
                        else if (part->period[np]-wtilde_htru<=0){
                                part->Smin_htru[np]=1e55;
                        }
			Fj             = fabs(gsl_ran_gaussian_ziggurat(part->r,0.2));
                        //while (pedot<0) {pedot  = 0.1+gsl_ran_gaussian_ziggurat(part->r, 0.07);}
                        pedot=0.125;
                        part->Fr[np]   =(9.0/sq(part->dist[np]))*pow(part->Edot[np]/1e29,pedot)*pow(10,Fj); // (From Johnston's paper) in mJy
                        if (part->Fr[np]>=0.1) count_comp_lorimer25+=1;

		}
		printf("Number of pulsars above 0.1 mJy : %d\n",count_comp_lorimer25);

 



	return(0);

}

int get_fomega(void *params){

	struct func_params *part= (struct func_params*)params;
	FILE *fomega_file=NULL;
	fomega_file=fopen("facteur_omega.dat","r");
	for(int i=0;i<91;i++){
		for(int j=0;j<91;j++){
			if (fscanf(fomega_file, "%lf", &part->fomega[i][j]) != 1) {
    				fprintf(stderr, "Failed read of values [%d][%d]\n", i, j);
    				fclose(fomega_file);
    				return EXIT_FAILURE;
			}

		}
	}
	fclose(fomega_file);
	return 0;

}


int gamma_flux(void *params){

    	struct func_params *part= (struct func_params*)params;
   	double cg; //scatter term
	cg=1.0;
	double lgamma;
	const double kpc2m=3.0856775807e19; //kpc to m
	//fac=0.3;
	long np=0;
	double pcst,pb;
	double pedot=-1.0;
	double alpha,xi;
	int alpha_int,xi_int;

           	for (np=0;np<part->Npulsars;np++){ 
		
			part->np=np;	

                        //lgamma = 1e15*pow(part->B[np]*1e4,0.11)*pow(part->Edot[np]*1e7,0.51); 
			//pcst=26.15+gsl_ran_gaussian_ziggurat(part->r, 2.6);
			//pb=0.11+gsl_ran_gaussian_ziggurat(part->r, 0.05);
			//while (pedot<0) {pedot=0.51+gsl_ran_gaussian_ziggurat(part->r, 0.09);}
			pcst=26.4;pb=0.11;pedot=0.51;
			//pb=0.11;pedot=0.51;
                        lgamma = pow(10,pcst)*pow(part->B[np]/1e8,pb)*pow(part->Edot[np]/1e26,pedot); //(W) from Kalapotharakos et al 2019
			part->Lgamma[np]=lgamma;
			if (part->alpha[np]<=M_PI/2) alpha=part->alpha[np];
                        else if (part->alpha[np]>M_PI/2) alpha=M_PI-part->alpha[np];
                        if (part->xi[np]<=M_PI/2) xi=part->xi[np]; // just to make it easier to read
                        else if (part->xi[np]>M_PI/2) xi=M_PI-part->xi[np];
			//printf("xi : %e\n",180*xi/M_PI);
		    
			alpha=(180.0/M_PI)*alpha;alpha_int=(int)round(alpha);
			xi=(180.0/M_PI)*xi;xi_int=(int)round(xi);
			cg=part->fomega[alpha_int][xi_int];
			//printf("%e,%d,%d,%ld\n",cg,alpha_int,xi_int,np);
			/*if (alpha<=-xi+70) {cg=1.9;}
			else if (alpha<=-xi+72.5) {cg=1.713;}
			else if (alpha<=-xi+75) {cg=1.527;}
			else if (alpha<=-xi+76) {cg=1.340;}
			else if (alpha<=-xi+77 || (alpha>=-(45.0/27.5)*xi+192.27)) {cg=1.153;}
			else if (alpha<=-xi+78 || (alpha>=-xi+130)) {cg=0.967;}
			else if (alpha<=-xi+79 || (alpha>=-(67.0/64.0)*xi+117.22)) {cg=0.780;}
			else if (alpha<=-xi+81 || (alpha>=-(75.0/70.0)*xi+111.43)) {cg=0.593;}
			else if ((alpha<=8 && xi >= 85) || (alpha>=85 && xi <= 7.5)) {cg=0.22;}
			else {cg=0.407;}*/
	        	//if(alpha<-xi+0.6109) cg=1.9; //the correction factor cg (or f_omega)
		   	//	else cg=1.;

      			part->Fg[np]=(1.0/(4*M_PI*cg*sq(part->dist[np]*kpc2m)))*lgamma; // gamma flux in W.m-2
			//if (part->accretion[np]==1) printf("Fg=%e,B=%e,Edot=%e\n",part->Fg[np],part->B[np],part->Edot[np]);
			//printf("Fg %e np %ld dist %e \n",part->Fg[np],np,part->dist[np]);		
                        //printf("lgamma %e Fg %e B %e Edot %e \n",lgamma,part->Fg[np],part->B[np],part->Edot[np]);	

      			//part->Fg[np]=fac*(part->Edot[np]*1e7/(sq(part->dist[np]*KPC2CM)))*sqrt(1e33/(part->Edot[np]*1e7)); // Petri 2011a 
 				
		} 
		printf("Fg computed for everyone\n");

return(0);


}
