#define _GNU_SOURCE
#include "birth_pulsars.h"
#include<gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<time.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>  
#include <sys/time.h>
#include <sched.h>

int count_lines(char *filename) {
    char command[256];
    snprintf(command, sizeof(command), "wc -l < %s 2>/dev/null", filename);
    FILE *pipe = popen(command, "r");
    if (!pipe) return -1; 
    int lines = 0;
    fscanf(pipe, "%d", &lines);
    pclose(pipe);
    return lines;
}

double initial_mass_function(double M1){
        return pow(M1,-2.3);
}

double initial_a_function(double a){
        double R_SUN=696340e3;
        double a0=10*R_SUN;
        double alpha_sep=0.07;
        double m=1.2;
        if (a<=a0) return alpha_sep*pow(a,m-1.0)*pow(1/a0,m);
        else if (a>a0) return alpha_sep/a;
}

double initial_e_function(double e){
	return pow(e,-0.42);
}

double initial_q_function(double q){
	return pow(q,-0.1);
}

void birth(gsl_rng *r,long index_sim){

	FILE *birth_file=NULL;
	char cmd[500];
	sprintf(cmd,"/home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/listBin%ld.dat",index_sim);
	birth_file=fopen(cmd,"w+");
	double M_NS;
	double M_WD;
	//double Z_M_WD=;
	double Z_M_NS=0.02;
	double spin_M_NS=3.52e-4;
	double M1;
	double M2;
	double Z2;
	double spin_2=2;
	int MS_lp=0;
	double a;
	double e;
	double age_pulsar;
	double step;
	double MSUN=1.98847e30;
	bool sample;
	double pdf_val;
	double comp_val;
	double alpha_sep=0.07;
	double R_SUN=696340e3;
	double a0=10*R_SUN;
	double log_a_binary;
	double q;
	//int type_WD;
	//int HEWD=0;
	//int COWD=1;
	//int ONEWD=2;
	//FILE *data_mass_et_type;
	//char cmd1[500];
	//sprintf(cmd1,"/home/mepietrin/sevn_gh/SEVN_run/data_mass_et_type%ld.txt",index_sim);
        //data_mass_et_type=fopen(cmd1,"w+");


	M_NS=(1.33+gsl_ran_gaussian_ziggurat(r,0.09));
	while (sample==false){
                M1=(8-0.01)*gsl_rng_uniform(r)+0.01;
                pdf_val=initial_mass_function(M1);
                comp_val=pow(8,-2.3)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
        sample=false;
	/*while (sample==false){
                q=gsl_rng_uniform(r);
                pdf_val=initial_q_function(q);
                comp_val=pow(1e-4,-0.1)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }*/
	//type_WD=rand()%3;
//	printf("le type est :%lf",type_WD);

	/*if(type_WD==0){
		M_WD=(0.4+gsl_ran_gaussian_ziggurat(r,0.140));
	}
	if(type_WD==1){
                M_WD=(1.0+gsl_ran_gaussian_ziggurat(r,0.140));

        }
	if(type_WD==2){*/
        M_WD=(1.2+gsl_ran_gaussian_ziggurat(r,0.140));

        


	//fprintf(data_type_WD,"\n");
	//fprintf(data_mass_et_type, "%e %d\n",M_WD,type_WD);

        while (sample==false){
                q=gsl_rng_uniform(r);
                pdf_val=initial_q_function(q);
                comp_val=pow(1e-4,-0.1)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
        sample=false;
	
	M2=q*M1;
	Z2=(log(0.02)-log(0.0002))*gsl_rng_uniform(r)+log(0.0002);
	Z2=pow(10.0,Z2);
	while (spin_2>1 || spin_2 <0){
		spin_2=gsl_rng_uniform(r);
		//spin_2=(5e-4+gsl_ran_gaussian_ziggurat(r,3e-4));
		//spin_2=(0.4+gsl_ran_gaussian_ziggurat(r,0.3));
	}
	while (sample==false){
        	log_a_binary=(log(5.75e6*R_SUN)-log(R_SUN))*gsl_rng_uniform(r)+log(R_SUN);
                a=exp(log_a_binary);
                pdf_val=initial_a_function(a);
                comp_val=(alpha_sep/a0)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
        sample=false;
	a=a/R_SUN;
	while (sample==false){
                e=gsl_rng_uniform(r);
                pdf_val=initial_e_function(e);
                comp_val=pow(1e-10,-0.42)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
	sample=false;
	//age_pulsar=8e8*1e-6;
	age_pulsar=((13.9e9-5e7)*gsl_rng_uniform(r)+5e7)*1e-6; //age in Myr
	step=age_pulsar/5000.0;
	MS_lp=100*gsl_rng_uniform(r);
	/*if(type_WD==0){
	fprintf(birth_file,"%eHEWD %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:1 %e %e %e %e\n",M_WD,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,a,e,age_pulsar,step);
	}

	if(type_WD==1){
        	fprintf(birth_file,"%eCOWD %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:1 %e %e %e %e\n",M_WD,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,a,e,age_pulsar,step);
        }*/

	//if(type_WD==2){
        fprintf(birth_file,"%eONEWD %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:1 %e %e %e %e\n",M_WD,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,a,e,age_pulsar,step);
        


	//fprintf(birth_file,"%eNS %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:1 %e %e %e %e\n",M_NS,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,a,e,age_pulsar,step);
	//fprintf(birth_file,"%eWD %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:1 %e %e %e %e\n",M_WD,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,a,e,age_pulsar,step);

	fclose(birth_file);
	//fclose(data_mass_et_type);

/*//ajout    
     FILE *info_birth;
     char cmd[500];
     char eonf[5];
     long end_name_file;
     end_name_file=s->get_svpar_num("end_of_name_file");
     sprintf(eonf,"%ld",end_name_file);
     sprintf(cmd,"/home/mepietrin/sevn_gh/SEVN_run/listBin%s.dat",eonf);
     info_birth=fopen(cmd,"r");
     double trash;char trash_str[100];
 
     fscanf(info_birth,"%le%s %le %le %s %le %le %le %le %s %s %le %le %le %le",&trash,trash_str,&trash,&trash,trash_str,&trash,&trash,&trash,&trash,trash_str,trash_str,&trash,&trash,&real_age,&trash);
 
     if (age(s)==real_age) {
         FILE *data_wd_mass;
         char cmd1[500];
         char eonf[5];
         long end_name_file;
         end_name_file=s->get_svpar_num("end_of_name_file");
         sprintf(eonf,"%ld",end_name_file);
         sprintf(cmd1, "/home/mepietrin/sevn_gh/SEVN_run/data_wd_mass%s.txt",eonf);
         data_wd_mass=fopen(cmd,"a+");
         fprintf(data_wd_mass,"%e\n",Mass);

     }
*/
}

/*int main(){
	
	printf("got here");
	gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
        struct timespec ts_real, ts_mono;
        clock_gettime(CLOCK_REALTIME, &ts_real);
        clock_gettime(CLOCK_MONOTONIC, &ts_mono);

        // Get a few bytes from /dev/urandom
        int urandom_fd = open("/dev/urandom", O_RDONLY);
        unsigned int urandom_seed;
        if (urandom_fd != -1) {
        read(urandom_fd, &urandom_seed, sizeof(urandom_seed));
        close(urandom_fd);
        } else {
        urandom_seed = 0;  // Fallback if /dev/urandom is not available
        }

        // Combine multiple sources of entropy
        unsigned long int seed = ts_real.tv_sec * 1000000000L + ts_real.tv_nsec +
                                ts_mono.tv_sec * 1000000000L + ts_mono.tv_nsec +
                                getpid() + sched_getcpu() + urandom_seed;

        gsl_rng_set(r, seed);
	birth(r);
	return 0;
}*/
