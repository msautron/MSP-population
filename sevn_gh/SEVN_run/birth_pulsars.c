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

double initial_mass_function(double M1){ //Kroupa (2001)
        return pow(M1,-2.3);
}

double initial_a_function(double a){ //Han et al. (2020) 
        double R_SUN=696340e3;
        double a0=10*R_SUN;
        double alpha_sep=0.07;
        double m=1.2;
        if (a<=a0) return alpha_sep*pow(a,m-1.0)*pow(1/a0,m);
        else if (a>a0) return alpha_sep/a;
}

double initial_e_function(double e){ //Sgalletta et al. (2023)
	return pow(e,-0.42);
}

double initial_q_function(double q){ //Sana et al. (2012)
	return pow(q,-0.1);
}

double sample_q(double qmin, double qmax, gsl_rng *r)
{
    double u = gsl_rng_uniform(r);
    double a = pow(qmin, 0.9);
    double b = pow(qmax, 0.9);
    double q = pow(a + u*(b - a), 1.0/0.9);
    return q;
}


void birth(gsl_rng *r,long index_sim){

	FILE *birth_file=NULL;
	char cmd[500];
	sprintf(cmd,"/home/matteo.sautron/sevn_gh/SEVN_run/listBin%ld.dat",index_sim);
	birth_file=fopen(cmd,"w+");
	double M_NS;
	double Z_M_NS=0.02;
	double spin_M_NS=3.52e-4;
	double M1;
	double M2;
	double Z2;
	double spin_2=2;
	int MS_lp=0;
	int MS_ls=1;
	double a;
	double e;
	double age_pulsar;
	double step;
	double MSUN=1.98847e30;
	double G=6.67430e-11; //Gravitational constant
	bool sample;
	double pdf_val;
	double comp_val;
	double alpha_sep=0.07;
	double R_SUN=696340e3;
	double a0=10*R_SUN;
	double log_a_binary;
	double q;
	double P_orb;
	M_NS=(1.33+gsl_ran_gaussian_ziggurat(r,0.09));
	while (sample==false){
                M1=(25-8)*gsl_rng_uniform(r)+8;
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
        }
        sample=false;*/
	q=sample_q(1.0/M1,1.0,r);
	M2=q*M1;
	//M2=1.5;//(3-1)*gsl_rng_uniform(r)+1;
	Z2=pow(10,(log10(0.02)-log10(0.0002))*gsl_rng_uniform(r)+log10(0.0002));
	while (spin_2>1 || spin_2 <0){
		spin_2=gsl_rng_uniform(r);
		//spin_2=(5e-4+gsl_ran_gaussian_ziggurat(r,3e-4));
		//spin_2=(0.4+gsl_ran_gaussian_ziggurat(r,0.3));
	}
	//Old consideration, as if both stars were main sequence stars
	while (sample==false){
        	log_a_binary=(log(5.75e6*R_SUN)-log(R_SUN))*gsl_rng_uniform(r)+log(R_SUN);
                a=exp(log_a_binary);
                pdf_val=initial_a_function(a);
                comp_val=(alpha_sep/a0)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
        sample=false;
	//log_a_binary=(log(5.75e6*R_SUN)-log(R_SUN))*gsl_rng_uniform(r)+log(R_SUN);
	//a=exp(log_a_binary);
	//a=((12-9)*gsl_rng_uniform(r)+9)*R_SUN;
	//Inspired by observations of El-Badry et al. (2024)
	/*P_orb=((1000-100)*gsl_rng_uniform(r)+100)*24*60*60; //P_orb in s
	a=pow((G*(M2+M_NS)*MSUN*pow(P_orb,2.0))/(4*pow(M_PI,2.0)),1.0/3.0);*/ //a in meter
	a=a/R_SUN;
	//Old consideration as if both stars were main sequence
	while (sample==false){
                e=gsl_rng_uniform(r);
                pdf_val=initial_e_function(e);
                comp_val=pow(1e-10,-0.42)*gsl_rng_uniform(r);
                if (comp_val<=pdf_val) {sample=true;}
        }
	sample=false;
	//e=gsl_rng_uniform(r);
	//e=(1-0.1)*gsl_rng_uniform(r)+0.1; //Inspired by observations of El-Badry et al. (2024)
	age_pulsar=((13e9-5e7)*gsl_rng_uniform(r)+5e7)*1e-6; //age in Myr
	//age_pulsar=pow(10,(log10(1.39)+10-log10(3)-8)*gsl_rng_uniform(r)+(log10(3)+8))*1e-6; //age in Myr
	step=age_pulsar/50000.0;
	MS_lp=(int)(gsl_rng_uniform(r)*100);
	//MS_ls=(int)(1.0+gsl_rng_uniform(r)*5);
	fprintf(birth_file,"%eNS %e %e delayed_gauNS 0 %e %e %e delayed_gauNS %%%d:%d %e %e %e %e\n",M_NS,Z_M_NS,spin_M_NS,M2,Z2,spin_2,MS_lp,MS_ls,a,e,age_pulsar,step);
	fclose(birth_file);
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
