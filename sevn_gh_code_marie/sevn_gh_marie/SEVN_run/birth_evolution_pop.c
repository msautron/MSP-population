#define _GNU_SOURCE
#include <stdio.h>   // Pour perror, printf
#include <stdlib.h>  // Pour strtod, NULL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include "birth_pulsars.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <sched.h>

int main(int argc, char **argv){

	//unsigned int sleep_time = getpid() % 100;
	//sleep(sleep_time);
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
        long Npulsars = (long) strtod(argv[1],(char **)NULL);
	long index_sim = (long) strtod(argv[2],(char **)NULL);
	char cmd[500];
	//char cmd1[500];
	char cmd2[500];
	char cmd3[500];
	char cmd4[500];
	char cmd5[500];
	//char cmd6[500];
	//char cmd7[500];
	double p_taud=gsl_rng_uniform(r);
	double taud;
	//system("rm /home/mepietrin/sevn_gh/SEVN_run/data_period_pulsars*.txt");
	//system("rm /home/mepietrin/sevn_gh/SEVN_run/data_alpha_pulsars.txt");
	//system("rm /home/mepietrin/sevn_gh/SEVN_run/data_beta_pulsars.txt");
	//system("rm /home/mepietrin/sevn_gh/SEVN_run/data_p_dot_pulsars*.txt");
	//sprintf(cmd1,"rm /home/mepietrin/sevn_gh/SEVN_run/data_WD_type%ld.txt",index_sim);
	//system(cmd1);
	//FILE *f=fopen("/home/mepietrin/sevn_gh/SEVN_run/data_period_pulsars.txt","r");
	if (p_taud<0.77) taud=0.05;
	else if (p_taud<0.885 && p_taud>=0.77) taud=1.4;
	else if (p_taud>=0.885) taud=2.3;
	sprintf(cmd,"bash /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/run.sh %ld %lf",index_sim,taud);
	sprintf(cmd2,"python3 /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/get_info_binary%ld.py",index_sim);
	sprintf(cmd3,"/home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_pop_ms%ld.txt",index_sim);
        sprintf(cmd4,"tail -n 1 /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_pop_ms%ld.txt >> /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_end_pulsars%ld.txt",index_sim,index_sim);
	//sprintf(cmd4,"tail -n 1 /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_pop_ms%ld.txt >> /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_end_pulsars%ld.txt",index_sim,index_sim);
	//sprintf(cmd6,"paste -d ' ' /home/mepietrin/sevn_gh/SEVN_run/data_end_pulsars%ld.txt /home/mepietrin/sevn_gh/SEVN_run/data_mass_et_type%ld.txt >> /home/mepietrin/sevn_gh/SEVN_run/data_end_pulsars%ld.txt",index_sim,index_sim,index_sim);
	sprintf(cmd5,"rm /home/matteo.sautron/sevn_gh_code_marie/sevn_gh_marie/SEVN_run/data_pop_ms%ld.txt",index_sim);
	//sprintf(cmd6,"rm /home/mepietrin/sevn_gh/SEVN_run/data_end_mass_type%ld.txt",index_sim);
	//sprintf(cmd7,"rm /home/mepietrin/sevn_gh/SEVN_run/data_mass_et_type%ld.txt",index_sim);
	printf("Npulsars = %ld\n",Npulsars);
       	for(long i=0;i<Npulsars;i++){
		birth(r,index_sim);
                system(cmd);
		system(cmd4);
		//system(cmd6);
		system(cmd5);
		//system(cmd7);
        }
	//system(cmd4);
	return 0;
	
}
