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
	char cmd_bis[500];
	char cmd2[500];
	char cmd3[500];
	char cmd_open[500];
	double p_taud=gsl_rng_uniform(r);
	double taud;
	if (p_taud<0.77) taud=0.05;
	else if (p_taud<0.885 && p_taud>=0.77) taud=1.4;
	else if (p_taud>=0.885) taud=2.3;
	sprintf(cmd,"bash /home/matteo.sautron/sevn_gh/SEVN_run/run.sh %ld %lf",index_sim,taud);
	sprintf(cmd_bis,"bash /home/matteo.sautron/sevn_gh/SEVN_run/run2.sh %ld %lf",index_sim,taud);
	sprintf(cmd2,"python3 /home/matteo.sautron/sevn_gh/SEVN_run/get_info_binary.py %ld",index_sim);
	sprintf(cmd3,"/home/matteo.sautron/sevn_gh/SEVN_run/data_pop_ms%ld.txt",index_sim);
	sprintf(cmd_open,"/home/matteo.sautron/sevn_gh/SEVN_run/listBin%ld.dat",index_sim);
        printf("Npulsars = %ld\n",Npulsars);
        for(long i=0;i<Npulsars;i++){
		FILE *listBin=NULL;
                listBin=fopen(cmd_open,"r");
                char trash1[64];char trash2[64];double trash3;double trash4;int trash5;double target;
                fscanf(listBin,"%63s %lf %lf %63s %d %lf",trash1,&trash3,&trash4,trash2,&trash5,&target);
                fclose(listBin);
                printf("Mass of the companion : %e\n",target);
		birth(r,index_sim);
		if (target>=2.0){
                	system(cmd);
			printf("run.sh was run\n");
		}
		else {
			system(cmd_bis);
			printf("run2.sh was run\n");
		}
        }
        return 0;
}
