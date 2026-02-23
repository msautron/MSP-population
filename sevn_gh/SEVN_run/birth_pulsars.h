#ifndef BIRTH_PULSARS_H
#define BIRTH_PULSARS_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

int count_lines(char *filename);
double initial_mass_function(double M1);
double initial_a_function(double a);
double initial_e_function(double e);
double initial_q_function(double q);
void birth(gsl_rng *r,long index_sim);
double sample_q(double qmin, double qmax, gsl_rng *r);

#endif
