#ifndef _MACRO_DEFINED
#define _MACRO_DEFINED
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include <stdbool.h>
#define cube(a) ((a)*(a)*(a))
#define maximum(a,b) ((a)>(b)?(a):(b))
#define minimum(a,b) ((a)<(b)?(a):(b))
#define N_MAX   10000000 // maximum number of pulsars
#define sq(a) ((a)*(a))

struct func_params{

	double *t_acc; //Stores duration of accretion
	double *M_NS_f; //Stores final value of the NS mass
	double *M_NS_i; //Stores initial value of the NS mass
	double *spin_up_or_not; //Stores 1 if accretion happened, 0 if not
	double *Mc_f; // Stores the final mass of the companion in solar mass
	double *Mc_i; // Stores the initial mass of the companion in solar mass
	double *a_bin; //Stores the semi major axis of the binary
	double *ecc_bin; // Stores the eccentricity of the binary
	double *ellipticity; //Stores the ellipticity of the NS
	double *comp_type; //Stores the type of the companion if it is not a remnant
	double *comp_rem_type; //Stores the type of the companion if it is a remnant
	double *DM; //Stores dispersion measure
	double *temp; //Stores sky temperature at (gl,gb) position
	double *n_omega_x;
	double *n_omega_y;
	double *n_omega_z;
	double *PA;
	double *taud_evol; //Stores the tau_d chosen randomly for evolution in SEVN 
        double R; //Radius of the neutron star
        double *alpha;//angle between the magnetic field and the rotation axis
        double *Binit; //table which returns the initial B
        double *B; //table which returns the initial B
        double *Pinit; //same for P
        double *age_pulsar; // stores the age of the pulsar (in s)
	double Bfield_const;
	long Bfield_var;
        double *x; //kpc  x coordinate in the Galctocentric frame
        double *y;//idem 
        double *z; // kpc height to the galactic plane
	double *gb; // galactic latitude
	double *gl; // galactic longitude
        double *z0; //kpc initial z0 in the Galactocentric frame
        double *x0; // idem 
        double *y0; //idem
	double *vx0; //table which stores the initial velocities on the x absciss, km/s
	double *vy0; //table which stores the initial velocities on the y absciss, km/s
	double *vz0; //table which stores the initial velocities on the z absciss, km/s
        double *vx;  // velocities on the x absciss, km/s
	double *vy;  // velocities on the y absciss, km/s
	double *vz;  // velocities on the z absciss, km/s
        double *err_rel_g; // relative error on the energy for the integration of the equation of movement 
        long   np; //pulsar number= index of tables that store the pulsars parameters) 
        double sigma_v; // sigma 1D
        double *period; //stores the actual period
        long    Npulsars; //total number of pulsar
	double Rexp; //in kpc parameter for the R distribution
	long v_young; //sigma 1D for the young pulsarvelocity
	long v_old; // old
        long *detec; // indicates if the pulsar is considered as detected
	long *detec_rad; // indicates if the pulsar is detected in radio
	long *detec_gam; //indicates if the pulsar is detected in gamma
	long *detec_rg; // indicates if the pulsar is detected in radio and gamma
        double zexp; // in kpc
	double *Smin_htru; //Minimum flux detectable by HTRU-mid
	double *Smin_pmps; //Minimum flux detectable by GBNCC
	double *Smin_fast; //Minimum flux detectable by FAST
	double *Smin_fermi; //Minimum flux detectable by Fermi/LAT
	double *Smin_SKAmid; //Minimux flux detectable by SKAmid
	double *Smin_SKAlow; //Minimux flux detectable by SKAlow
	double *beta; //Spin orbit angle
	double *beta0; //Initial spin orbit angle
	double **fomega; //f_omega values taken from table 
	double Fmin;
	double Kr;
	double *Nb_orb; //Estimation of the number of orbits done by the NS
	double *delta; //Stores the gamma-ray peak separation
	double *Edot;
	double *Pdot;
	double *w_r_htru; //width of the radio profile observed with htru_pks
	double *xi; //angle in radians
	double *rho; //witdh of the beam in radian
	double *w_int; //Width of the radio profile (just geometry)
	double *w_r_fast; //width of the radio profile observed with fast 
	double *w_r_pmps; //width of the radio profile observed with pmps
	double *Fr; //radio flux table
	double *Fg; //gamma
	double *Lgamma; //L gamma
	double *subms;  //Did the pulsar went below the ms ? 1 or 0 for true or false
        double *alpha0; //a0= initial inclination angle
	double *dist; //disance from us
        _Bool vacuum,vacuum_evol,ff_evol;	
	const gsl_rng_type * T;
    	gsl_rng * r;  /* global generator */
};

#define SI_C 2.997924858e8 /* speed of light in units of m/s */
#define TMILKY 13.5e9 /* Age of the Milky Way in years */
#define KPC2CM 3.0856775807e21 /* kpc in cm */
#define SI_I 1e38   /*Moment of Inertia in kg.m2 */ 
#define R_NS 12000   /*Moment of Inertia in kg.m2 */ 
#define M_PI 3.14159265358979323846 /* pi */
#define SI_mu0 1.25663706212e-6 /* vacuum permeability in H/m */
#define G_grav 6.67430e-20 //gravitational constant km^-3 kg^-1 s^-2
#define MSUN 1.98847e30 //Solar mass kg
#define SI_eps0 8.85418782e-12 //Vacuum permittivity in F/m
#endif

