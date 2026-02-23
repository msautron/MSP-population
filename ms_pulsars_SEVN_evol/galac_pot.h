void distrib_vinit(void *params);
__device__ double phi_tot(void *params,long np);
__device__ double tot_energy(void *params,long np);
__device__ void phi_1(void *params,double phi1[3],long np);
__device__ void phi_2(void *params,double phi2[3],long np);
__device__ void phi_NFW(void *params,double phi_NFW[3],long np);
__device__ void phi_kepler(void *params,double phi_k[3],long np);
void send_data(void *params);
void nb_orbit(void *params);
