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
#include"galac_pot.h"

void distrib_vinit(void *params){ //Give an initial speed to a pulsar

    struct func_params *part= (struct func_params*)params;
    long np;
    double vx,vy,vz;
    double v=-1;
    double age_pulsar_yr;
    double cos_theta,phi;
    const double yr_sec=365*24*3600;

    for(np=0;np<part->Npulsars;np++){
       
       age_pulsar_yr    =  part->age_pulsar[np]/yr_sec;
       if (age_pulsar_yr<3*1e6) part->sigma_v=part->v_young;
       else part->sigma_v=part->v_old;

       //Computation of the velocity if we ignore the fact that pulsars are going in the same direction as the rotation axis when they are born
       /*vx = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v); //km/s
       vy = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
       vz = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);*/

       //Computation of the velocity if we do not ignore the fact above 
       cos_theta   =  2*gsl_rng_uniform(part->r)-1;
       phi         =  2*M_PI*gsl_rng_uniform(part->r);
       part->n_omega_x[np]=sqrt(1-sq(cos_theta))*cos(phi);
       part->n_omega_y[np]=sqrt(1-sq(cos_theta))*sin(phi);
       part->n_omega_z[np]=cos_theta;
       while (v<0){
	       v=sqrt(8.0/M_PI)*part->sigma_v+gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
	       //v=gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
       }
       //v=sqrt(8.0/M_PI)*part->sigma_v+gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
       vx=v*part->n_omega_x[np];
       vy=v*part->n_omega_y[np];
       vz=v*part->n_omega_z[np];


       part->vx0[np]=vx;part->vx[np]=vx;
       part->vy0[np]=vy;part->vy[np]=vy;
       part->vz0[np]=vz;part->vz[np]=vz;
       v=-1;

    }

}
/*
double phi_tot(void *params,long np){ //Compute the gravitational potential felt by a pulsar

      struct func_params *part= (struct func_params*)params;
      const double kpc2km=3.0856775807e16;
      double L0=1.0;
      double x=part->x[np];double y=part->y[np];double z=part->z[np];
      double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));double r=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km));
      double a1=0,b1=0.267*kpc2km,M1=1.02e10*MSUN;
      double a2=4.4*kpc2km,b2=0.308*kpc2km,M2=6.50535e10*MSUN;
      //double rc=6*kpc2km,Mc=5e10*MSUN;
      double Mh=12474*2.325*1e7*MSUN;double ah=7.7*kpc2km;
      double phi;
      double phi_1;
      double phi_2;
      double phi_k;
      double phi_NFW;
      phi_1=-(G_grav*M1)/(sqrt(sq(R)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1)))));
      phi_2=-(G_grav*M2)/(sqrt(sq(R)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2)))));
      phi_NFW=-(G_grav*Mh/(r))*log(1.0+(r/ah));
      phi_k=-(G_grav*4e6*MSUN)/r;
      phi=phi_1+phi_2+phi_NFW+phi_k;
      return phi;

}

double tot_energy(void *params,long np){ //Compute the sum of Kinetic and potential energy of a pulsar

      struct func_params *part= (struct func_params*)params;
      double vx,vy,vz;
      double v0=100.0;
      vx=part->vx[np]*v0;vy=part->vy[np]*v0;vz=part->vz[np]*v0;
      double phi=phi_tot(part,np);
      double E_tot=0.5*(sq(vx)+sq(vy)+sq(vz))+phi;
      return E_tot;

}

//Functions that compute the gradient of the different potentials

void phi_1(void *params,double phi1[3],long np){ //Potential for the bulge of the Milky Way

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double a1=0.0,b1=0.267*kpc2km,M1=1.02e10*MSUN;
       double L0=1.0;double v0=100.0;
       double T0=(L0*kpc2km)/v0;
       double x=part->x[np];double y=part->y[np];double z=part->z[np];
       double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));
       double phi1_x=-((G_grav*M1*x)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)))*(sq(T0)/L0);
       double phi1_y=-((G_grav*M1*y)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)))*(sq(T0)/L0);
       double phi1_z=-((G_grav*M1*z)/((pow(sq(R)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)*sqrt(sq(z*L0*kpc2km)+sq(b1)))))*(sq(T0)/L0)*(a1+sqrt(sq(z*L0*kpc2km)+sq(b1)));
       phi1[0]=phi1_x;phi1[1]=phi1_y;phi1[2]=phi1_z;
}


void phi_2(void *params,double phi2[3],long np){ //Potential for the disk of the Milky Way

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double a2=4.4*kpc2km,b2=0.308*kpc2km,M2=6.50535e10*MSUN;
       double L0=1.0;double v0=100.0;
       double T0=(L0*kpc2km)/v0;
       double x=part->x[np];double y=part->y[np];double z=part->z[np];
       double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));
       double phi2_x=-((G_grav*M2*x)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)))*(sq(T0)/L0);
       double phi2_y=-((G_grav*M2*y)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)))*(sq(T0)/L0);
       double phi2_z=-((G_grav*M2*z)/((pow(sq(R)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)*sqrt(sq(z*L0*kpc2km)+sq(b2)))))*(sq(T0)/L0)*(a2+sqrt(sq(z*L0*kpc2km)+sq(b2)));
       phi2[0]=phi2_x;phi2[1]=phi2_y;phi2[2]=phi2_z;
}

void phih(void *params,double phih[3],long np){ //Potential for Dark Matter halo but did not work

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double rc=6*kpc2km,Mc=5e10*MSUN;double L0=1.0;double v0=100.0;
       double T0=(L0*kpc2km)/v0;
       double x=part->x[np];double y=part->y[np];double z=part->z[np];
       double r=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km));
       double phih_x=-((G_grav*Mc/rc)*((x/(sq(rc)+sq(r)))-((rc/(sq(r)*r))*atan(r/rc)*x)+(x/sq(r))*(1/(1+sq(r/rc)))))*(sq(T0)/L0);
       double phih_y=-((G_grav*Mc/rc)*((y/(sq(rc)+sq(r)))-((rc/(sq(r)*r))*atan(r/rc)*y)+(y/sq(r))*(1/(1+sq(r/rc)))))*(sq(T0)/L0);
       double phih_z=-((G_grav*Mc/rc)*((z/(sq(rc)+sq(r)))-((rc/(sq(r)*r))*atan(r/rc)*z)+(z/sq(r))*(1/(1+sq(r/rc)))))*(sq(T0)/L0);
       phih[0]=phih_x;phih[1]=phih_y;phih[2]=phih_z;

}

void phi_NFW(void *params,double phi_NFW[3],long np){ //Potential for the dark matter halo of the galaxy 

	struct func_params *part= (struct func_params*)params;
        const double kpc2km=3.0856775807e16;
        double Mh=12474*2.325*1e7*MSUN;double ah=7.7*kpc2km;double L0=1.0;double v0=100.0;
        double T0=(L0*kpc2km)/v0;
        double x=part->x[np];double y=part->y[np];double z=part->z[np];
        double r=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km));
        double phi_x=-(G_grav*Mh*x/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        double phi_y=-(G_grav*Mh*y/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        double phi_z=-(G_grav*Mh*z/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        phi_NFW[0]=phi_x;phi_NFW[1]=phi_y;phi_NFW[2]=phi_z;

}

void phi_kepler(void *params,double phi_k[3],long np){ //Potential for the supermassive black hole at the center of the Milky Way

	struct func_params *part= (struct func_params*)params;
        const double kpc2km=3.0856775807e16;
        double Mc=4e6*MSUN;double L0=1.0;double v0=100.0;
        double T0=(L0*kpc2km)/v0;
        double x=part->x[np];double y=part->y[np];double z=part->z[np];
        double r_2=sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km);
        double phi_x=-(G_grav*Mc*x*sq(T0))/(pow(r_2,1.5)*L0);
        double phi_y=-(G_grav*Mc*y*sq(T0))/(pow(r_2,1.5)*L0);
        double phi_z=-(G_grav*Mc*z*sq(T0))/(pow(r_2,1.5)*L0);
        phi_k[0]=phi_x;phi_k[1]=phi_y;phi_k[2]=phi_z;

}

void evol_galac_pot_verlet(void *params){ //Verlet integration scheme

    struct func_params *part= (struct func_params*)params;
    const double kpc2km=3.0856775807e16;
    const double yr_sec=365*24*3600;long np;
    double L0=1.0;double v0=100.0;
    double T0=(L0*kpc2km)/v0;
    double step=1e4*yr_sec;
    double step_wd=step/T0;
    double step_fshift=step_wd*0.5;
    FILE *file=NULL;
    file=fopen("x_y_err_verlet.txt","w+");
    for(np=0;np<part->Npulsars;np++){

       part->x[np]=part->x[np]/L0;part->y[np]=part->y[np]/L0;part->z[np]=part->z[np]/L0;
       part->vx[np]=part->vx[np]/v0;part->vy[np]=part->vy[np]/v0;part->vz[np]=part->vz[np]/v0;

       double E0=tot_energy(part,np);
       for(double t=TMILKY*yr_sec-part->age_pulsar[np];t<TMILKY*yr_sec;t+=step){

         //double gphih[3];
	 double gphi_1[3],gphi_2[3];
	 double gphi_NFW[3];
         double grad_phi[3];
	 int i;

	 //Computation of the gradient with the space coordinates
         phi_1(part,gphi_1,np);phi_2(part,gphi_2,np);phi_NFW(part,gphi_NFW,np);//phih(part,gphih,np);//phi_kepler(pos,gphi_k);
         for(i=0;i<3;i++){

              grad_phi[i]=gphi_1[i]+gphi_NFW[i]+gphi_2[i];//gphi_k[i];

         }

	 //First shift of the velocities
	 part->vx[np]+=grad_phi[0]*step_fshift;
	 part->vy[np]+=grad_phi[1]*step_fshift;
	 part->vz[np]+=grad_phi[2]*step_fshift;

	 //First shift of the coordinates
         part->x[np]+=part->vx[np]*step_wd;
         part->y[np]+=part->vy[np]*step_wd;
         part->z[np]+=part->vz[np]*step_wd;

	 //Computation of the gradient with the space coordinates
         phi_1(part,gphi_1,np);phi_2(part,gphi_2,np);phi_NFW(part,gphi_NFW,np);//phih(part,gphih,np);//phi_kepler(pos,gphi_k);
         for(i=0;i<3;i++){

              grad_phi[i]=gphi_1[i]+gphi_NFW[i]+gphi_2[i];//gphi_k[i];

         }

	 //Second shift of the velocities 
	 part->vx[np]+=grad_phi[0]*step_fshift;
         part->vy[np]+=grad_phi[1]*step_fshift;
         part->vz[np]+=grad_phi[2]*step_fshift;


       }

     double Ef=tot_energy(part,np);
     part->x[np]=part->x[np]*L0;part->y[np]=part->y[np]*L0;part->z[np]=part->z[np]*L0;
     part->vx[np]=part->vx[np]*v0;part->vy[np]=part->vy[np]*v0;part->vz[np]=part->vz[np]*v0;

     part->x_s[np]= part->x[np]; // shift center of the Galaxy to the Sun
     part->y_s[np]= part->y[np]-8.5;
     part->z_s[np]= part->z[np]-0.015;
     part->dist[np]= sqrt(sq(part->x_s[np])+sq(part->y_s[np])+sq(part->z_s[np]));

     part->err_rel_g[np]=fabs((E0-Ef)/E0)*100;
     fprintf(file,"%e|%e|%e|\n",part->x_s[np],part->y_s[np],part->err_rel_g[np]);

    }
    
  fclose(file);

}

void evol_galac_PEFRL(void *params){ //PEFRL integration scheme

    struct func_params *part= (struct func_params*)params;
    double RAD = 180/M_PI;double gbr;double glr;
    const double yr_sec=365*24*3600;long np;
    double step;
    FILE *file=NULL;
    const double kpc2km=3.0856775807e16;
    double L0=1.0;double v0=100.0;
    double T0=(L0*kpc2km)/v0;
    //double step_wd=step/T0;
    double step_wd;
    double xsi=0.1786178958448091;double lambda=-0.2123418310626054;
    double chi=-0.6626458266981849e-01;
    file=fopen("x_y_err_PEFRL.txt","w+");
    for(np=0;np<part->Npulsars;np++){

       part->x[np]=part->x[np]/L0;part->y[np]=part->y[np]/L0;part->z[np]=part->z[np]/L0;
       part->vx[np]=part->vx[np]/v0;part->vy[np]=part->vy[np]/v0;part->vz[np]=part->vz[np]/v0;

       double E0=tot_energy(part,np);
       for(double t=TMILKY*yr_sec-part->age_pulsar[np];t<TMILKY*yr_sec;t+=step){

	 if(part->age_pulsar[np]<=1e2*yr_sec) {step=10*yr_sec;}
	 else if(part->age_pulsar[np]>1e2*yr_sec && part->age_pulsar[np]<=1e3*yr_sec) {step=25*yr_sec;}
	 else if(part->age_pulsar[np]>1e3*yr_sec && part->age_pulsar[np]<=1e4*yr_sec) {step=5e1*yr_sec;}
	 else if(part->age_pulsar[np]>1e4*yr_sec && part->age_pulsar[np]<=1e5*yr_sec) {step=5e2*yr_sec;}
	 else if(part->age_pulsar[np]>1e5*yr_sec && part->age_pulsar[np]<=1e6*yr_sec) {step=2.5e3*yr_sec;}
	 else if(part->age_pulsar[np]>1e6*yr_sec && part->age_pulsar[np]<=1e7*yr_sec) {step=2.5e4*yr_sec;}
	 else if(part->age_pulsar[np]>1e7*yr_sec && part->age_pulsar[np]<=1e8*yr_sec) {step=1e5*yr_sec;}
	 else if(part->age_pulsar[np]>1e8*yr_sec) {step=1e5*yr_sec;}

	 step_wd=step/T0;

         double gphi_1[3];
         double gphi_2[3];
         double gphi_NFW[3];
         double grad_phi[3];
	 double gphi_k[3];
         int i;

         //First shift of space coordinates
         part->x[np]+=part->vx[np]*step_wd*xsi;
         part->y[np]+=part->vy[np]*step_wd*xsi;
         part->z[np]+=part->vz[np]*step_wd*xsi;

         //Computation of the gradient with the new space coordinates
         phi_1(part,gphi_1,np);
         phi_2(part,gphi_2,np);
         phi_NFW(part,gphi_NFW,np);
	 phi_kepler(part,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //Shift of the velocities
         part->vx[np]+=grad_phi[0]*step_wd*0.5*(1-2*lambda);
         part->vy[np]+=grad_phi[1]*step_wd*0.5*(1-2*lambda);
         part->vz[np]+=grad_phi[2]*step_wd*0.5*(1-2*lambda);

         //second shift of space coordinates
         part->x[np]+=part->vx[np]*step_wd*chi;
         part->y[np]+=part->vy[np]*step_wd*chi;
         part->z[np]+=part->vz[np]*step_wd*chi;

         //Computation of the gradient with the new space coordinates
         phi_1(part,gphi_1,np);
         phi_2(part,gphi_2,np);
         phi_NFW(part,gphi_NFW,np);
	 phi_kepler(part,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //Second shift of velocities
         part->vx[np]+=grad_phi[0]*step_wd*lambda;
         part->vy[np]+=grad_phi[1]*step_wd*lambda;
         part->vz[np]+=grad_phi[2]*step_wd*lambda;

         //third shift of space coordinates
         part->x[np]+=part->vx[np]*step_wd*(1-2*(chi+xsi));
         part->y[np]+=part->vy[np]*step_wd*(1-2*(chi+xsi));
         part->z[np]+=part->vz[np]*step_wd*(1-2*(chi+xsi));

         //Computation of the gradient with the new space coordinates
         phi_1(part,gphi_1,np);
         phi_2(part,gphi_2,np);
         phi_NFW(part,gphi_NFW,np);
	 phi_kepler(part,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //third shift of velocities
         part->vx[np]+=grad_phi[0]*step_wd*lambda;
         part->vy[np]+=grad_phi[1]*step_wd*lambda;
         part->vz[np]+=grad_phi[2]*step_wd*lambda;

         //fourth shift of space coordinates
         part->x[np]+=part->vx[np]*step_wd*chi;
         part->y[np]+=part->vy[np]*step_wd*chi;
         part->z[np]+=part->vz[np]*step_wd*chi;

         //Computation of the gradient with the new space coordinates
         phi_1(part,gphi_1,np);
         phi_2(part,gphi_2,np);
         phi_NFW(part,gphi_NFW,np);
	 phi_kepler(part,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //fourth shift of velocities
         part->vx[np]+=grad_phi[0]*step_wd*0.5*(1-2*lambda);
         part->vy[np]+=grad_phi[1]*step_wd*0.5*(1-2*lambda);
         part->vz[np]+=grad_phi[2]*step_wd*0.5*(1-2*lambda);

         //final shift space coord
         part->x[np]+=part->vx[np]*step_wd*xsi;
         part->y[np]+=part->vy[np]*step_wd*xsi;
         part->z[np]+=part->vz[np]*step_wd*xsi;
       }

      //Save values in a file + computation energy and error
      double Ef=tot_energy(part,np);
      part->err_rel_g[np]=fabs((E0-Ef)/E0)*100;
      part->x[np]=part->x[np]*L0;part->y[np]=part->y[np]*L0;part->z[np]=part->z[np]*L0;
      part->vx[np]=part->vx[np]*v0;part->vy[np]=part->vy[np]*v0;part->vz[np]=part->vz[np]*v0;

      part->x_s[np]= part->x[np]; // shift center of the Galaxy to the Sun
      part->y_s[np]= part->y[np]-8.5;
      part->z_s[np]= part->z[np]-0.015;
      part->dist[np]= sqrt(sq(part->x_s[np])+sq(part->y_s[np])+sq(part->z_s[np]));
      double r= sqrt(sq(part->x_s[np])+sq(part->y_s[np]));
      if(part->dist[np]<1e-15){
            gbr=0;
      }else{
            gbr=asin(part->z_s[np]/part->dist[np]);
           }
      part->gb[np]=gbr*RAD;

      if(r<1e-15){
            glr=0;
                 }
      else{
            if(part->x_s[np]>=0){
                glr=acos(-part->y_s[np]/r);
                                }
	    else{
                glr=acos(part->y_s[np]/r)+M_PI;
                }
            }

      part->gl[np]=glr*RAD;
      fprintf(file,"%e|%e|%e|%e|%e|%e|%e|\n",part->x[np],part->y[np],part->z[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np]);

      
    }
    fclose(file);
}

void send_data(void *params){

     struct func_params *part= (struct func_params*)params;
     long np;
     FILE *file=NULL;
     file=fopen("data_for_py.txt","w+");
     const double kpc2km=3.0856775807e16;
     for(np=0;np<part->Npulsars;np++){

	 double R=sqrt(sq(part->x[np])+sq(part->y[np]));
	 double phi=atan(part->y[np]/part->x[np]);
	 double vr=(part->x[np]*part->vx[np]+part->y[np]*part->vy[np])/(sqrt(sq(part->x[np])+sq(part->y[np])));
	 double vphi=-(1.0/(1+sq(part->y[np]/part->x[np])))*((part->y[np]*(part->vx[np]/kpc2km)/(2.0*sq(part->x[np])))+(part->x[np]*(part->vy[np]/kpc2km)/(2.0*sq(part->y[np]))));
	 fprintf(file,"%e|%e|%e|%e|%e|%e|\n",R,vr,vphi,part->z[np],part->vz[np],phi);
     }
     fclose(file);


}*/

void nb_orbit(void *params){

        struct func_params *part= (struct func_params*)params;
        long np;
        const double kpc2km=3.0856775807e16;
        const double yr_sec=365*24*3600;
        double Mh=12474*2.325*1e7*MSUN;
        double ah=7.7*kpc2km;
        double r_km;
        double r2_dphi_dr;
        double P_orb;
        for(np=0;np<part->Npulsars;np++){
                r_km=sqrt(sq(part->x0[np])+sq(part->y0[np])+sq(part->z0[np]))*kpc2km;
                r2_dphi_dr=G_grav*Mh*(((1.0/r_km)*log(1.0+r_km/ah))-1.0/(r_km+ah));
                P_orb=(M_PI*2*sqrt(sq(r_km)/r2_dphi_dr))/(yr_sec);
                part->Nb_orb[np]=(part->age_pulsar[np]/yr_sec)/P_orb;
        }
}
