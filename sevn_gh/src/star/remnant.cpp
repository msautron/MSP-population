//
// Created by iorio on 7/5/21.
//

#include <remnant.h>
#include <star.h>


Staremnant::Staremnant(Star *s, double Mremnant) : born_time(s->getp(Worldtime::ID)),Mremnant_at_born(Mremnant) {}

double Staremnant::get(Star *s, size_t ID) {

    if (ID==Mass::ID)
        return Mass(s);
    else if(ID==Radius::ID)
        return Radius(s);
    else if(ID== Luminosity::ID)
        return Luminosity(s);
    else if(ID== OmegaRem::ID)
        return OmegaRem(s);
    else if(ID== Inertia::ID)
        return Inertia(s);
    else if(ID== Bmag::ID)
        return Bmag(s);
    else if(ID == Xspin::ID)
        return Xspin(s);
    else
        svlog.critical("Property with ID="+ utilities::n2s(ID,__FILE__,__LINE__)+" not implemented in Staremnant",__FILE__,__LINE__,sevnstd::notimplemented_error());

    return 0.;
}

double Staremnant::InertiaSphere(Star *s) const {
    return 0.4*s->getp(Mass::ID)*s->getp(Radius::ID)*s->getp(Radius::ID);
}

double Staremnant::age(Star *s) const {
    return s->getp(Worldtime::ID) - born_time;
}

double Staremnant::Mass(Star *s) const {
    return s->getp(Mass::ID);
}

/****BH****/
double BHrem::Radius(Star *s) const {
    return utilities::R_Schwarzschild(s->getp(Mass::ID));
}

double BHrem::Xspin(Star *s) const {

    if (s -> get_svpar_str("xspinmode")=="accretion"){
        return XspinAccretion(s);
    }

    return xspin;
}

double BHrem::estimate_Xspin(Star *s) const {

    const auto &spin_mode = s -> get_svpar_str("xspinmode");


    if (spin_mode == "disabled")
            return std::nan("");
    
    const double &z0 = s -> get_Z();
    const double &mco = s -> getp_0(MCO::ID); //Last value before the SN explosion, the SN explosion has already reset MCO at this point

    if (spin_mode == "geneva")
        return XspinGeneva(z0, mco);
    else if (spin_mode == "mesa")
        return XspinMESA(z0, mco);
    else if (spin_mode == "fuller")
        return XspinFuller();
    else if (spin_mode == "maxwellian"){
        const double &sigma_xspin = s -> get_svpar_num("xspin_sigma_maxwell");
        return XspinMaxwellian(sigma_xspin);
    }
    else if (spin_mode == "zeros"){
        return XspinZeros();
    }
    else if (spin_mode == "accretion"){
        return 0.0;
    }
    else{
        svlog.critical("xspinmode "+spin_mode+" not allowed",__FILE__,__LINE__,sevnstd::params_error(""));
    }

    return std::nan("");
}

void BHrem::default_initialiser(Star *s) {
    remnant_type=Lookup::Remnants::BH;
    set_Xspin(estimate_Xspin(s));
}

int BHrem::apply_Bavera_correction_to_Xspin(double period, double mass_wr) {
    double alpha, beta, period_days,_xspin;
    period_days = period * 365.25;


    if (period_days > 1.0)
        _xspin=xspin;
    else {
        /* These values are used if period and mass_wr are taken at C depletion. (Here we use the instant previous to the explotion.) */
        alpha = 0.029928 + exp(-0.282998 * mass_wr);
        alpha = -0.051237 / alpha;
        beta = 0.010905 + exp(-0.422213 * mass_wr);
        beta = -0.027090 / beta;
        _xspin = log10(period_days) * (alpha * log10(period_days) + beta);
        _xspin = std::min(std::max(0.0,_xspin),1.0); //Force xspin between 0 and 1
    }

    //Set the new value of _xspin
    set_Xspin(_xspin);

    return EXIT_SUCCESS;
}

double BHrem::XspinAccretion(Star *s) const {

    double Mratio = get_Mremnant_at_born() / s->getp(Mass::ID);
    double _xspin=1;
    if (Mratio>0.999999){
        _xspin=xspin;
    }
    else if (1/Mratio<=std::sqrt(6)){
            _xspin = std::sqrt(2./3.) * Mratio * (4. - std::sqrt(18*Mratio*Mratio-2));
    }

    return _xspin;
}


/****NS****/
double NSrem::Luminosity(Star *s) const {
    //cooling curve  LNS=0.02*M^(2/3)/(max(t,0.1)^2) from Eq. 93 using Hurley, 2000
    double Dt = std::max(age(s),0.1); //Time from the remnant creation
    return 0.02*std::pow(s->getp(Mass::ID),2./3.)/(Dt*Dt);
}

double NSrem::Inertia(Star *s) const {
    //Simple the ienrtia of a sphere 2/5 M R^2
    return InertiaSphere(s);
}

/*double NSrem::OmegaRem(_UNUSED Star *s) const {
    //CE: Fill the function, the return of this function is the Omega (in s) after DT from the NS birth
    //Useful constants
    const double yr = 3.1557600e7; //yr in seconds
    const double rsun = 6.95700e10; //rsun in cm
    const double msun = 1.98892e33; //msun in g

    //The first time it is called it is just Omega0
    if(age(s)==0){
        return Omega0;
    }
    //return s->getp(OmegaRem::ID);

    //TODO Radius and Inertia should be the one at the beginning of the evolution not their last values
    //In any case it does not matter a lot since Radius and Inertia will be constant during SSE
    const double R = s->getp(Radius::ID)*rsun; //Stellar radius at the beginning of the evolution
    const double I = s->getp(Inertia::ID)*rsun*rsun*msun; //Inertia  at the beginning of the evolution

    const double &Bi =  s->getp_0(Bmag::ID); //Initial B in Gauss in this timestep, Notice B should have been already evolved
    const double &Bf =  s->getp(Bmag::ID);//Final B in Gauss, Notice B should have been already evolved, we can also use Bmag(s)
    const double &Omegai = s->getp(OmegaRem::ID); //Initial value for Omega, Notice we use getp because we are evolving it now, so the  value of the previous step is still in getp
    const double & dt = s->getp(Timestep::ID)*1e6*yr; //Timestep
    const double taud = tau_magnetic*1e6*yr; //Magnetic field decay time scale in Myr;
    const double c_cm_s = utilities::c*rsun/yr; //Speed of light in cm/s

    double factor1 = 2*std::pow(R, 6)*sinalpha*sinalpha/(3*std::pow(c_cm_s, 3)*I);
    //I think it should be 8*M_PI*std::pow(R, 6)*sinalpha*sinalpha/(3*std::pow(c_cm_s, 3)*I);

    double factor2 = Bmin*Bmin*dt-taud*Bmin*(Bf-Bi)-taud/2*(Bf*Bf-Bi*Bi);

    double Omega = std::pow(pow(Omegai, -2)+2*factor1*factor2, -0.5);

    return Omega;
}*/

double func_angle_mhd(double x, void *params){ // function for which we want to find the roots Eq. 20 from Phillipov 2014

      param_MHD *part= (param_MHD*)params;
      double result;
      double k;
      double tau_d=part->tau_d;
      //double age = age(s);
      const double& dt = part->dt;
      double cosa0 = part->cosalpha;
      double t_used;
      double alpha_d=part->alpha_d;
      double tau_MHD=part->tau_MHD;
      double sin2a0=1.-pow(cosa0,2.0);

      t_used  =   (alpha_d*tau_d)*(pow(1.+dt/tau_d,1.-2/alpha_d)-1)/(alpha_d-2.); //Eq.19 from Dirson et al. 2022

      //t_used=dt;

      k  = t_used/tau_MHD + 0.5/sin2a0 + log(sqrt(sin2a0));

      return result= 0.5 /(x*x) + log(x) - k;

}

double evol_inc_angle(double cos2a0,param_MHD param_MHD_v,double *Pdot,double *Edot,double *alpha,double alpha0){

	double x_lo = 1.0e-17, x_hi = std::sqrt(1-cos2a0);
        double sina;
	double k1 = 7.388120797268824e-33; // en 1/s (8*M_PI*pow(R,6))/(3*mu0*cube(C)*I);
        double Omega;
        int status;
        double omega_dot;
        int iter=0;
        const gsl_root_fsolver_type *R2;
        gsl_root_fsolver *s2;
        if(fabs(func_angle_mhd(x_lo,&param_MHD_v))<1e-15){

                sina  = x_lo;

                } else if(fabs(func_angle_mhd(x_hi,&param_MHD_v))<1e-15){

                        sina  = x_hi;

                } else{

                        gsl_function  F;
                        F.function         =    &func_angle_mhd;
                        F.params           =    &param_MHD_v;
                        R2                 =    gsl_root_fsolver_brent;
                        s2                 =    gsl_root_fsolver_alloc(R2);

                        gsl_set_error_handler_off();
                        status=gsl_root_fsolver_set(s2,&F,x_lo,x_hi);

                        for(iter = 0; iter < 100; iter++) { // What happens if iter=100 is reached and no solutions found???

                                status = gsl_root_fsolver_iterate(s2);

                                double left_int  = gsl_root_fsolver_x_lower(s2);
                                double right_int = gsl_root_fsolver_x_upper(s2);

                                status = gsl_root_test_interval(left_int, right_int, 1.0e-10, 1.0e-15);
                                if(status != GSL_CONTINUE) {
                                        sina=(left_int+right_int)/2.;
                                        break;
                                                           }
                        }

                                gsl_root_fsolver_free (s2);

                        }

                if(status!=0) sina=std::sqrt(1-cos2a0);

                if (alpha0>M_PI/2) *alpha=M_PI-asin(sina);
		else if (alpha0<=M_PI/2) *alpha=asin(sina);

                Omega          = param_MHD_v.Omega0*(cos2a0*sina)/((std::sqrt(1-cos2a0))*pow(cos(asin(sina)),2.0));
                omega_dot      = 1.5*k1*pow(param_MHD_v.B*1e-4,2.0)*pow(Omega,3.0)*(1+pow(sina,2.0));
                *Pdot          = 2*M_PI*omega_dot/pow(Omega,2.0);
                *Edot          = param_MHD_v.I*1e-3*pow(1e-2,2.0)*Omega*omega_dot;

return Omega;

}

double omega_dot_calc(double omega,double alpha,const double dt,double B0,double M_NS,double taud,double Bmin,double ell){

        double MSUN=1.98847e30;
        double R_NS=12000; //Radius of a NS
        double G=6.67430e-11;
        double C=2.997924858e8;
        const double alpha_d = 1.5;
        const double mu_0=1.25663706212e-6;
        double B=B0*pow(1.0+(dt/taud),-1.0/alpha_d)+Bmin;
	//double B=B0;
        double mu=pow(4*M_PI*pow(B,2.0)*pow(R_NS,6.0)/mu_0,0.5);
        double r_lc=C/omega;
        double n3=-pow(mu,2.0)/(pow(r_lc,3.0));
        double I=(2.0/5.0)*M_NS*MSUN*pow(R_NS,2.0);
        double omega_dot;
	double n_gw=-(4.0/5.0)*((G*pow(I*ell,2.0)*pow(omega,5.0))/(pow(M_PI,2.0)*pow(C,5.0)));
        omega_dot=(1.0/I)*(n3*(1+pow(sin(alpha),2.0))+n_gw);
        return omega_dot;
}

double alpha_dot_calc(double omega,double alpha,const double dt,double B0,double M_NS,double taud,double Bmin){

        double R_NS=12000; //Radius of a NS
        double MSUN=1.98847e30;
        double C=2.997924858e8;
        const double mu_0=1.25663706212e-6;
        const double alpha_d = 1.5;
	//double B=B0;
        double B=B0*pow(1.0+(dt/taud),-1.0/alpha_d)+Bmin;
        double I=(2.0/5.0)*M_NS*MSUN*pow(R_NS,2.0);
        double mu=pow(4*M_PI*pow(B,2.0)*pow(R_NS,6.0)/mu_0,0.5);
        double r_lc=C/omega;
        double n3=-pow(mu,2.0)/(pow(r_lc,3.0));
        double alpha_dot;
        alpha_dot=(n3*sin(alpha)*cos(alpha))/(I*omega);
        return alpha_dot;
}

double evo_iso_pulsar(double omega,double alpha,const double dt,double B0,double M_NS,double taud,double Bmin,double ell,double *Pdot,double *Edot,double *alpha_end){

    //RK4 to make evolve alpha,beta and omega
    double k1,k2,k3,k1_alpha,k2_alpha,k3_alpha,k4,k4_alpha;
    double omega_dot;
    double MSUN=1.98847e30;
    double R_NS=12000; //NS radius in m
    double I=(2.0/5.0)*M_NS*MSUN*pow(R_NS,2.0);

    //Computation of k1 for each variable
    k1=omega_dot_calc(omega,alpha,0,B0,M_NS,taud,Bmin,ell);
    k1_alpha=alpha_dot_calc(omega,alpha,0,B0,M_NS,taud,Bmin);

    //Computation of k2 for each variable
    k2=omega_dot_calc(omega+k1*0.5*dt,alpha+k1_alpha*0.5*dt,dt/2.0,B0,M_NS,taud,Bmin,ell);
    k2_alpha=alpha_dot_calc(omega+k1*0.5*dt,alpha+k1_alpha*0.5*dt,dt/2.0,B0,M_NS,taud,Bmin);

    //Computation of k3 for each variable
    k3=omega_dot_calc(omega+k2*0.5*dt,alpha+k2_alpha*0.5*dt,dt/2.0,B0,M_NS,taud,Bmin,ell);
    k3_alpha=alpha_dot_calc(omega+k2*0.5*dt,alpha+k2_alpha*0.5*dt,dt/2.0,B0,M_NS,taud,Bmin);

    //Computation of K4 for each variable
    k4=omega_dot_calc(omega+k3*dt,alpha+k3_alpha*dt,dt,B0,M_NS,taud,Bmin,ell);
    k4_alpha=alpha_dot_calc(omega+k3*dt,alpha+k3_alpha*dt,dt,B0,M_NS,taud,Bmin);

    //RK4 scheme
    omega+=(dt/6.0)*(k1+2*k2+2*k3+k4);
    alpha+=(dt/6.0)*(k1_alpha+2*k2_alpha+2*k3_alpha+k4_alpha);

    //std::cout << "After isolated evolution omega is: " << omega << " Hz" << std::endl;
    //std::cout << "After isolated evolution alpha is: " << alpha*180/M_PI << " degrees " << std::endl;
    //std::cout << "Time step is: " << (dt/(365*24*3600))*1e-6 << " Myr" << std::endl;

    //Important quantities to compute
    omega_dot      = k1;
    *Pdot          = -2*M_PI*omega_dot/pow(omega,2.0);
    *Edot          = -I*omega*omega_dot;
    *alpha_end=alpha;
    return omega;

}

double NSrem::OmegaRem(_UNUSED Star *s) {

	//CE: Fill the function, the return of this function is the Omega (in s) after DT from the NS birth
    	//Useful constants
    	const double yr = 3.1557600e7; //yr in seconds
    	const double rsun = 6.95700e10; //rsun in cm
    	const double msun = 1.98892e33; //msun in g
	double Omega;
	//Useful parameters
	double cos2a0 = cosalpha*cosalpha;
	double alpha0=acos(cosalpha);
	double Pdot,Edot;
	double alpha_end,alpha,alpha_send;
	const double & dt = s->getp(Timestep::ID)*1e6*yr; //Timestep in s
	/*FILE *for_PPdotplot; //To comment every line related to it if simulating the population
	for_PPdotplot=fopen("/home/matteo.sautron/sevn_gh/SEVN_run/to_plotPPdot.txt","a+");*/
	/*int save1=0;
	int save2=0;
	int save3=0;
	int save4=0;
	int save5=0;*/
	//Params params = {s,this};
	param_MHD val_param_MHD;
	val_param_MHD.tau_d=tau_magnetic;
	val_param_MHD.dt=s->getp(Timestep::ID);
	val_param_MHD.cosalpha=cosalpha;
	val_param_MHD.alpha_d=alpha_d;
	val_param_MHD.tau_MHD=tau_MHD;
	val_param_MHD.Omega0=Omega0;
	val_param_MHD.B=Bmag(s);

    	//The first time it is called it is just Omega0
    	if(age(s)==0){
		M0=s->getp(Mass::ID);
        	return Omega0;
    		}
    	//return s->getp(OmegaRem::ID);

    	//TODO Radius and Inertia should be the one at the beginning of the evolution not their last values
    	//In any case it does not matter a lot since Radius and Inertia will be constant during SSE
    	const double R = s->getp(Radius::ID)*rsun; //Stellar radius at the beginning of the evolution
    	const double I = s->getp(Inertia::ID)*rsun*rsun*msun; //Inertia  at the beginning of the evolution
	val_param_MHD.I=I;
	
	//Computation of alpha thanks to usual FFE differential equations system 
	if (Omega0>=960) {Omega=evo_iso_pulsar(Omega0,alpha0,dt,Bmag(s)*1e-4,s->getp(Mass::ID),tau_magnetic*1e6*yr,Bmin*1e-4,ell,&Pdot,&Edot,&alpha_end);alpha_send=alpha_end;}
	//Computation of alpha thanks to equation 19 from Dirson et al. 2022
	else if (Omega0<960) {Omega=evol_inc_angle(cos2a0,val_param_MHD,&Pdot,&Edot,&alpha,alpha0); alpha_send=alpha;}
	if (acc_or_not==false) {set_cosalpha(cos(alpha_send));Omega0=Omega;}

	double P=2*M_PI/(Omega0);
	set_Pdot(Pdot);
	if (fabs(P)<1e-3) {sub_ms_or_not=1;}
	/*if (age(s)>=1e-1 && save[1000]==0) {P_saved[1000]=fabs(P);P_dot_saved[1000]=fabs(Pdot);save[1000]=1;}
	for(int j=1;j<nb_saved;j++){
		if (age(s)>=j*real_age/nb_saved && age(s)<(j+1)*real_age/nb_saved && save[j]==0){
                P_saved[j]=fabs(P);P_dot_saved[j]=fabs(Pdot);age_saved[j]=age(s)*1e6*365*24*3600;
                save[j]=1;
        	}
	}*/
	//std::cout << "P is : " << P << " in s, and alpha is : " << alpha*180/M_PI << " in degrees, and Pdot is : " << Pdot << " in s/s, Age = " << age(s) << " in Myr " << std::endl;
	acc_activate(false);
	if (age(s)==real_age and !isnan(P)) {
		FILE *data_ms_pulsars;
		FILE *data_taud;
		FILE *data_beta;
		FILE *data_beta_init;
		FILE *data_subms;
		/*FILE *for_PPdotplot; //To comment every line related to it if simulating the population
        	for_PPdotplot=fopen("/home/matteo.sautron/sevn_gh/SEVN_run/to_plotPPdot.txt","a+");*/
		double spin_up_val;
		char eonf[5];
		long end_name_file;
		char cmd[500];
		char cmd2[500];
		char cmd3[500];
		char cmd4[500];
		char cmd5[500];
		char cmd6[500];
		if(spin_up==true && acc_or_not==false) {spin_up_val=1.0;}
		else if (spin_up==true && acc_or_not==true) {spin_up_val=2.0;}
		else {spin_up_val=0.0;}
		end_name_file=s->get_svpar_num("end_of_name_file");
		sprintf(eonf,"%ld",end_name_file);
		sprintf(cmd,"/home/matteo.sautron/sevn_gh/SEVN_run/data_pop_ms%s.txt",eonf); 
		sprintf(cmd2,"python3 /home/matteo.sautron/sevn_gh/SEVN_run/get_info_binary.py %ld",end_name_file);
		sprintf(cmd3,"/home/matteo.sautron/sevn_gh/SEVN_run/data_taud%s.txt",eonf);
		sprintf(cmd4,"/home/matteo.sautron/sevn_gh/SEVN_run/data_beta%s.txt",eonf);
		sprintf(cmd5,"/home/matteo.sautron/sevn_gh/SEVN_run/data_beta_init%s.txt",eonf);
		sprintf(cmd6,"/home/matteo.sautron/sevn_gh/SEVN_run/data_subms%s.txt",eonf);
        	data_ms_pulsars=fopen(cmd,"a+");
		data_taud=fopen(cmd3,"a+");
		data_beta=fopen(cmd4,"a+");
		data_beta_init=fopen(cmd5,"a+");
		data_subms=fopen(cmd6,"a+");
		/*for(int i=1;i<nb_saved+1;i++){
			fprintf(for_PPdotplot,"%e %e %e\n",P_saved[i],P_dot_saved[i],age_saved[i]);
		}*/
		fprintf(data_ms_pulsars,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e ",fabs(P),fabs(Pdot),Edot,acos(cosalpha),Bmag(s)*1e-4,age(s)*1e6*365*24*3600,Pinit,B0*1e-4,alpha0,tps_acc,spin_up_val,s->getp(Mass::ID),M0,ell);
		fprintf(data_taud,"%e\n",tau_magnetic);
		fprintf(data_beta,"%e\n",beta);
		fprintf(data_beta_init,"%e\n",beta0*180.0/M_PI);
		fprintf(data_subms,"%d\n",sub_ms_or_not);
		fclose(data_ms_pulsars);
		fclose(data_taud);
		fclose(data_beta);
		fclose(data_beta_init);
		fclose(data_subms);
		//fclose(for_PPdotplot);
		system(cmd2);
	}
	return Omega;
}

double NSrem::Bmag(_UNUSED Star *s) const {

    //CE: Fill the function, the return of this function is the Bmag (in Gauss) after DT from the NS birth
    //The first time it is called it is just B0
    if(age(s)==0){
        return B0;
    }
    //return s->getp(Bmag::ID);


    //Estimate vaule for this timestep
    const double& dt = s->getp(Timestep::ID);
    const double& Bold = s->getp(Bmag::ID); //Old value of the magnetic field
    const double taud = tau_magnetic; //Magnetic field decay time scale in Myr;

    //double B=Bold;
    double B = Bold*pow(1.0+(dt/taud),-1.0/alpha_d)+Bmin; //Decay as a power law as in Sautron et al. (2024) 
    //double B = (Bold-Bmin)*std::exp(-dt/taud)+Bmin;
    //std::cout << "Magnetic field is  : " << B*1e-4 << " T (normal decay)" << std::endl;

    return B;
}


void NSrem::default_initialiser(Star *s) {
    remnant_type=Lookup::Remnants::NS_CCSN;
    FILE *info_birth;
    /*FILE *for_PPdotplot; //To comment every line related to it if simulating the population
    for_PPdotplot=fopen("/home/matteo.sautron/sevn_gh/SEVN_run/to_plotPPdot.txt","a+");*/
    double C=2.997924858e8;
    double mu_0=1.25663706212e-6;
    char cmd[500];
    char eonf[5];
    long end_name_file;
    end_name_file=s->get_svpar_num("end_of_name_file");
    sprintf(eonf,"%ld",end_name_file);
    sprintf(cmd,"/home/matteo.sautron/sevn_gh/SEVN_run/listBin%s.dat",eonf);
    info_birth=fopen(cmd,"r");
    double trash;char trash_str[100];
    fscanf(info_birth,"%le%s %le %le %s %le %le %le %le %s %s %le %le %le %le",&trash,trash_str,&trash,&trash,trash_str,&trash,&trash,&trash,&trash,trash_str,trash_str,&trash,&trash,&real_age,&trash);
    //Drawn a random alpha (we store the sinalpha)
    //Assuming an isotropic magnetic axis direction we have to uniformly sample in cos(alpha) between -1 and 1.
    cosalpha = generate_uniform(-1,1);
    alpha0=acos(cosalpha);
    double cosbeta = generate_uniform(-1,1);
    std::cout << "cosalpha : " << cosalpha << std::endl;  
    //sinalpha = std::sqrt(1-cosalpha*cosalpha);
    double alpha = acos(cosalpha);
    beta = acos(cosbeta);
    beta0= beta;
    //fprintf(data_beta_init,"%e\n",beta*180.0/M_PI);
    //fclose(data_beta_init);


    //TODO Different options for the B0 and Omega0 initilisation are now hardcoded.
    //we should add a new set of parameters to do it at runtime or a more sophisticated method to add alternative NS models

    //Initiliase B0
    //Generate uniform (see Sgalletta+23)
    //double B0min = 1e10; //CE: this is the minimum value for the uniform distribution of B0 in Gauss
    //double B0max = 1e13; //CE: this is the maximum value for the uniform distribution of B0 Gauss
    Bmin = generate_uniform(1e3,1e4)*1e4; //Bmin in Gauss
    //B0 = generate_uniform(B0min,B0max);

    //Generate B0 from log-normal distribution (Igoshev et al. (2022))
    double Bmean=275422870.33381635;
    //double Bmean=3e8;
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::normal_distribution<double> distributionB(0,0.5);
    B0=pow(10,log10(Bmean)+distributionB(gen))*1e4;
    //std::cout << "Initial magnetic field of the NS is : " << B0 << " G" << std::endl;


    //Generate flat-in-log (see Sgalletta+23)
    //double logB0 = generate_uniform(10.0, 13.0);
    //B0 = std::pow(10.0, logB0);
    //Generate Faucher-Giguere: log normal distribution
    //double logB0 =  (std::normal_distribution<>{12.65, 0.55})(utilities::mtrand);
    //B0 = std::pow(10.0, logB0);

    //Initialise Omega0
    //double Period0min = 0.01; //CE: this is the minimum value for the uniform distribution of Omega0 in s  (10 ms)
    //double Period0max = 0.1; //CE: this is the maximum value for the uniform distribution of Omega0 in s (100 ms)
    //double Period0 = generate_uniform(Period0min,Period0max);
    //Generate Faucher-Giguere: normal distrution with mean 0.3 s and variance 0.15 s (see Sgalletta+23)
    //double Period0 = (std::normal_distribution<>{0.3, 0.15})(utilities::mtrand);
    
    //Generate P0 from log-normal distribution (Igoshev et al. (2022))
    double Pmean=129e-3;
    std::normal_distribution<double> distributionP(0,0.45);
    double Period0=pow(10,log10(Pmean)+distributionP(gen));
    Pinit=Period0;
    double Pdot_init,omega_dot_init;
    omega_dot_init=-pow(2.0*M_PI/Pinit,3.0)*(4*M_PI*pow(12000,6.0)*pow(B0*1e-4,2.0)*(1.0+pow(sin(alpha),2.0)))*(1.0/(mu_0*pow(C,3.0)*1e38));
    Pdot_init=-omega_dot_init*(pow(Pinit,2.0))/(2*M_PI);
    /*fprintf(for_PPdotplot,"%e %e\n",Pinit,Pdot_init);
    for(int i=0;i<nb_saved;i++){
	    save[i]=0;
    }*/

    //Generate ellipticity of the NS
    //double ellmean=9.0;
    //std::normal_distribution<double> distributionell(ellmean,1.2);
    //ell=pow(10,-distributionell(gen));
    ell=pow(10,-generate_uniform(5,13));

    std::cout << " Ellipticity of the NS : " << ell << std::endl;
    std::cout << " Periode initiale : " << Period0 << " s " << std::endl;

    Omega0 = 2*M_PI/Period0;

    fclose(info_birth);
    //fclose(for_PPdotplot);
    //Initialise
    double k1    = 7.388120797268824e-33; // en 1/s (8*M_PI*pow(R,6))/(3*mu0*cube(C)*I);
    tau_magnetic = s->get_svpar_num("ns_magnetic_tscale"); //Magnetic field decay time scale in Myr
    //alpha_d      = s->get_svpar_num("ns_alpha_d"); //value of parameter alpha_d for magnetic field decay power law
    tau_0        = (1.0/k1)*pow(1.0/(B0*1e-4),2.0)*pow(0.5*Period0/M_PI,2.0); //Useful to compute alpha during spin down in s
    tau_MHD      = (1.0/(1e6*365.0*24.0*3600.0))*(2.0/3.0)*tau_0*((1.0-pow(cos(alpha),2.0))/(pow(cos(alpha),4.0))); //Useful to compute alpha during spin down in Myr
}

void NSrem::print_log_message(Star *s) {
    std::string w = utilities::log_print("NS",s,remnant_type,get_Mremnant_at_born(),B0,Omega0,cosalpha);
    s->print_to_log(w);
}


/****WD****/
double WDrem::Radius(Star *s) const {
    //Eq.91 from Hurley+00
    const double RNS = NSrem::Rns; //Neutron star radius 10 km in Rsun, from Hurley, 2000
    const double &Mch = s->get_svpar_num("sn_Mchandra"); //Mchandra
    const double &Mass = s->getp(Mass::ID);

    double RWD = 0.0115*std::sqrt(pow(Mch/Mass,0.6666666667) -  pow(Mass/Mch,0.6666666667));

    return std::max(RNS,RWD);
}

double WDrem::Luminosity(Star *s) const {
    //Eq. 90 in Hurley+00
    double num = 635*s->getp(Mass::ID)*pow(s->get_Z(),0.4);
    double den = A_luminosity* pow(age(s) + 0.1,1.4);

    return num/den;
}


void Zombierem::initialise(Star *s) {
    remnant_type=Lookup::Remnants::Zombie;
    luminosity = s->getp(Luminosity::ID);
    radius = s->getp(Radius::ID);
}

