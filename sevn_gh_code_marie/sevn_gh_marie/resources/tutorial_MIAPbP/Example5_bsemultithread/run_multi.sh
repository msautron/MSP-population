#!/usr/bin/env bash
cmdargument=$1

RED="\033[1;31m"
GREEN="\033[1;32m"
NC="\033[0m" # No Color






#########PARAMETERS#############################

#-------------------------------
#SEVNPATH
#-------------------------------
SEVN="/Users/giulianoiorio/Documents/SEVNp" #Complete path to the SEVN folder
SEVNEXE="${SEVN}/build/exe/sevnB.x" #Complete path to SEVN BSE executable

#-------------------------------
#RUN OPTIONS
#-------------------------------
NTHREADS="3" #Number of OpenMP threads (1 means no parallel threads, sequential execution)
NCHUNK="1000" #Evolve Nchunk at time
DTOUT="list"      #If list use the dtout reported in the input list, otherwise use this value for all the stars and binaries (Can be a number in Myr (e.g. 10), a colon separated sequence in Myr (e.g. 10:100:10 goes from 10 Myr to 100 Myr with 10 Myr step, or a list of comma separated numbers in Myr inside curly brackets (e.g. {10,50,40}))


#-------------------------------
#Input and Tables
#-------------------------------
if [ "$cmdargument" = "" ]; then
  LISTBIN="listBin.dat" #Complete path to input file (list of binaries or single stars)
else
  LISTBIN=$cmdargument
fi
IBMODE="new" #Input file format for binaries [new*] [legacy] [sevn1]
TABLES="${SEVN}/tables/SEVNtracks_parsec_ov05_AGB" #Complete path to look-up tables
TABLESHE="${SEVN}/tables/SEVNtracks_parsec_pureHe36" #Complete path to look-up tables for pure-He stars
TEND="list"
TSTART="list"
RSEED="false"




#-------------------------------
#OUTPUT
#-------------------------------
OUTPATH="sevn_output/" #Complete path to the output folder (the folder will be automatically created or cleaned if it already exists)
OMODE="csv" #Format for output files [h5] [ascii*] [csv]
NAMEPREX=""  #prefix to add to the name of the systems
LOGLEVEL="critical"  #Log output level: [debug] [info] [warning] [error] [critical]
LITPHASES="false" #Use literal phases instead of numbers in output [true] - [false*]
LOGFILE="true" #If true produce the logfile output  [true] - [false*]
SCOL="Mass:MHE:MCO:Radius:Luminosity:Temperature:Phase:RemnantType" #Additional columns to print in the output file for single stellar evolution runs. Default is empty, but any property of single stars can be added (check names in the Property class)
BCOL="BWorldtime:Semimajor:Eccentricity:BEvent" #Additional columns to print in the output file for binary stellar evolution runs. Default is empty, but any property of binary stars can be added (check names in the BinaryProperty class)


#-------------------------------
#METALLICITY
#-------------------------------
Z="list"  #Stellar metallicity - [*list][number]. If list use the Z in the input file otherwise overwrite all the Zs.

#-------------------------------
#PRESCRIPTIONS
#-------------------------------
WINDSMODE="hurley" #Prescriptions for wind accretion and the associated orbital changes - [hurley*] [disabled]
TIDES="tides_simple" #Prescriptions for tides - [simple*] [disabled]
GWMODE="peters" #Prescriptions for gravitational-wave decay - [peters*] [disabled]
RLMODE="hurley_rl" #Prescriptions for Roche-Lobe overflow and mass transfer/accretion - [hurley_rl*] [hurley_bse] [disabled]
CIRCMODE="periastron" #Prescriptions for orbit circularisation at onset of RLO - [periastron*] [periastron_full] [angmomg] [semimajor] [disabled]
CEMODE="energy" #Prescriptions for common-envelope evolution - [energy*] [disabled]
MIXING="simple" #Prescriptions  for mixing (merger) - [simple*] [disabled]
COLLMODE="hurley" #Prescriptions for collision at periastron - [hurley*] [disabled]
HARDMODE="disabled" #Prescriptions for hardening in stellar clusters [disabled*] [fastcluster]
SNORBCHANGE="hurley" #Prescriptions for orbital changes after SN kicks - [hurley*] [disabled]
SNKICKS="unified" #Prescriptions for SN kicks - [unified*] [hobbs] [hobbs_pure] [cc15] [ecus30] [ec15cc265] [zeros]
SNPISN="mapelli20" #Prescription for pair instability SN - [mapelli20*] [iorio22] [iorio22_limited] [farmer19] [disabled]
SNNML="lattimer89" #Prescription for neutrino mass loss in SN - [lattimer89*] [disabled]
SUPERNOVA="list" #Prescription for the SN explosion mechanism - [list*] [rapid] [rapid_gauNS] [delayed] [delayed_gauNS] [compact] [deathmatrix] [directcollapse]
BHXSPIN="disabled" #Prescription for the BH spin - [disabled*] [geneva] [mesa] [fuller] [maxwellian] [accretion]

#---------OPTIONS---------------
SQHE="false" #If true enable the Quasi Homogeneous Evolution  after a RLO mass transfer following Elrdige&Stanway11
TABCONV="true" #If true estimate the properties of the convective envelope using the tables (xxxconv.dat)
TABXSUP="false"  #If true use the surface abundance tables (xxxsup.dat)
TABINERTIA="false" #If true estimate the properties of the stellar inerita using the tables (inertia.dat)
INERTIAMODE="Hurley" #option for inertia estimate when tabuse_inertia is false [*Hurley][DeMink][hspherecore][hsphere]
TABRHE="true" #If true estimate the properties of the HE core radius using the tables (rhe.dat)
TABRCO="true" #If true estimate the properties of the CO core radius using the tables (rhe.dat)
THGHURLEY="false" #If true estimate the HG time from the Hurley+00 functional forms instead of using the convective envelope
OPTIMISTIC="false" #If true allow the star in the HG (Hurley phase 2) to start a CE
#-------------------------------
#Parameters
#-------------------------------
#-------GW-------#
GWTSHOLD="1" #Enable GW decay if GW_time_decay < GWTSHOLD*Hubble_time
GWONLYBCO="false" #If true activate the GW orbital decay on for binary compact objects
#-------RLO-------#
RLOEPSNOVA="0.001" #Fraction of accreted matter retained in nova eruption
RLOMACCR="0.5" #Fraction of the mass lost by the primary that is accreted onto the secondary during RLO
RLOGAM="-2" #Angular momentum lost during RLO. [-1]: from the primary, [-2]: from the secondary, [>0]: fraction lost from the system
RLOSTABILITY="qcrit_Hradiative_stable"  #"Option for RLO mass transfer stability
RLONTMAX="5"    #Max value of the mass to use in the normalisation of the nuclear mass transfer (Eq. 59 Hurley+02)
RLOCOLLISION="false"  #If true allow collision at periastron during RLO
RLOSMTMS="true"   #If true mass transfer from radiative MS and pureHE MS are always stable
#-------ACCRETION-------#
EDDF="1" #Eddington factor to limit accretion on a compact object (>1 means super-Eddington)
#-------SN-------#
MCHANDRA="1.44" #Chandrasekar mass limit for WD formation
SNLOWECSN="1.38" #Minimum value for the CO mass to go ECSN
SNLOW="1.44" #Minimum CO value for the CO mass to go SN (i.e. max CO mass for ECSN)
SNLOWECSNHE="-1"  #Minimum value for the CO mass to go ECSN for pureHe star, if -1 use the same value as H star
SNLOWHE="-1" #Minimum CO value for the CO mass to go SN (i.e. max CO mass for ECSN) for pureHe star, if -1 use the same value as H star
SNC25TS="0.35" #csi25 parameter threshold for explosion/implosion decision, if -1 use a stochastic threshold based on the results of  Patton&Sukhbold20
SNCOMPFB="0.9" #Fallback fraction for implosions in the compact SN option
SNMINVKICK="0.0" #Minimum SN Kick after all the corrections
SNVKICKSTD="265.0" #Standard deviation  of the Maxwellian distribution of kick velocity (Used in the Hobbs and Unified SN kick model)
#-------WINDS-------#
WALPHA="1.5" #alpha factor to tune the amount of wind accretion (Eq.6 Hurley+02)
WBETA="0.125" #beta factor to tune wind velocity (Eq.9 in Hurley+02)
#-------CE-------#
CEALPHA="3" #alpha in binding energy (Eq. 73 in Hurley02)
CELAM="-1" #if >0 Constant Lambda in binding energy (Eq. 69 in Hurley02). If -1 same Lambda as in BSE (-11 and -12 are other BSE-like variations).
#-4 Lambda interpolated from Klencki21 (-41 is a variation in which the Lambda is not interpolated but it is "quantised" in bins)
#-5 Lambda interpolated from  Xu&Li10  (-51 is a variation in which the Lambda is not interpolated but it is "quantised" in bins)
CELAMHE="0.5" #Constant Lambda in binding energy used for pureHe stars(Eq. 69 in Hurley02).
# Notice: some Lambda model have their own estimate of Lamdba_he (for example option -4), in this case this value is not considered
CELAMFTH="1"  #Fraction of internal energy that goes to the binding energy. Used only if star_lambda<0. Notice that some Lambda models do not include an option for the fraction of internal energu (e.g. option -5)
CEKCE="1"  #Fraction of non core mass  participating to the CE (e.g. envelope of giants) retained after the CE coalescence.If -1, use a rescaled version of eq. 77 In Hurley
CEKNCE="1" #Fraction of non core mass not participating to the CE (e.g. a MS star) retained after the CE coalescence.If -1, use the eq. 77 in Hurley 2002 (ce_kce is ignored)
#-------Hardening-------#
HARDRHOC="39900" #central density of the cluster in Msun/pc^3
HARDSIGMA="5" #3D velocity dispersion of the cluster in km/s
HARDXI="3" #xi parameter for fastcluster hardening option (Eq. 11 Mapelli+21)
HARDKAPPA="0.1" #kappa parameter for fastcluster hardening option (Eq. 13 Mapelli+21)
HARDMASS="1" #Average mass of the perturber stars in the environment producing the hardening in Msun
#-------NS-------#
NSMAX="3.0" #Maximum NS mass
NSMAGTSCALE="1000" #Magnetic field decay timescale in Myr
NSMAGMSCALE="0.15" #Magnetic field decay mass-scale in Msun
NSMASSMEAN="1.33" #NS masses are drawn from a Gaussian with this mean. Notice, not all the SNMODE options allows to use it
NSMASSSTD="0.09" #NS masses are drawn from a Gaussian with this std. Notice, not all the SNMODE options allows to use it
#-------BH-------#
MAXWSDXSPIN="0.1" # Standard deviation of the Maxwellian distribution for Xspin - default: 0.1.
BAVERAXSPIN="false" # Bavera correction for the black-hole spin - default: false.


# -------------------------------
#Advanced parameters - DO NOT CHANGE THEM UNLESS YOU ARE AWARE OF WHAT YOU ARE DOING
#-------------------------------

# -------------------------------
#Parameters for jumping onto new tracks
#-------------------------------
JTEMAX="0.005" #Maximum relative error in mass when jumping on a new track
JTMAXDM="1.2" #Maximum new ZAMS tested when jumping on a new track (Mzams_new_max = Mzams_old + JTMAXDM*DM_accreted_or_donated)
JTMINDM="0" #Minimum new ZAMS tested when jumping on a new track (Mzams_new_min = Mzams_old + JTMINDM*DM_accreted_or_donated)
JTMAXITER="10" #Maximum number of iterations for reaching convergence
JTDMTSHOLD="0.01" #Maximum relative change in total mass for not changing track
# -------------------------------
#Other parameters
#-------------------------------
MAXREP="50" #Maximum number of repetitions allowed in the sse and bse. If we reach this number an error is raised
NAKEDTS="1E-4" #Mass difference threshold (Msun) between envelope and core to set a star as nakedHe or nakedCO.
WRTS="0.021" #Relative difference threshold between envelope (Mass-MHE) and total mass to define a star as Wolf Rayet
INITERRSTOP="false" #If true terminate the run when a error on a system initialisation is thrown
SMAXCO="false" #If true the first time a star develops a CO core, we set the maximum CO core Mass for SSE as the last value of the interpolating tracks
SMINHE="false" #If true the first time a star develops a CO core, we set the minimum HE core Mass for SSE as the last value of the interpolating tracks
# -------------------------------
#Interpolation
#-------------------------------
INTW="linear" #Option for setting the weights in the Mass interpolation [*linear][rational][log], Notice for Radius, Inertia and Luminosity it set to log
INTWLOG="log" #Option to set weights for mass interpolation for log properties (i.e. Radius,Luminosity,Inertia) [linear][rational][*log]
INTWPHASE="rational" #Option to set weights for estimating the stellar phases times [linear][*rational][log]
# -------------------------------
#Timestep
#-------------------------------
TSMAXVAR="0.05"  #Relative maximum variation of stellar and binary properties used in the adaptive time step
TSMIN="-1" #Force the adaptive timestep to be larger than this value, it will it has the priority on any other option, -1 means that the option is disabled
TSMAX="-1" #Force the adaptive timestep to be smaller than this value, it has the priority on any other option, -1 means that the option is disabled
TSSPIN="false" #If true take into account the variation (SSE only) of OmegaSpin in the adaptive timestep
TSSPINBIN="false" #If true take into account the variation (BSE only) of OmegaSpin in the adaptive timestep
TSNSSPIN="false" #If true take into account the variation of OmegaRem for NS in the adaptive timestep. It should be set to true if interested on pulsars
CKSTALLING="true" #If true check stalling stars. If the elapsed evolution time is larger than 20s an error is thrown

PSTRING="-nthreads $NTHREADS \
-list $LISTBIN \
-ibmode $IBMODE \
-tables $TABLES \
-tables_HE $TABLESHE \
-omode $OMODE \
-o $OUTPATH \
-wmode $WINDSMODE \
-tmode $TIDES \
-gwmode $GWMODE \
-rlmode $RLMODE \
-circmode $CIRCMODE \
-cemode $CEMODE \
-mixmode $MIXING \
-collmode $COLLMODE \
-hardmode $HARDMODE \
-kmode $SNORBCHANGE \
-sn_kicks $SNKICKS \
-sn_pairinstability $SNPISN \
-sn_neutrinomaloss $SNNML \
-snmode $SUPERNOVA \
-xspinmode $BHXSPIN \
-xspin_sigma_maxwell $MAXWSDXSPIN \
-xspin_bavera $BAVERAXSPIN \
-gw_tshold $GWTSHOLD \
-gw_onlybco $GWONLYBCO \
-jtrack_h_err_rel_max $JTEMAX \
-jtrack_max_dm_factor $JTMAXDM \
-jtrack_max_iteration $JTMAXITER \
-jtrack_min_dm_factor $JTMINDM \
-jtrack_tshold_dm_rel $JTDMTSHOLD \
-eddington_factor $EDDF \
-rlo_eps_nova $RLOEPSNOVA \
-rlo_f_mass_accreted $RLOMACCR \
-rlo_gamma_angmom $RLOGAM \
-rlo_stability $RLOSTABILITY \
-rlo_max_nuclearmt $RLONTMAX \
-rlo_mtstable_ms $RLOSMTMS \
-rlo_enable_collision $RLOCOLLISION \
-sn_Mchandra $MCHANDRA \
-sn_co_lower_ecsn $SNLOWECSN \
-sn_co_lower_sn $SNLOW \
-sn_co_lower_ecsn_pureHe $SNLOWECSNHE \
-sn_co_lower_sn_pureHe $SNLOWHE \
-sn_max_ns_mass $NSMAX \
-sn_min_vkick $SNMINVKICK \
-sn_kick_velocity_stdev $SNVKICKSTD \
-sn_Mremnant_average_NS $NSMASSMEAN \
-sn_Mremnant_std_NS $NSMASSSTD \
-w_alpha $WALPHA \
-w_beta $WBETA \
-io_literal_phases $LITPHASES \
-log_level $LOGLEVEL \
-dtout "\"${DTOUT}\"" \
-ev_Nchunk $NCHUNK \
-ev_max_repetitions $MAXREP \
-tf $TEND \
-tini $TSTART \
-rseed $RSEED \
-rlo_QHE $SQHE \
-tabuse_envconv $TABCONV \
-tabuse_Xsup $TABXSUP \
-tabuse_inertia $TABINERTIA \
-inertiamode $INERTIAMODE \
-tabuse_rco  $TABRCO \
-tabuse_rhe  $TABRCO \
-ce_alpha $CEALPHA \
-ce_kce $CEKCE \
-ce_knce $CEKNCE \
-star_lambda $CELAM \
-star_lambda_pureHe $CELAMHE \
-star_lambda_fth $CELAMFTH \
-star_tshold_WR_envelope $WRTS \
-ns_magnetic_tscale $NSMAGTSCALE \
-ns_magnetic_mscale $NSMAGMSCALE \
-sn_compact_csi25_tshold $SNC25TS \
-sn_compact_fallback $SNCOMPFB \
-ev_naked_tshold $NAKEDTS \
-Z $Z  \
-io_logfile $LOGFILE \
-initerror_stop $INITERRSTOP \
-ev_set_maxCO $SMAXCO \
-ev_set_minHE $SMINHE \
-use_thg_hurley $THGHURLEY \
-ts_maximum_variation $TSMAXVAR \
-ts_min_dt $TSMIN \
-ts_max_dt $TSMAX \
-ts_check_spin $TSSPIN \
-ts_check_spin_bin $TSSPINBIN \
-ts_check_NSspin $TSNSSPIN \
-optimistic_scenario_hg $OPTIMISTIC \
-hard_rhoc $HARDRHOC \
-hard_sigma $HARDSIGMA \
-hard_xi $HARDXI \
-hard_kappa $HARDKAPPA \
-hard_mass_average $HARDMASS \
-ev_setwM $INTW \
-ev_setwM_log $INTWLOG \
-ev_setwM_tphase $INTWPHASE \
-check_stalling $CKSTALLING"

#Check SCOL, BCOL and NAMEPREX, if they are empty do not add it,
#otherwise they will generate errors (SEVN expect a non-empty value after an option)
#the [[:blank:]] is a command to remove all the white spaces.

if [ "${SCOL//[[:blank:]]/}" != ""   ]; then
  PSTRING="$PSTRING -scol $SCOL"
fi

if [ "${BCOL//[[:blank:]]/}" != "" ]; then
  PSTRING="$PSTRING -bcol $BCOL"
fi

if [ "${NAMEPREX//[[:blank:]]/}" != "" ]; then
   PSTRING="$PSTRING -name_prefix $NAMEPREX"
fi

LAUNCH="$SEVNEXE $PSTRING"
DATETIME=$(date '+%^B %d, %Y @ %H:%M:%S');
mkdir -p $OUTPATH

#Copy the executable that will be used to run
cp ${SEVNEXE}  ${OUTPATH}


#This is where we have the actual launch command
echo -e "${GREEN}SEVN is running...${NC}"
eval $LAUNCH
retVal=$?

if [ $retVal -eq 0 ]; then
    echo -e "${GREEN}SEVN execution was SUCCESSFUL!"
    echo -e "The output can be found in $OUTPATH${NC}"
else
    echo -e "${RED}SEVN execution has FAILED!${NC}"
    exit 1
fi

echo "$DATETIME $HOSTNAME [SEVN launch string]: $LAUNCH" > $OUTPATH/launch_line.txt
