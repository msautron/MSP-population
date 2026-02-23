import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d
import collections
from scipy.stats import kstest,ks_2samp,linregress
from scipy.optimize import curve_fit
import pandas as pd
from astropy.table import Table
from matplotlib.colors import LogNorm

#Dictionnary of status
dict_star_stage = {
        0:"Pre-Main Sequence",
        1:"Main Sequence",
        2:"Terminal Main Sequence",
        3:"Shell H Burning",
        4:"Core He Burning",
        5:"Terminal Core He Burning",
        6:"Shell He Burning",
        7:"Remnant"
}

dict_remnant = {
        0:"NotARemnant",
        1:"HeWD",
        2:"COWD",
        3:"ONeWD",
        4:"NS ECSN",
        5:"NS CCSN",
        6:"BH",
        -1:"Isolated"
}

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota,P_orba,Ecca,type_compa=[],[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out (all) 
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2,P_orb2,Ecc2,type_comp2=[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (all)
Pa3,P_dota3,da3,za3,xa3,ya3,agea3,E_dota3,P_orba3,Ecca3=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out (radio)
P3,P_dot3,d3,z3,x3,y3,age3,E_dot3,latitude3,longitude3,P_orb3,Ecc3,type_comp3=[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio)
Pa4,P_dota4,da4,za4,xa4,ya4,agea4,E_dota4,P_orba4,Ecca4=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out (gamma)
P4,P_dot4,d4,z4,x4,y4,age4,E_dot4,latitude4,longitude4,P_orb4,Ecc4,type_comp4=[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (gamma)
Pa5,P_dota5,da5,za5,xa5,ya5,agea5,E_dota5,P_orba5,Ecca5=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out (radio-gamma)
P5,P_dot5,d5,z5,x5,y5,age5,E_dot5,latitude5,longitude5,P_orb5,Ecc5,type_comp5=[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio-gamma)
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
G=6.67430e-11
test_l=[]
Ba,B2,Ba3,Ba4,Ba5,B3,B4,B5=[],[],[],[],[],[],[],[]
B0=3.2e15 #constant to compute the surface magnetic field of the observations
P_selected,Pdot_selected=[],[]
#B_init=[]

#Put the data at the right place
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("[-+]?\d*[.]\d*[Ee]*[-+]*\d*")
reg_4 = re.compile(r"-?\d+\.\d+")
reg_5=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
reg_6=re.compile("-?\d+[.]\d+")
reg_comp=re.compile(r'\b[A-Z][a-zA-Z]{0,3}\b')
reg_survey=re.compile(r"\b(?!NULL\b)[\w,]{3,}\b")
reg_5_2=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}|\d+")
reg_comp_mass=re.compile(r"-?\d+\.\d+(?:[eE][+-]?\d+)?")

with open("data_ATNF_fastfermipmps.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("data_ATNF_fastfermipmps.txt","r") as f:
    data2_comp=re.findall(reg_comp,f.read())

with open("data_ATNF_fastfermipmps.txt","r") as f:
    data_survey=re.findall(reg_survey,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data=re.findall(reg_1,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data_type=re.findall(reg_2,f.read())

with open("data_taud_detec.txt","r") as f: 
    data_taud=re.findall(reg_5,f.read())

with open("wr.txt","r") as f:
    data_wr=re.findall(reg_5,f.read())

with open("data_beta_detec_all.txt","r") as f:
    data_beta=re.findall(reg_5,f.read())

with open("data_beta_detec_init_all.txt","r") as f:
    data_beta_init=re.findall(reg_5,f.read())

with open("data_beta_all.txt","r") as f:
    data_beta2=re.findall(reg_5,f.read())

with open("gamma_peak_sep.txt","r") as f:
    g_peak_sep_data=re.findall(reg_5,f.read())

with open("w10_atnf.txt","r") as f:
    w10_atnf=re.findall(reg_6,f.read())

with open("mpulsar_obs.txt","r") as f:
    mpulsar_obs=re.findall(reg_6,f.read())

with open("xi_rho_data.txt","r") as f:
    xi_rho_data=re.findall(reg_5,f.read())

with open("ATNF_canonical_pulsars.txt","r") as f:
    data_cano=re.findall(reg_3,f.read())

with open("gl_gb_ATNF.txt","r") as f:
    data_gl_gb_ATNF=re.findall(reg_6,f.read())

with open("pos_init.txt","r") as f:
    data_pos_init=re.findall(reg_5,f.read())

with open("save_coord.txt","r") as f:
    data_x_y=re.findall(reg_5,f.read())

with open("pos_all.txt","r") as f:
    data_pos_all=re.findall(reg_5,f.read())

with open("subms_info_complete.txt","r") as f:
    data_subms=re.findall(reg_5_2,f.read())

with open("wint.txt","r") as f:
    data_wint=re.findall(reg_5,f.read())

with open("med_comp_mass_ATNF.txt","r") as f:
    data_comp_mass_obs=re.findall(reg_comp_mass,f.read())

with open("AIC_info.txt","r") as f:
    AIC_data=re.findall(reg_5,f.read())

#Get the data of the initial position of all pulsars
#x_pulsar,y_pulsar=[],[]
#for i in range(int(len(data_x_y)/2)):
#    x_pulsar.append(float(data_x_y[2*i]))
#    y_pulsar.append(float(data_x_y[2*i+1]))

#plt.figure(356)
#plt.xlim(-20,20)
#plt.ylim(-20,20)
#plt.scatter(x_pulsar[:100000],y_pulsar[:100000],c='blue',marker='*',s=0.0000001,label='positions of pulsars')
#plt.axis('equal')
#plt.xlabel('x (kpc)')
#plt.ylabel('y (kpc)')
#plt.savefig("initial_positions_pulsars.pdf",dpi=300)
#plt.close()

#Get the median mass of companion observed
comp_mass_obs=[]
for i in range(int(len(data_comp_mass_obs))):
    comp_mass_obs.append(float(data_comp_mass_obs[i]))

#Get the info about the subms pulsars
P_subms_not_sorted,Pdot_subms_not_sorted,Mi_subms_not_sorted,Mf_subms_not_sorted,detec_r_not_sorted,detec_g_not_sorted,detec_rg_not_sorted,subms_or_not_not_sorted=[],[],[],[],[],[],[],[]
for i in range(int(len(data_subms)/8)):
    P_subms_not_sorted.append(float(data_subms[8*i+0]))
    Pdot_subms_not_sorted.append(float(data_subms[8*i+1]))
    Mf_subms_not_sorted.append(float(data_subms[8*i+2]))
    Mi_subms_not_sorted.append(float(data_subms[8*i+3]))
    subms_or_not_not_sorted.append(float(data_subms[8*i+4]))
    detec_r_not_sorted.append(float(data_subms[8*i+5]))
    detec_g_not_sorted.append(float(data_subms[8*i+6]))
    detec_rg_not_sorted.append(float(data_subms[8*i+7]))

P_subms_sorted,Pdot_subms_sorted,Mi_subms_sorted,Mf_subms_sorted,subms_or_not_sorted=[],[],[],[],[]
for i in range(len(P_subms_not_sorted)):
    if (detec_r_not_sorted[i]==1 or detec_g_not_sorted[i]==1 or detec_rg_not_sorted[i]==1) and subms_or_not_not_sorted[i]==1:
        P_subms_sorted.append(P_subms_not_sorted[i])
        Pdot_subms_sorted.append(Pdot_subms_not_sorted[i])
        Mi_subms_sorted.append(Mi_subms_not_sorted[i])
        Mf_subms_sorted.append(Mf_subms_not_sorted[i])
        subms_or_not_sorted.append(subms_or_not_not_sorted[i])

sum_subms_detec=sum(subms_or_not_sorted)
sum_subms_nodetec=sum(subms_or_not_not_sorted)

print(f'Number of pulsars that went under the ms : {sum_subms_nodetec}\nNumber of pulsars that went under the ms and were detected : {sum_subms_detec}')

#print(P_subms_sorted)
#print(Pdot_subms_sorted)
#print(Mi_subms_sorted)
#print(Mf_subms_sorted)

#Get the position of every MSP born
x_MSP,y_MSP,z_MSP,Lg_all_MSPborn,longitude_allMSP_born,dist_rad_allMSP,dist_allMSP=[],[],[],[],[],[],[]
for i in range(int(len(data_pos_all)/5)):
    x_MSP.append(float(data_pos_all[5*i+0]))
    y_MSP.append(float(data_pos_all[5*i+1]))
    z_MSP.append(float(data_pos_all[5*i+2]))
    Lg_all_MSPborn.append(float(data_pos_all[5*i+3]))
    longitude_allMSP_born.append(float(data_pos_all[5*i+4]))
    dist_rad_allMSP.append(np.sqrt(x_MSP[i]**2+y_MSP[i]**2))
    dist_allMSP.append(np.sqrt(x_MSP[i]**2+y_MSP[i]**2+z_MSP[i]**2))

#Get the data from the ATNF catalog for longitude and latitude
latitude_ATNF,longitude_ATNF=[],[]
for i in range(int(len(data_gl_gb_ATNF)/2)):
    longitude_ATNF.append(float(data_gl_gb_ATNF[2*i]))
    latitude_ATNF.append(float(data_gl_gb_ATNF[2*i+1]))

#Get the initial positions of detected objects
x_init,y_init,z_init=[],[],[]
for i in range(int(len(data_pos_init)/3)):
    x_init.append(float(data_pos_init[3*i+0]))
    y_init.append(float(data_pos_init[3*i+1]))
    z_init.append(float(data_pos_init[3*i+2]))

#Get the data of the ATNF catalogue
P_canoa,P_dot_canoa,d_canoa,z_canoa,x_canoa,y_canoa,age_canoa,E_dot_canoa,B_canoa=[],[],[],[],[],[],[],[],[]
for i in range(int(len(data_cano)/9)):
    P_canoa+=[float(data_cano[i*9])] #Period of the rotation of the pulsar in seconds
    P_dot_canoa+=[float(data_cano[i*9+1])] #Period derivative of the rotation of the pulsar no units
    d_canoa+=[float(data_cano[i*9+2])] #Distance to us in kpc
    z_canoa+=[float(data_cano[i*9+3])] #Z position in the galactocentric frame in kpc
    x_canoa+=[float(data_cano[i*9+4])] #X position in the galactocentric frame in kpc
    y_canoa+=[float(data_cano[i*9+5])] #age of the pulsar in years
    age_canoa+=[float(data_cano[i*9+6])] #age of the pulsars in seconds
    E_dot_canoa+=[float(data_cano[i*9+7])]  #Spin down power of the pulsar in ergs/s
    B_canoa+=[float(data_cano[i*9+8])] #surface magnetic field in T

P_cano2,P_dot_cano2,d_cano2,z_cano2,x_cano2,y_cano2,age_cano2,E_dot_cano2,B_cano2=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dot_canoa)):
    if P_dot_canoa[i]!=0.0 and d_canoa[i]!=0.0 and E_dot_canoa[i]!=0 and d_canoa[i]<=25:
        P_cano2+=[P_canoa[i]]
        P_dot_cano2+=[P_dot_canoa[i]]
        d_cano2+=[d_canoa[i]]
        z_cano2+=[z_canoa[i]]
        x_cano2+=[x_canoa[i]]
        y_cano2+=[y_canoa[i]]
        age_cano2+=[age_canoa[i]]
        E_dot_cano2+=[E_dot_canoa[i]]
        B_cano2+=[B_canoa[i]]

#Get rho et xi from simulation
rho_sim,xi_sim=[],[]
for i in range(int(len(xi_rho_data)/2)):
    xi_sim.append(float(xi_rho_data[2*i]))
    rho_sim.append(float(xi_rho_data[2*i+1]))

#Get P and Porb of the saved AIC simulation
P_AIC_saved,P_orb_AIC_saved=[],[]
for i in range(int(len(AIC_data)/2)):
    P_AIC_saved.append(float(AIC_data[2*i]))
    P_orb_AIC_saved.append(float(AIC_data[2*i+1]))

#Get the data about the gamma-ray peak separation
g_peak_sep_sim=[]
for i in range(int(len(g_peak_sep_data))):
    g_peak_sep_sim.append(float(g_peak_sep_data[i]))

#Get the measured mass of pulsar in binaries
mpulsar_obs2=[]
for i in range(len(mpulsar_obs)):
    mpulsar_obs2.append(float(mpulsar_obs[i]))

#Get data from 3PC and filter the MSPs/Young pulsar
df=pd.read_excel('/home/matteo.sautron/Documents/cuda/canonical_pulsars/pulsar_population_repo/3PC_Catalog_20230803.xls')
data_3PC=Table.from_pandas(df)
data_3PC['B_S']=pd.to_numeric(data_3PC['B_S'],errors='coerce')
mask=data_3PC['B_S'] < 6e10
data_3PC_filtered=data_3PC[mask]

size_3PC_filtered=len(data_3PC_filtered['P0'])

print(f'Number of recycled pulsars in 3PC {size_3PC_filtered}')

name_to_exclude=['J0023+0923','J0340+4130','J0605+3757','J0653+4706','J1125-6014','J1142+0119','J1301+0833','J1312+0051','J1455-3330','J1544+4937','J1630+3734','J1730-2304','J1741+1351','J1745+1017','J1824-2452A','J1843-1113','J1921+0137','J1921+1929','J1946+3417','J1959+2048','J2017+0603','J2042+0246','J2215+5135','J2234+0944']

mask_exclusion = ~np.isin(data_3PC_filtered['PSRJ'], name_to_exclude)
data_3PC_filtered_complete = data_3PC_filtered[mask_exclusion]

#Get the data from tau_d file for each MSPs
taud_MSP=[]
for i in range(int(len(data_taud))):
    taud_MSP.append(float(data_taud[i]))

#Get the data from the wr file for each MSPs
wr_MSP=[]
for i in range(int(len(data_wr))):
        wr_MSP.append(float(data_wr[i]))

wint_MSP=[]
for i in range(int(len(data_wint))):
    wint_MSP.append(float(data_wint[i]))

#Get the spin-orbit angle from the simulations
spin_orb_angle=[]
for i in range(int(len(data_beta))):
    spin_orb_angle.append(float(data_beta[i]))

spin_orb_angle_all=[]
for i in range(int(len(data_beta2))):
    spin_orb_angle_all.append(float(data_beta2[i]))
    spin_orb_angle_all[i]=min(spin_orb_angle_all[i]%360,360-spin_orb_angle_all[i]%360)

spin_orb_angle_init=[]
for i in range(int(len(data_beta_init))):
    spin_orb_angle_init.append(float(data_beta_init[i]))
    spin_orb_angle_init[i]=min(spin_orb_angle_init[i]%360,360-spin_orb_angle_init[i]%360)

#List to store the binary info
Mc_i,Mc_f,a_bin,ecc_bin,comp_type,comp_rem_type=[],[],[],[],[],[]
stellar_phase_str,remnant_type_str=[],[]

#Get the data of the ATNF catalogue
for i in range(len(data2_comp)):
    type_comp2+=[data2_comp[i]]

for i in range(int(len(data2)/11)):
    Pa+=[float(data2[i*11])] #Period of the rotation of the pulsar in seconds
    P_dota+=[float(data2[i*11+1])] #Period derivative of the rotation of the pulsar no units
    P_orba+=[float(data2[i*11+2])] #Orbital period of the binary
    Ecca+=[float(data2[i*11+3])] #Eccentricity of the binary 
    da+=[float(data2[i*11+4])] #Distance to us in kpc
    za+=[float(data2[i*11+5])] #Z position in the galactocentric frame in kpc
    xa+=[float(data2[i*11+6])] #X position in the galactocentric frame in kpc
    ya+=[float(data2[i*11+7])] #age of the pulsar in years
    agea+=[float(data2[i*11+8])] #age of the pulsars in seconds
    Ba+=[float(data2[i*11+9])*1e-4] #surface magnetic field in T
    E_dota+=[float(data2[i*11+10])*1e-7] #Spin down power of the pulsar in W

for i in range(len(P_dota)):
    P2+=[Pa[i]]
    P_dot2+=[P_dota[i]]
    d2+=[da[i]]
    z2+=[za[i]]
    x2+=[xa[i]]
    y2+=[ya[i]]
    age2+=[agea[i]]
    E_dot2+=[E_dota[i]]
    B2+=[Ba[i]]
    if (Ecca[i]!=0 or (Ecca[i]==0 and (Pa[i]==0.0025731519721683 or Pa[i]==0.00621853194840048 or Pa[i]==0.00248392757283615 or Pa[i]==0.00323373735406827 or Pa[i]==0.00229472709570208 or Pa[i]==0.00335433608290620 or Pa[i]==0.00160740168480632))):
        Ecc2+=[Ecca[i]]
    elif (Ecca[i]==0 and (Pa[i]!=0.0025731519721683 or Pa[i]!=0.00621853194840048 or Pa[i]!=0.00248392757283615 or Pa[i]!=0.00323373735406827 or Pa[i]!=0.00229472709570208 or Pa[i]!=0.00335433608290620 or Pa[i]!=0.00160740168480632)):
        Ecc2+=[-1]
    if P_orba[i]!=0:
        P_orb2+=[P_orba[i]]
    elif P_orba[i]==0:
        P_orb2+=[-1]

#Computation of latitude and longitude for the ATNF data
for i in range(len(z2)):
    if d2[i]<1e-15:
        lat=0
        latitude2+=[lat]
    else:
        lat=np.arcsin(z2[i]/d2[i])*RAD
        latitude2+=[lat]
    r=(x2[i]**2+y2[i]**2)**0.5
    if x2[i]>=0:
        longi=np.arccos(-y2[i]/r)*RAD
        longitude2+=[longi]
    else:
        longi=(np.arccos(-y2[i]/r)+np.pi)*RAD
        longitude2+=[longi]

#Handle the type of pulsar
type_pulsar_obs=[]
for i in range(len(data_survey)):
    if ("fast" in data_survey[i] or "pksmb" in data_survey[i] or "htru_pks" in data_survey[i]) and not "Fermi" in data_survey[i]:
        type_pulsar_obs.append(1) #Pulsar radio
    elif not ("fast" in data_survey[i] and not "pksmb" in data_survey[i] and not "htru_pks" in data_survey[i]) and "Fermi" in data_survey[i]:
        type_pulsar_obs.append(2) #Pulsar gamma
    elif ("fast" in data_survey[i] or "pksmb" in data_survey[i] or "htru_pks" in data_survey[i]) and "Fermi" in data_survey[i]:
        type_pulsar_obs.append(3) #Pulsar rg

#print(len(P2))
#print(len(type_pulsar_obs))
#print(w10_atnf)
w10_r,w10_rg=[],[]
for i in range(len(type_pulsar_obs)):
    if type_pulsar_obs[i]==1:
        P3+=[P2[i]]
        P_dot3+=[P_dot2[i]]
        d3+=[d2[i]]
        z3+=[z2[i]]
        x3+=[x2[i]]
        y3+=[y2[i]]
        age3+=[age2[i]]
        E_dot3+=[E_dot2[i]]
        B3+=[B2[i]]
        Ecc3+=[Ecc2[i]]
        P_orb3+=[P_orb2[i]]
        latitude3+=[latitude2[i]]
        longitude3+=[longitude2[i]]
        w10_r+=[(float(w10_atnf[i])*(1e-3)*(360/P2[i]))%360]
    elif type_pulsar_obs[i]==2:
        P4+=[P2[i]]
        P_dot4+=[P_dot2[i]]
        d4+=[d2[i]]
        z4+=[z2[i]]
        x4+=[x2[i]]
        y4+=[y2[i]]
        age4+=[age2[i]]
        E_dot4+=[E_dot2[i]]
        B4+=[B2[i]]
        Ecc4+=[Ecc2[i]]
        P_orb4+=[P_orb2[i]]
        latitude4+=[latitude2[i]]
        longitude4+=[longitude2[i]]
    elif type_pulsar_obs[i]==3:
        P5+=[P2[i]]
        P_dot5+=[P_dot2[i]]
        d5+=[d2[i]]
        z5+=[z2[i]]
        x5+=[x2[i]]
        y5+=[y2[i]]
        age5+=[age2[i]]
        E_dot5+=[E_dot2[i]]
        B5+=[B2[i]]
        Ecc5+=[Ecc2[i]]
        P_orb5+=[P_orb2[i]]
        latitude5+=[latitude2[i]]
        longitude5+=[longitude2[i]]
        w10_rg+=[(float(w10_atnf[i])*(1e-3)*(360/P2[i]))%360]

print(f'Number of observed radio pulsars (FAST GPPS + PMPS): {len(P3)}')
print(f'Number of observed gamma pulsars (Fermi): {len(P4)}')

#Log calculation for age, P and Pdot for the ATNF data
log_age_selec=[]
for i in range(len(age2)):
    log_age2+=[(np.log(age2[i]))/(np.log(10))]

log_Pa=[]
for i in range(len(Pa)):
    log_Pa+=[np.log10(Pa[i])]

for i in range(len(P2)):
    log_P2+=[(np.log(P2[i]))/(np.log(10))]
    log_Pdot2+=[(np.log(P_dot2[i]))/(np.log(10))]

#get the data of the P_Pdot_positions file (simulation)
for i in range(len(data_type)):
    var2+=data_type[i][1]
    type_pulsar+=[int(var2)]
    var2=''


for i in range(len(data)):
    for j in range(len(data[i])-1):
        var+=data[i][j]
    data[i]=float(var)
    var=''

t_acc,spin_up_info,M_NS_i,M_NS_f,NS_ell=[],[],[],[],[]
for i in range(int(len(data)/31)):
    P+=[data[31*i]]    #Period of rotation of the pulsar in seconds
    P_dot+=[(data[31*i+1])]  #Period derivative of rotation of the pulsar, no units 
    x+=[data[31*i+2]]#Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y+=[data[31*i+3]]#Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    age+=[data[31*i+4]] #Age of the pulsar in seconds 
    error+=[data[31*i+5]] #error on the positions of the pulsars (error computed with the energy)
    distance+=[data[31*i+6]] #Distance to the galactic center in kpc
    longitude+=[data[31*i+7]] #galactic latitude in degrees
    latitude+=[data[31*i+8]] #galactic longitude in degrees
    cos_alpha0+=[np.cos(data[31*i+9])] #cosinus of the initial inclination angle
    cos_alpha+=[np.cos(data[31*i+10])] #cosinus of the inclination angle after all the evolution
    Bf+=[data[31*i+11]] #Magnetic field of the pulsar today, in tesla
    z+=[data[31*i+12]] #Position on the z-absciss relative to the sun in the galactocentric frame in kpc
    vx+=[data[31*i+13]] #Velocity on the x-absciss of the pulsar today, in km/s
    vy+=[data[31*i+14]] #Velocity on the y-absciss of the pulsar today, in km/s
    vz+=[data[31*i+15]] #Velocity on the z-absciss of the pulsar today, in km/s
    vx0+=[data[31*i+16]] #Velocity on the x-absciss of the pulsar initially, in km/s
    vy0+=[data[31*i+17]] #Velocity on the y-absciss of the pulsar initially, in km/s
    vz0+=[data[31*i+18]] #Velocity on the z-absciss of the pulsar initially, in km/s
    PA+=[(180/np.pi)*np.arccos(data[31*i+19])]
    t_acc+=[data[31*i+20]] #Duration of spin up episode in s
    spin_up_info+=[data[31*i+21]] #=1 if spin up happened =0 otherwise
    M_NS_i+=[data[31*i+22]] #Initial mass of the NS in solar mass
    M_NS_f+=[data[31*i+23]] #Final mass of the NS in solar mass
    Mc_i+=[data[31*i+24]] #Initial mass of the companion star in solar mass
    Mc_f+=[data[31*i+25]] #Final mass of the companion star in solar mass
    a_bin+=[data[31*i+26]] #Final semi major axis in Rsun
    ecc_bin+=[data[31*i+27]] #Eccentricity of the binary
    comp_type+=[int(data[31*i+28])] #Type of the companion if not a remnant (int)
    comp_rem_type+=[int(data[31*i+29])] #Type of the companion if it is a remnant (int)
    NS_ell+=[data[31*i+30]] #Ellipticity of the NS

print(f'Check recovery of the data: {P[-1]}\n')

num_spin_up=0
P_spinup,P_dot_spinup,x_spinup,y_spinup,age_spinup,error_spinup,distance_spinup,longitude_spinup,latitude_spinup,cos_alpha0_spinup,cos_alpha_spinup=[],[],[],[],[],[],[],[],[],[],[]
Bf_spinup,z_spinup,vx_spinup,vy_spinup,vz_spinup,vx0_spinup,vy0_spinup,vz0_spinup,PA_spinup,t_acc_spinup,M_NS_i_spinup,M_NS_f_spinup,Mc_i_spinup,Mc_f_spinup,a_bin_spinup,ecc_bin_spinup,comp_type_spinup,comp_rem_type_spinup=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
type_pulsar_spinup,NS_ell_spinup,taud_spinup,wr_spinup=[],[],[],[]
spin_orb_angle_spinup,spin_orb_angle_init_spinup=[],[]
xi_sim_spinup,rho_sim_spinup=[],[]
w_geometry_spinup=[]
for i in range(len(spin_up_info)):
    if spin_up_info[i]==1 and P[i]>1e-3:
        P_spinup.append(P[i])
        P_dot_spinup.append(P_dot[i])
        x_spinup.append(x[i])
        y_spinup.append(y[i])
        age_spinup.append(age[i])
        error_spinup.append(error[i])
        distance_spinup.append(distance[i])
        longitude_spinup.append(longitude[i])
        latitude_spinup.append(latitude[i])
        cos_alpha0_spinup.append(cos_alpha0[i])
        cos_alpha_spinup.append(cos_alpha[i])
        Bf_spinup.append(Bf[i])
        z_spinup.append(z[i])
        vx_spinup.append(vx[i])
        vy_spinup.append(vy[i])
        vz_spinup.append(vz[i])
        vx0_spinup.append(vx0[i])
        vy0_spinup.append(vy0[i])
        vz0_spinup.append(vz0[i])
        PA_spinup.append(PA[i])
        t_acc_spinup.append(t_acc[i])
        M_NS_i_spinup.append(M_NS_i[i])
        M_NS_f_spinup.append(M_NS_f[i])
        Mc_i_spinup.append(Mc_i[i])
        Mc_f_spinup.append(Mc_f[i])
        a_bin_spinup.append(a_bin[i])
        ecc_bin_spinup.append(ecc_bin[i])
        comp_type_spinup.append(comp_type[i])
        comp_rem_type_spinup.append(comp_rem_type[i])
        type_pulsar_spinup.append(type_pulsar[i])
        NS_ell_spinup.append(NS_ell[i])
        taud_spinup.append(taud_MSP[i])
        wr_spinup.append(wr_MSP[i])
        w_geometry_spinup.append(wint_MSP[i])
        spin_orb_angle_spinup.append(min(spin_orb_angle[i]%360,360-spin_orb_angle[i]%360))
        spin_orb_angle_init_spinup.append(spin_orb_angle_init[i])
        xi_sim_spinup.append(180*xi_sim[i]/np.pi)
        rho_sim_spinup.append(180*rho_sim[i]/np.pi)
        num_spin_up+=1

#Check the r_KG/r_L
r_kg_spinup,r_L_spinup,ratio_rkG_rL=[],[],[]
for i in range(len(rho_sim_spinup)):
    r_kg=((rho_sim_spinup[i]/1.24)**2*P_spinup[i])*12000
    r_L=((3e8*P_spinup[i])/(2*np.pi))
    r_kg_spinup.append(r_kg)
    r_L_spinup.append(r_L)
    ratio_rkG_rL.append(r_kg/r_L)

#count_nb_belowPdot=0
#for i in range(len(P_spinup)):
#    if P_dot_spinup[i]<2e-17:
#        count_nb_belowPdot+=1

#print(f'Number of ONeWD below Pdot = 2e-17 : {count_nb_belowPdot}\n')

#Mass of NS analysis
three_biggest=sorted(M_NS_f_spinup, reverse=True)[:10]
print(f'Mean mass NS : {np.mean(M_NS_f_spinup)}\n')
print(f'Maximum mass of NS : {three_biggest}\n')

#Count number of aligned spin-orbit angle
count_below10deg=sum(angle < 10 for angle in spin_orb_angle_spinup)
prct_below10deg=100*(count_below10deg/len(spin_orb_angle_spinup))
print(f'Percentage of binaries with the spin orbit angle below 10 ° : {prct_below10deg} %')

#Compute gamma flux to compare to observations
flux_gamma_all=[]
for i in range(len(P_spinup)):
    #if (type_pulsar[i]>1):
    Edot_spinup=4*np.pi**2*1e38*P_dot_spinup[i]*(P_spinup[i])**(-3)
    fgval=(10**(26.4)*((Bf_spinup[i]/1e8)**(0.11))*((Edot_spinup/1e26)**(0.51)))/(4*np.pi*(distance_spinup[i]*3.086e19)**2)
    fgval=(10**(26.4)*((Bf_spinup[i]/1e8)**(0.11))*((Edot_spinup/1e26)**(0.51)))
    flux_gamma_all.append(fgval)

#flux_gamma_all_sorted=sorted(flux_gamma_all)[:20]
#print(f'Minimum flux detected among gamma pulsars: {flux_gamma_all_sorted}') 
#print(f'Minimum Lum all gamma pulsar detected:')
#print(np.min(flux_gamma_all))
#print(f'Maximum Lum all gamma pulsar detected:')
#print(np.max(flux_gamma_all))
#print(f'Mean Lum all gamma pulsar detected:')
#print(np.mean(flux_gamma_all))

P_spinup2,P_dot_spinup2,x_spinup2,y_spinup2,age_spinup2,error_spinup2,distance_spinup2,longitude_spinup2,latitude_spinup2,cos_alpha0_spinup2,cos_alpha_spinup2=[],[],[],[],[],[],[],[],[],[],[]
Bf_spinup2,z_spinup2,vx_spinup2,vy_spinup2,vz_spinup2,vx0_spinup2,vy0_spinup2,vz0_spinup2,PA_spinup2,t_acc_spinup2,M_NS_i_spinup2,M_NS_f_spinup2,Mc_i_spinup2,Mc_f_spinup2,a_bin_spinup2,ecc_bin_spinup2,comp_type_spinup2,comp_rem_type_spinup2=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
type_pulsar_spinup2,NS_ell_spinup2,taud_spinup2,wr_spinup2=[],[],[],[]
g_peak_sep_sim2=[]
spin_orb_angle_spinup2=[]
for i in range(len(spin_up_info)):
    if spin_up_info[i]==1 and P[i]>1e-3:
        if (type_pulsar[i]==1) or (type_pulsar[i]==3) or (type_pulsar[i]==2 and (P[i]<=6e-3 and P_dot[i]<=1e-19)):
            P_spinup2.append(P[i])
            P_dot_spinup2.append(P_dot[i])
            x_spinup2.append(x[i])
            y_spinup2.append(y[i])
            age_spinup2.append(age[i])
            error_spinup2.append(error[i])
            distance_spinup2.append(distance[i])
            longitude_spinup2.append(longitude[i])
            latitude_spinup2.append(latitude[i])
            cos_alpha0_spinup2.append(cos_alpha0[i])
            cos_alpha_spinup2.append(cos_alpha[i])
            Bf_spinup2.append(Bf[i])
            z_spinup2.append(z[i])
            vx_spinup2.append(vx[i])
            vy_spinup2.append(vy[i])
            vz_spinup2.append(vz[i])
            vx0_spinup2.append(vx0[i])
            vy0_spinup2.append(vy0[i])
            vz0_spinup2.append(vz0[i])
            PA_spinup2.append(PA[i])
            t_acc_spinup2.append(t_acc[i])
            M_NS_i_spinup2.append(M_NS_i[i])
            M_NS_f_spinup2.append(M_NS_f[i])
            Mc_i_spinup2.append(Mc_i[i])
            Mc_f_spinup2.append(Mc_f[i])
            a_bin_spinup2.append(a_bin[i])
            ecc_bin_spinup2.append(ecc_bin[i])
            comp_type_spinup2.append(comp_type[i])
            comp_rem_type_spinup2.append(comp_rem_type[i])
            type_pulsar_spinup2.append(type_pulsar[i])
            NS_ell_spinup2.append(NS_ell[i])
            taud_spinup2.append(taud_MSP[i])
            wr_spinup2.append(wr_MSP[i])
            spin_orb_angle_spinup2.append(min(spin_orb_angle[i]%360,360-spin_orb_angle[i]%360))
            g_peak_sep_sim2.append(g_peak_sep_sim[i])

count_spinup=0
for i in range(len(spin_up_info)):
    if spin_up_info[i]==1:
        count_spinup+=1
print(f'Number of Pulsars that got spun up : {count_spinup}')

#Get the initial mass of companion whose pulsars ended up with P>10^-1 s
Mc_bigP,Mc_lowP=[],[]
for i in range(len(Mc_i_spinup)):
    if P_spinup[i]>=0.02:
        Mc_bigP.append(Mc_i_spinup[i])
    else:
        Mc_lowP.append(Mc_i_spinup[i])


P_orb=[]
RSUN=6.957e8
MSUN=1.98847e30
for i in range(len(P_spinup)):
    if comp_rem_type_spinup[i]==-1:
        P_orb.append(0)
    else: 
        P_orb_val=(np.sqrt(4*((np.pi)**2)*(a_bin_spinup[i]*RSUN)**3/(G*MSUN*(Mc_f_spinup[i]+M_NS_f_spinup[i]))))/(24*60*60) #Orbital period of the binary in days
        P_orb.append(P_orb_val)
        #if comp_rem_type_spinup[i]==4 or comp_rem_type_spinup[i]==5:
            #print(f'The charac : {a_bin_spinup[i]}, {Mc_f_spinup[i]}, {M_NS_f_spinup[i]}\n')
            #print(f'It has a Porb of : {P_orb_val}\n')

#Get the mass of WD companions at the end of simulation
M_WD_only=[]
for i in range(len(P_spinup)):
    if comp_rem_type_spinup[i]==1 or comp_rem_type_spinup[i]==2 or comp_rem_type_spinup[i]==3:
        M_WD_only.append(Mc_f_spinup[i])

M_NS_accr,Mdot_mean=[],[]
for i in range(len(M_NS_i_spinup)):
    M_NS_accr.append(M_NS_f_spinup[i]-M_NS_i_spinup[i]) #Solar mass
    if t_acc_spinup[i]!=0:
        Mdot_mean.append(M_NS_accr[i]/(t_acc_spinup[i]/(365*24*3600))) #Solar mass per year
    else:
        Mdot_mean.append(0.0)

for num in comp_type_spinup:
    if num in dict_star_stage:
        stellar_phase_str.append(dict_star_stage[num])

for num in comp_rem_type_spinup:
    if num in dict_remnant:
        remnant_type_str.append(dict_remnant[num])

stellar_phase_str2,remnant_type_str2=[],[]

for num in comp_type_spinup2:
    if num in dict_star_stage:
        stellar_phase_str2.append(dict_star_stage[num])

for num in comp_rem_type_spinup2:
    if num in dict_remnant:
        remnant_type_str2.append(dict_remnant[num])

print(f"Number of pulsars that had a spin-up episode : {num_spin_up}")

#Compute the number of possible transients MSPs formula from Cui & Li (2024), taken from Dubus et al. (1999) 
M_dot_cri_trans,P_trans,Pdot_trans=[],[],[]
count_transients=0
for i in range(len(M_NS_f_spinup)):
    if P_orb[i]!=0 and Mc_f_spinup[i]!=-1:
        M_dot_cri_trans.append(3.2e-9*(((M_NS_f_spinup[i]/1.4)**0.5)*((Mc_f_spinup[i])**(-0.2))*((P_orb[i])**1.4)))
        if (M_dot_cri_trans[i]<(Mdot_mean[i]/0.5)):
            count_transients+=1
            P_trans.append(P_spinup[i])
            Pdot_trans.append(P_dot_spinup[i])

    else:
        M_dot_cri_trans.append(0)

print(f'Number of possible transients MSPs : {count_transients}')

c_light=2.997924858e8 #velocity of light in m/s

#Relation between M_WD and P_orb check (Tauris & Savonije (1999))
M_WD_sim_g2l100,M_WD_sim_g100,M_WD_rel1,M_WD_rel2=[],[],[],[]
ecc_rel1,ecc_rel2=[],[]
a_wd=4.75
b_wd=1.1e5
c_wd=0.115
count_below_5pc,count_below_10pc,count_below_20pc=0,0,0
for i in range(len(P_spinup)):
    if remnant_type_str[i]=='HeWD' and Mc_f[i]>0.18 and Mc_f[i]<0.46 and P_orb[i]>=2 and P_orb[i]<=100:
        M_WD_sim_g2l100.append(Mc_f[i])
        val_int_mwd=((P_orb[i]/b_wd)**(1.0/a_wd))+c_wd
        M_WD_rel1+=[val_int_mwd] #Porb must be in days and M_WD will be in solar mass
        ecc_rel1.append((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd))
        if ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.05):
            count_below_5pc+=1
            count_below_10pc+=1
            count_below_20pc+=1
        elif ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.1):
            count_below_10pc+=1
            count_below_20pc+=1
        elif ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.2):
            count_below_20pc+=1
    elif remnant_type_str[i]=='HeWD' and Mc_f[i]>0.18 and Mc_f[i]<0.46 and P_orb[i]>100:
        M_WD_sim_g100.append(Mc_f[i])
        val_int_mwd=((P_orb[i]/b_wd)**(1.0/a_wd))+c_wd
        M_WD_rel2+=[val_int_mwd] #Porb must be in days and M_WD will be in solar mass
        ecc_rel2.append((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd))
        if ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.05):
            count_below_5pc+=1
            count_below_10pc+=1
            count_below_20pc+=1
        elif ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.1):
            count_below_10pc+=1
            count_below_20pc+=1
        elif ((np.abs(Mc_f[i]-val_int_mwd)/val_int_mwd)<0.2):
            count_below_20pc+=1


#print(f'Percentage of He WD companion of NS with a relative error of less than 5% while computing its mass with the formula of Tauris & Savonije (1999) : {(100*count_below_5pc)/(len(ecc_rel1)+len(ecc_rel2))}%')
#print(f'Percentage of He WD companion of NS with a relative error of less than 10% while computing its mass with the formula of Tauris & Savonije (1999) : {(100*count_below_10pc)/(len(ecc_rel1)+len(ecc_rel2))}%')
#print(f'Percentage of He WD companion of NS with a relative error of less than 20% while computing its mass with the formula of Tauris & Savonije (1999) : {(100*count_below_20pc)/(len(ecc_rel1)+len(ecc_rel2))}%')

#Edot computation 
Edot=[]
Inertia=1e38
count_rad,count_gam,count_radgam,count_radgam_E_big,count_gam_E_big,count_rad_E_big,count_radgam_E_bigbig,count_gam_E_bigbig,count_rad_E_bigbig=0,0,0,0,0,0,0,0,0
for i in range(len(P_spinup)):
    Edot+=[4*np.pi**2*Inertia*P_dot_spinup[i]*(P_spinup[i]**(-3))]

for i in range(0,len(age_spinup)):
    logage=(np.log(age_spinup[i]/(365*24*60*60)))/(np.log(10))
    log_age+=[logage]

log_charac_age=[]
for i in range(len(P_spinup)):
    log_P+=[(np.log(P_spinup[i]))/(np.log(10))]
    log_Pdot+=[(np.log(P_dot_spinup[i]))/(np.log(10))]
    log_charac_age.append(np.log10((P_spinup[i]/(2*P_dot_spinup[i]))/(365*24*3600)))

#Prep death line Ruderman & Sutherland 1975
R_NS=12000
mu_0=1.25663706212e-6 
P_dot_death,P_dot_death2=[],[]
P_death=[np.log10(i) for i in np.arange(1e-2,3e1,0.01)]
P_death2=[10**(P_death[i]) for i in range(len(P_death))]
for i in range(len(P_death)):
    P_dot_death+=[3*P_death[i]+np.log10(((16*(np.pi**3)*(R_NS**6)*(1+((np.sin(45*np.pi/180)**2))))*(17000000**2))/(Inertia*mu_0*(c_light**3)))]
    P_dot_death2+=[10**(P_dot_death[i])]

#B(P,Pdot) computation in order to compare with the decaying Bf
B_ppdot,v_norm,v0_norm,err_rel_B=[],[],[],[]
for i in range(len(P_spinup)):
    B_ppdot+=[((Inertia*mu_0*(c_light**3)*P_dot_spinup[i]*P_spinup[i])/(16*(np.pi**3)*(R_NS**6)*(1+(np.sin(np.arccos(cos_alpha_spinup[i]))**2))))**0.5]
    v_norm+=[(vx_spinup[i]**2+vy_spinup[i]**2+vz_spinup[i]**2)**0.5]
    v0_norm+=[(vx0_spinup[i]**2+vy0_spinup[i]**2+vz0_spinup[i]**2)**0.5]
    err_rel_B+=[(np.abs(B_ppdot[i]-Bf_spinup[i]))/np.abs(B_ppdot[i])]

#Prep plot period old pulsars and plot old pulsars spin-velocity angle
P_old=[]
PA_old,PA_young=[],[]
period_old,period_young,pdot_young,pdot_old=[],[],[],[]
for i in range(len(log_age)):
    if log_age[i]>7.5 and log_age[i]<9:
        P_old+=[log_P[i]]

for i in range(len(PA_spinup)):
    if log_age[i]>=9.3:
        PA_old.append(PA_spinup[i])
        period_old.append(P_spinup[i])
        pdot_old.append(P_dot_spinup[i])
    elif log_age[i]<9.3:
        PA_young.append(PA_spinup[i])
        period_young.append(P_spinup[i])
        pdot_young.append(P_dot_spinup[i])


#Lists depending on the pulsar emission type
P_radio,P_dot_radio,x_radio,y_radio,age_radio,error_radio,distance_radio=[],[],[],[],[],[],[]
P_gamma,P_dot_gamma,x_gamma,y_gamma,age_gamma,error_gamma,distance_gamma=[],[],[],[],[],[],[]
P_radio_gamma,P_dot_radio_gamma,x_radio_gamma,y_radio_gamma,age_radio_gamma,error_radio_gamma,distance_radio_gamma=[],[],[],[],[],[],[]
wr_r_or_rg,P_r_or_rg=[],[]
taud_r,taud_g,taud_rg=[],[],[]
flux_gamma_only=[]
xi_r,xi_g,xi_rg=[],[],[]
alpha_r,alpha_g,alpha_rg=[],[],[]
Bf_r,Bf_g,Bf_rg=[],[],[]
z_r,z_g,z_rg=[],[],[]
Fr_gamma=[]
Edot_g=0
count_gamma_fr_good=0
w_geometry_r_or_rg=[]
ratio_rkG_rL_radio,ratio_rkG_rL_gamma,ratio_rkG_rL_rg=[],[],[]
rho_r,rho_g,rho_rg=[],[],[]
r_kg_r,r_kg_g,r_kg_rg=[],[],[]
for i in range(len(P_spinup)):
    if type_pulsar_spinup[i]==1:
        P_radio+=[P_spinup[i]]
        P_dot_radio+=[P_dot_spinup[i]]
        x_radio+=[x_spinup[i]]
        y_radio+=[y_spinup[i]]
        age_radio+=[age_spinup[i]]
        distance_radio+=[distance_spinup[i]]
        wr_r_or_rg.append(wr_spinup[i])
        w_geometry_r_or_rg.append(w_geometry_spinup[i])
        P_r_or_rg.append(P_spinup[i])
        taud_r.append(taud_spinup[i])
        xi_r.append(xi_sim_spinup[i])
        alpha_r.append(min(180*np.arccos(cos_alpha_spinup[i])/np.pi,180-180*np.arccos(cos_alpha_spinup[i])/np.pi))
        Bf_r.append(Bf_spinup[i])
        z_r.append(z_spinup[i])
        rho_r.append(rho_sim_spinup[i])
        ratio_rkG_rL_radio.append(ratio_rkG_rL[i])
        r_kg_r.append(r_kg_spinup[i]/12000)
    elif type_pulsar_spinup[i]==2:
        P_gamma+=[P_spinup[i]]
        P_dot_gamma+=[P_dot_spinup[i]]
        x_gamma+=[x_spinup[i]]
        y_gamma+=[y_spinup[i]]
        age_gamma+=[age_spinup[i]]
        distance_gamma+=[distance_spinup[i]]
        taud_g.append(taud_spinup[i])
        flux_gamma_only.append(flux_gamma_all[i])
        xi_g.append(xi_sim_spinup[i])
        alpha_g.append(min(180*np.arccos(cos_alpha_spinup[i])/np.pi,180-180*np.arccos(cos_alpha_spinup[i])/np.pi))
        Bf_g.append(Bf_spinup[i])
        z_g.append(z_spinup[i])
        Edot_g=4*np.pi**2*1e38*(P_dot_spinup[i])*(P_spinup[i])**(-3)
        Fj=np.random.normal(loc=0.0,scale=0.2)
        Fr_gamma_value=9*1e3*(distance_spinup[i])**(-2)*(Edot_g/1e29)**(0.125)*10**Fj
        Fr_gamma.append(Fr_gamma_value)
        ratio_rkG_rL_gamma.append(ratio_rkG_rL[i])
        rho_g.append(rho_sim_spinup[i])
        r_kg_g.append(r_kg_spinup[i]/12000)
        if (Fr_gamma_value>=30):
            count_gamma_fr_good+=1
    elif type_pulsar_spinup[i]==3:
        P_radio_gamma+=[P_spinup[i]]
        P_dot_radio_gamma+=[P_dot_spinup[i]]
        x_radio_gamma+=[x_spinup[i]]
        y_radio_gamma+=[y_spinup[i]]
        age_radio_gamma+=[age_spinup[i]]
        distance_radio_gamma+=[distance_spinup[i]]
        wr_r_or_rg.append(wr_spinup[i])
        w_geometry_r_or_rg.append(w_geometry_spinup[i])
        P_r_or_rg.append(P_spinup[i])
        taud_rg.append(taud_spinup[i])
        xi_rg.append(xi_sim_spinup[i])
        alpha_rg.append(min(180*np.arccos(cos_alpha_spinup[i])/np.pi,180-180*np.arccos(cos_alpha_spinup[i])/np.pi))
        Bf_rg.append(Bf_spinup[i])
        z_rg.append(z_spinup[i])
        ratio_rkG_rL_rg.append(ratio_rkG_rL[i])
        rho_rg.append(rho_sim_spinup[i])
        r_kg_rg.append(r_kg_spinup[i]/12000)

#Estimate gamma luminosity in the galactic center
L_gamma_selected=[]
for i in range(len(P_gamma)):
    if (np.sqrt(x_gamma[i]**2+y_gamma[i]**2+z_g[i]**2)<=0.7):
        E_dot=4*np.pi**2*Inertia*P_dot_gamma[i]*(P_gamma[i])**-3
        L_gamma_selected.append(10**(26.4)*(Bf_g[i]/1e8)**(0.11)*(E_dot/1e26)**(0.51))

for i in range(len(P_radio_gamma)):
    if (np.sqrt(x_radio_gamma[i]**2+y_radio_gamma[i]**2+z_rg[i]**2)<=0.7):
        E_dot=4*np.pi**2*Inertia*P_dot_radio_gamma[i]*(P_radio_gamma[i])**-3
        L_gamma_selected.append(10**(26.4)*(Bf_rg[i]/1e8)**(0.11)*(E_dot/1e26)**(0.51))

L_tot_GC=sum(L_gamma_selected)
L_mean=np.mean(L_gamma_selected)
#print(L_gamma_selected)
#print(f'In the Galactic center:')
#print(f'Min L_gamma value : {min(L_gamma_selected)}, Max L_gamma value : {max(L_gamma_selected)}')
#print(f'Total luminosity in a 0.7 kpc radius around the Galactic center : {L_tot_GC} W\nMean luminosity of a pulsar in the GC : {L_mean} W\nActual number of pulsar in the GC : {len(L_gamma_selected)}\n')

count_gpulsar=0
count_gpulsar2=0
#count_AIC_recycled=0
P_gamma_comp,P_dot_gamma_comp=[],[]
Fr_gamma_area_obs,xi_g2,rho_g2,alpha_g2=[],[],[],[]
for i in range(len(P_gamma)):
    if P_dot_gamma[i]<3e-20:
        count_gpulsar+=1
        P_gamma_comp.append(P_gamma[i])
        P_dot_gamma_comp.append(P_dot_gamma[i])
        xi_g2.append(xi_g[i])
        rho_g2.append(rho_g[i])
        alpha_g2.append(alpha_g[i])
        Fr_gamma_area_obs.append(Fr_gamma[i])
    if Bf_g[i]<=4e4:
        count_gpulsar2+=1

#for i in range(len(P_spinup)):
#    if Bf_spinup[i]<5e6:
#        count_AIC_recycled+=1

#print(f'Number of simulated WD which became recycled pulsars : {count_AIC_recycled}')

print(f'Number of simulated gamma pulsars detected below Pdot = 3.10^-20 {count_gpulsar}\n')
print(f'Number of simulated gamma pulsars detected below B = 4.10^4 T {count_gpulsar2}\n')

#Plot the death line of the article from Mitra et al. (2019)
T_6=2
T_6_max=2.8
T_6_min=1.9
eta=0.15
alpha_l=45*np.pi/180
alpha_l_max=65*np.pi/180
alpha_l_min=0*np.pi/180
b=40
b_min=60
b_max=30
const=(3.16e-4*(T_6**4)*1e-15)/((eta)**2*b*(np.cos(alpha_l))**2)
const_min=(3.16e-4*(T_6_min**4)*1e-15)/((eta)**2*b_min*(np.cos(alpha_l_min))**2)
const_max=(3.16e-4*(T_6_max**4)*1e-15)/((eta)**2*b_max*(np.cos(alpha_l_max))**2)
P_line=[10**(np.log10(i)) for i in np.arange(1e-3,1e2,0.001)]
Pdot_line=[10**np.log10(const*(i**2)) for i in np.arange(1e-3,1e2,0.001)]
Pdot_line4=[Pdot_line[i]*10**(-0.55) for i in range(len(Pdot_line))]
Pdot_line5=[Pdot_line[i]*10**(1.15) for i in range(len(Pdot_line))]
Pdot_line2=[10**np.log10(const_min*(i**2)) for i in np.arange(1e-3,1e2,0.001)]
Pdot_line3=[10**np.log10(const_max*(i**2)) for i in np.arange(1e-3,1e2,0.001)]

#Plot the death line of the article of Chen & Ruderman (1993)
const_CR93=10**((43.8-16)/2)*16*(np.pi**3)*(R_NS**6)*(1+(np.sin(alpha_l)**2))/(Inertia*mu_0*(c_light**3))
P_dot_line_CR93=[10**np.log10(const_CR93*(i**(2))) for i in np.arange(1e-3,1e2,0.001)]

#Plot the line with the critical magnetic field 4.4e9 T, distinguishing magnetar from canonical pulsars
B_crit=4.4e9
B_linecrit=[10**np.log10(((B_crit/B0)**2)*(i**(-1))) for i in np.arange(1e-3,1e2,0.001)]

#Plot the Edot lines
P_dot_Edot1e23W=[(1e23*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e28W=[(1e28*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e25W=[(1e25*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e31W=[(1e31*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]

#Plot the B lines
P_dot_B5e6=[(5e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e3=[(1e3/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e4=[(1e4/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e5=[(1e5/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e6=[(1e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e7=[(1e7/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e8=[(1e8/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]

#Display of the statistics on all pulsar
for i in range(len(Edot)):
    if Edot[i] > 1e31 and type_pulsar_spinup[i]==1:
        count_rad_E_bigbig+=1
    if Edot[i] > 1e31 and type_pulsar_spinup[i]==2:
        count_gam_E_bigbig+=1
    if Edot[i] > 1e31 and type_pulsar_spinup[i]==3:
        count_radgam_E_bigbig+=1
    if Edot[i] > 1e28 and type_pulsar_spinup[i]==1:
        count_rad_E_big+=1
    if Edot[i] > 1e28 and type_pulsar_spinup[i]==2:
        count_gam_E_big+=1
    if Edot[i] > 1e28 and type_pulsar_spinup[i]==3:
        count_radgam_E_big+=1
    if type_pulsar_spinup[i]==1:
        count_rad+=1
    if type_pulsar_spinup[i]==2:
        count_gam+=1
    if type_pulsar_spinup[i]==3:
        count_radgam+=1

print(f"Number of radio pulsars with Edot > 1e31 W : {count_rad_E_bigbig}\nNumber of gamma pulsars with Edot > 1e31 W : {count_gam_E_bigbig}")
print(f"Number of radio-gamma pulsars with Edot > 1e31 W : {count_radgam_E_bigbig}")
print(f"Number of radio pulsars with Edot > 1e28 W : {count_rad_E_big}")
print(f"Number of gamma pulsars with Edot > 1e28 W : {count_gam_E_big}")
print(f"Number of radio-gamma pulsars with Edot > 1e28 W : {count_radgam_E_big}")
print(f"Number of radio pulsars : {count_rad}")
print(f"Number of gamma pulsars : {count_gam}")
print(f"Number of radio-gamma pulsars : {count_radgam}")

count_geo_flux_r_good,count_geo_flux_r_good2=0,0
#Check both radio flux and geometry for pulsation in gamma-ray pulsars
for i in range(len(P_gamma)):
    if ((np.abs(xi_g[i]-alpha_g[i])<=rho_g[i] or np.abs(xi_g[i]-(180-alpha_g[i]))<= rho_g[i]) and Fr_gamma[i]>=30):
        count_geo_flux_r_good+=1

for i in range(len(Fr_gamma_area_obs)):
    if ((np.abs(xi_g2[i]-alpha_g2[i])<=rho_g2[i] or np.abs(xi_g2[i]-(180-alpha_g2[i]))<= rho_g2[i]) and Fr_gamma_area_obs[i]>=30):
        count_geo_flux_r_good2+=1

print(f'Number of detected gamma-ray pulsars with a radio flux greater than 30 µJy : {count_gamma_fr_good}\nRatio of gamma-ray pulsars with a radio flux greater than 30 µJy : {(count_gamma_fr_good)/(len(P_gamma))}')
print(f'Number of detected gamma-ray pulsars with a radio flux greater than 30 µJy and good geometry for radio pulsation: {count_geo_flux_r_good}\nRatio of gamma-ray pulsars with a radio flux greater than 30 µJy and good geometry for pulsation: {count_geo_flux_r_good/(len(P_gamma))}')
print(f'Number of detected gamma-ray pulsars with a radio flux greater than 30 µJy and good geometry for radio pulsation (area of observed gamma-ray pulsars only): {count_geo_flux_r_good2}\nRatio of gamma-ray pulsars with a radio flux greater than 30 µJy and good geometry for pulsation (area of observed gamma-ray pulsars only): {count_geo_flux_r_good2/(len(Fr_gamma_area_obs))}')

#Compute the gamma-ray peak separation
data_3PC_filtered['PKSEP']=pd.to_numeric(data_3PC_filtered['PKSEP'],errors='coerce')
for i in range(len(data_3PC_filtered['PKSEP'])):
    if data_3PC_filtered['PKSEP'][i] >0.5 or data_3PC_filtered['PKSEP'][i] < 0:
        data_3PC_filtered['PKSEP'][i]=1.0-data_3PC_filtered['PKSEP'][i]


#Compute the number of MSP in the Galactic center
dist_MSP_GC=[]
nb_MSP_in_GC=0
L_gamma_allMSPborn_inGC=0
for i in range(len(x_MSP)):
    dist_MSPborn=np.sqrt((x_MSP[i])**2+(y_MSP[i])**2+(z_MSP[i])**2)
    if (dist_MSPborn <= 0.7):
        L_gamma_allMSPborn_inGC+=Lg_all_MSPborn[i]
        nb_MSP_in_GC+=1

print(f'Number of MSP in the GC (detected or not) : {nb_MSP_in_GC}\nLuminosity of MSP in the GC : {L_gamma_allMSPborn_inGC}\n')

#Axisymetric distribution in the disc or not ?
theta_disc=np.arctan2(y_MSP,x_MSP)
theta_disc=(180/np.pi)*(theta_disc % (2*np.pi))
#theta_disc=[]
#for i in range(len(x_MSP)):
    #theta_disc.append((180/np.pi)*(np.arccos(x_MSP[i]/(np.sqrt((x_MSP[i])**2+(y_MSP[i])**2)))))
    #theta_disc.append((180/np.pi)*(np.arctan2(y_MSP[i]/x_MSP[i])%2*np.pi))

#KS test 1D python (all the data)
KS_test=kstest(P_dot_spinup2,P_dot2)
test_stat_all1=KS_test.statistic
p_value_all1=KS_test.pvalue
print("----ALL THE DATA----\n")
print(f"d_value of Pdot KS test = {test_stat_all1}")
print(f"p_value of Pdot={p_value_all1}")

KS_test=kstest(P_spinup2,P2)
test_stat_all2=KS_test.statistic
p_value_all2=KS_test.pvalue
print(f"d_value of P KS test = {test_stat_all2}")
print(f"p_value of P={p_value_all2}")

##KS test 1D python (gamma-ray only population)
KS_test=kstest(P_dot_gamma_comp,data_3PC_filtered_complete['P1'])
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----GAMMA ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of Pdot for the gamma only pulsars={p_value}")

KS_test=kstest(P_gamma_comp,data_3PC_filtered_complete['P0'])
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of P for the gamma only pulsars={p_value}")

#KS test 1D python (radio/gamma-ray population)
KS_test=kstest(P_dot_radio_gamma,P_dot5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO/GAMMA PULSARS----\n")
print(f"d_value of Pdot KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of Pdot for the radio/gamma pulsars={p_value}")

KS_test=kstest(P_radio_gamma,P5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of P for the radio/gamma pulsars={p_value}")

#KS test 1D python (radio population)
KS_test=kstest(P_dot_radio,P_dot3)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the radio only pulsars= {test_stat}")
print(f"p_value of Pdot for the radio only pulsars={p_value}")

KS_test=kstest(P_radio,P3)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio only pulsars= {test_stat}")
print(f"p_value of P for the radio only pulsars={p_value}")

#KS test 1D python (radio population, P<100ms)
P_rad_trunc_sim,Pdot_rad_trunc_sim,P_rad_trunc_obs,Pdot_rad_trunc_obs=[],[],[],[]
for i in range(len(P_radio)):
    if P_radio[i]<0.1:
        P_rad_trunc_sim.append(P_radio[i])
        Pdot_rad_trunc_sim.append(P_dot_radio[i])

for i in range(len(P3)):
    if P3[i]<0.1:
        P_rad_trunc_obs.append(P3[i])
        Pdot_rad_trunc_obs.append(P_dot3[i])

KS_test=kstest(Pdot_rad_trunc_sim,Pdot_rad_trunc_obs)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO ONLY PULSARS WITH P<100 ms----\n")
print(f"d_value of Pdot KS test for the radio only pulsars= {test_stat}")
print(f"p_value of Pdot for the radio only pulsars={p_value}")

KS_test=kstest(P_rad_trunc_sim,P_rad_trunc_obs)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio only pulsars= {test_stat}")
print(f"p_value of P for the radio only pulsars={p_value}")

#print(len(P_rad_trunc_sim))
#print(len(P_rad_trunc_obs))

#histSIM,xsim_edges,ysim_edges=np.histogram2d(log_P,log_Pdot,bins=(50,50))
#histobs,xobs_edges,yobs_edges=np.histogram2d(log_P2,log_Pdot2,bins=(50,50))
#histSIM=np.rot90(histSIM)
#histSIM=np.rot90(histSIM)
#histobs=np.rot90(histobs)
#histobs=np.rot90(histobs)

#hist_diff=((histSIM/len(log_P))-(histobs/len(log_P2)))/(((histSIM/len(log_P))+(histobs/len(log_P2)))**1.0)

#Make the plots

condition = [Pdot2 < Pdot3 for Pdot2, Pdot3 in zip(Pdot_line2,Pdot_line3)]

#P-Pdot observed MSP
plt.figure(0)
#plt.scatter(P,P_dot,c='red',marker='o',s=5,label='Simulation data',zorder=2)
plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF & 3PC data') #whole pop
plt.scatter(data_3PC_filtered['P0'],data_3PC_filtered['P1'],c='blue',marker='o',s=5) #whole pop
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[15]*1.0,P_dot_Edot1e28W[15]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(fontsize='small')
plt.savefig('P_Pdot_plot_observation_only.pdf',dpi=300)
plt.close()

#P-Pdot plot all pulsars
plt.figure(1)
#plt.scatter(period_young,pdot_young,c='red',marker='o',s=5,label='Young < 1e9 yr',zorder=2)
#plt.scatter(period_old,pdot_old,c='blue',marker='o',s=5,label='Old > 1e9 yr',zorder=2)
plt.scatter(P_spinup2,P_dot_spinup2,c='red',marker='o',s=5,label='Simulation',zorder=2)
plt.scatter(P_cano2,P_dot_cano2,c='purple',marker='o',s=5,label='ATNF data - normal pulsars') #cano pop
plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data - recycled pulsars') #recycled pop
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=35,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1000]*1.0,P_dot_Edot1e28W[1000]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=35,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=35,zorder=3)
plt.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.5,P_dot_Edot1e31W[1]*0.25,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=35,zorder=3)
#Plot the B lines
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-15,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[7000]*0.6,P_dot_B1e5[7000]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-15,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[7000]*0.6,P_dot_B1e6[7000]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-15,zorder=3)
plt.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[7000]*0.6,P_dot_B1e7[7000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-15,zorder=3)
plt.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[7000]*0.6,P_dot_B1e8[7000]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-15,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(fontsize='x-small')
plt.savefig('P_Pdot_plot.pdf',dpi=300)
plt.close()

#P-Pdot plot radio pulsars only 
plt.figure(2)
plt.scatter(P_radio,P_dot_radio,c='red',marker='o',s=10,label='Simulation',zorder=2)
plt.scatter(P3,P_dot3,c='blue',marker='o',s=10,label='ATNF data',zorder=1)
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )
#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[15]*1.0,P_dot_Edot1e28W[15]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.xlim(1e-3,1)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_r.pdf')
plt.close()

#P-Pdot plot gamma pulsars only
plt.figure(3)
plt.scatter(P_gamma,P_dot_gamma,c='green',marker='o',s=10,label='Simulation',zorder=3)
#plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )
#plt.scatter(P4,P_dot4,c='yellow',marker='o',s=10,label='ATNF data')
plt.scatter(data_3PC_filtered_complete['P0'],data_3PC_filtered_complete['P1'],c='blue',marker='o',s=10,label='3PC data',zorder=2)
#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[15]*1.0,P_dot_Edot1e28W[15]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.xlim(1e-3,1)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_g.pdf')
plt.close()

#P-Pdot plot radio gamma pulsars
plt.figure(4)
plt.scatter(P5,P_dot5,c='blue',marker='o',s=10,label='ATNF data')
plt.scatter(P_radio_gamma,P_dot_radio_gamma,c='purple',marker='o',s=10,label='Simulation')
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )
#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[15]*1.0,P_dot_Edot1e28W[15]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.xlim(1e-3,1)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_rg.pdf')
plt.close()

#Positions plot
plt.figure(5)
plt.scatter(x_spinup,y_spinup,s=2,c='red',alpha=0.4,label='Simulation data',zorder=2)
#plt.scatter(x_radio,y_radio,s=2,c='green',alpha=0.4,label='r-Simulation data',zorder=2)
#plt.scatter(x_radio_gamma,y_radio_gamma,s=2,c='purple',alpha=0.4,label='rg-Simulation data',zorder=2)
#plt.scatter(y_dirson22,x_dirson22,s=2,c='green',alpha=0.7,label='Dirson et al.(2022)')
plt.scatter(x2,y2,s=2,c='blue',alpha=0.5,label='ATNF data',zorder=1) #Whole pop
#plt.scatter(x_init,y_init,s=2,c='blue',alpha=0.6,label='Initial simulation data position',zorder=3)
plt.scatter([0],[8.5],c='yellow',marker='o',s=10,label='The Sun',zorder=3) #position of the sun
plt.xlim(-20,20)
plt.ylim(-20,20)
#plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend()
plt.savefig('Positions_detected_pulsars.pdf',dpi=300)
plt.close()

#Distance histogram
plt.figure(6)
plt.hist(distance_spinup,bins=20,range=(0,25),edgecolor='red',color='red',alpha=0.5,label='Simulation',zorder=3,histtype='step')
plt.hist(d2,bins=20,range=(0,25),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step') #Whole pop
#plt.hist(dist_dirson22,bins=20,range=(0,25),edgecolor='black',color='green',alpha=0.8,label='Dirson et al. (2022)')
plt.legend()
plt.yscale('log')
plt.xlabel('d (kpc)')
plt.ylabel('Frequency')
#plt.title('Distance to earth of the pulsars')
plt.savefig('histo_dist.pdf',dpi=300)
plt.close()

#Age histogram
plt.figure(7)
plt.hist(log_charac_age,bins=20,range=(6,11),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step')
plt.hist(log_age2,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step') #Whole pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
plt.yscale('log')
#plt.title('Histogram of the age of the detected pulsars')
plt.savefig('histo_charac_age.pdf',dpi=300)
plt.close()

#Latitude histogram
plt.figure(8)
plt.hist(latitude_spinup,bins=20,range=(-60,60),edgecolor='red',color='red',alpha=0.5,label='Simulation',zorder=3,histtype='step')
#plt.hist(lat_dirson22,bins=20,range=(-60,60),edgecolor='black',color='green',alpha=0.5,label='Dirson et al.(2022)',zorder=1)
plt.hist(latitude_ATNF,bins=20,range=(-60,60),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',zorder=2,histtype='step') #Whole pop
plt.yscale('log')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the latitude of the detected pulsars')
plt.savefig('histo_latitude.pdf',dpi=300)
plt.close()

#Log(P) histogram
plt.figure(9)
plt.hist(log_P,bins=20,range=(-3,0),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step')
plt.hist(log_Pa,bins=20,range=(-4,1.5),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step') #Whole pop 
#plt.hist(P_old,bins=20,range=(-2,1.5),edgecolor='black',color='green',alpha=0.5,label='Old pulsars')
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period of the detected pulsars')
plt.savefig('histo_period.pdf',dpi=300)
plt.close()

#Log(P_dot) histogram
plt.figure(10)
plt.hist(log_Pdot,bins=20,range=(-23,-15),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step')
plt.hist(log_Pdot2,bins=20,range=(-24,-16),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step') #Whole pop
plt.legend()
plt.xlabel(r'Log ($\dot{P}$) ($\ \dot{P}$ in $s.s^{-1})$')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period derivative of the detected pulsars')
plt.savefig('histo_pdot.pdf',dpi=300)
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(11)
plt.hist(cos_alpha_spinup,bins=20,range=(-1,1),edgecolor='red',color='red',alpha=0.5,density=True,label=r'cos($\alpha$)',histtype='step')
#plt.hist(cos_alpha0,bins=20,range=(-1,1),edgecolor='blue',color='blue',alpha=0.5,label='cos(alpha0)',histtype='step')
plt.legend()
plt.xlabel(r'cos($\alpha$)')
plt.ylabel('p.d.f')
#plt.title('Histogram of the inclination angle of the detected pulsars')
plt.savefig('histo_cosalpha.pdf')
plt.close()

#Velocities
plt.figure(14)
plt.hist(v_norm,bins=20,range=(0,500),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Simulation')
plt.legend()
plt.xlabel('Velocities of the pulsars in km/s')
plt.ylabel('Frequency')
#plt.title("Histogram of the velocities of the detected pulsars in the simulation")
plt.savefig('histo_velocities.pdf')
plt.close()

#Spin-velocity angles histogram
plt.figure(16)
plt.hist(PA_spinup,bins=20,range=(0,180),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Frequency')
plt.yscale('log')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle.pdf',dpi=300)
plt.close()

#Spin-velocity angle (old pulsars) histogram
plt.figure(17)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees for pulsars older than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_old.pdf')
plt.close()

#Spin-velocity angle (young pulsars) histogram
plt.figure(18)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Simulation')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees for pulsars younger than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_young.pdf')
plt.close()

#Spin-velocity angle (young pulsars and old) histogram
plt.figure(19)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=0.5,histtype='step',label='Age < 10Myr',zorder=2)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Age > 10Myr')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_young_and_old.pdf',dpi=300)
plt.close()

#Plot age=f(spin-vel angle)
plt.figure(20)
plt.scatter(PA_spinup,log_age,c='red',marker='o',s=5,label='Simulation data')
plt.xlim(0,180)
plt.ylim(0,12)
#plt.yscale('log')
#plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Log(age) in yr')
plt.legend()
plt.savefig('spin_vel_age_plot.pdf')
plt.close()

#Initial velocities
plt.figure(21)
plt.hist(v0_norm,bins=20,range=(0,500),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Simulation')
plt.legend()
plt.xlabel('Birth kick velocities of the pulsars in km/s')
plt.ylabel('Frequency')
#plt.title("Histogram of the velocities of the detected pulsars in the simulation")
plt.savefig('histo_BK_velocities.pdf')
plt.close()

#Duration of accretion
t_acc_no_zero=[x for x in t_acc_spinup if x!=0]
for i in range(len(t_acc_no_zero)):
    t_acc_no_zero[i]=np.log10(t_acc_no_zero[i]/(365*24*3600))

t_acc_no_zero2=[x for x in t_acc_spinup2 if x!=0]
for i in range(len(t_acc_no_zero2)):
    t_acc_no_zero2[i]=np.log10(t_acc_no_zero2[i]/(365*24*3600))

plt.figure(22)
plt.hist(t_acc_no_zero,bins=20,range=(4.8,10),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel('Duration of accretion in yr (log scale)')
plt.ylabel('Frequency')
plt.savefig('histo_duration_accretion.pdf',dpi=300)
plt.close()

#Initial mass of NS
plt.figure(23)
plt.hist(M_NS_i_spinup,bins=20,range=(0,3),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Initial mass of the NS in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_initmass_NS.pdf',dpi=300)
plt.close()

#Final mass of NS
plt.figure(24)
plt.hist(M_NS_f_spinup,bins=20,range=(0,3),edgecolor='red',color='red',alpha=1,histtype='step',label='Simulation')
plt.axvline(x=2.28, color='black', linestyle='--', linewidth=2, label='Theoretical maximum mass\nreachable by a NS from\nRuiz et al. (2018)')
#plt.hist(mpulsar_obs2,bins=20,range=(0,3),edgecolor='blue',color='blue',density=True,alpha=1,histtype='step',label='Measurements obtained in\nbinary pulsars system')
plt.xlabel(r'Final mass of the NS in $M_{\odot}$')
plt.ylabel('Number of pulsars')
plt.legend(fontsize='small')
plt.savefig('histo_finalmass_NS.pdf',dpi=300)
plt.close()

#Accreted mass
mean_MNS_accr=np.mean(M_NS_accr)
md_MNS_accr=np.median(M_NS_accr)
print(f'Mean accreted mass = {mean_MNS_accr} M_odot')
print(f'Median accreted mass = {md_MNS_accr} M_odot')
plt.figure(25)
plt.hist(M_NS_accr,bins=20,range=(0,1.2),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Accreted mass by the NS in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_accretedmass_NS.pdf',dpi=300)
plt.close()

#Mdot mean of NS
plt.figure(26)
plt.hist(np.log10(Mdot_mean),bins=20,range=(-11.5,-7),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'$\dot{M}_{\rm mean}$ of NS')
plt.ylabel('Frequency')
plt.savefig('histo_Mdotmean_NS.pdf',dpi=300)
plt.close()

#Semi major axis
plt.figure(27)
plt.hist(a_bin_spinup,bins=20,range=(0,100),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Semi major axis in $R_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_abin.pdf',dpi=300)
plt.close()

#Eccentricity
plt.figure(28)
plt.hist(ecc_bin_spinup,bins=20,range=(0,0.01),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Eccentricity')
plt.ylabel('Frequency')
plt.savefig('histo_ecc.pdf',dpi=300)
plt.close()

#Stellar type of the companion 
count_stellar_type = collections.Counter(stellar_phase_str)
categories = list(count_stellar_type.keys())
frequencies = list(count_stellar_type.values())
plt.figure(29, figsize=(6, 6))
x_positions = range(len(categories))
plt.bar(x_positions, frequencies, width=0.6, align='center')
plt.xticks(x_positions, categories, rotation=45, ha='right')
plt.xlabel("Stellar type of the companion")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig("histo_stellar_type.pdf", dpi=300)
plt.close()

#Remnant type of the companion
count_remnant_type = collections.Counter(remnant_type_str)
categories = list(count_remnant_type.keys())
frequencies = list(count_remnant_type.values())
plt.figure(30)
#plt.figure(30, figsize=(6, 6))
x_positions = range(len(categories))
plt.bar(x_positions, frequencies, width=0.6, align='center')
plt.xticks(x_positions, categories, rotation=45, ha='right')  
plt.xlabel("Type of the companion")
plt.ylabel("Number of pulsars")
plt.yscale('log')
plt.tight_layout()  
plt.savefig("histo_remnant_type.pdf", dpi=300)
plt.close()

#Final mass of the companion
plt.figure(32)
plt.hist(Mc_f_spinup,bins=20,range=(0,4),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Final mass of the companion in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_finalcompmass.pdf',dpi=300)
plt.close()

#Initial mass of the companion
plt.figure(33)
plt.hist(Mc_i_spinup,bins=20,range=(0,25),edgecolor='red',color='red',alpha=1,histtype='step')
plt.xlabel(r'Initial mass of the companion in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_initcompmass.pdf',dpi=300)
plt.close()

three_lowest_comp_mass=sorted(Mc_i_spinup, reverse=False)[:50]

print(f'The lowest companion mass : {three_lowest_comp_mass}')

P_orb_hr=[]
for i in range(len(P_orb)):
    P_orb_hr.append(P_orb[i]*24.0)

#P_orb=f(M_2) spiders
plt.figure(34)
category_to_group = {
    "Main Sequence": "MS",
    "Terminal Main Sequence": "MS",
    "Shell H Burning": "Giant",
    "Core He Burning": "Giant",
    "Terminal Core HeBurning": "Giant",
    "Shell He Burning": "AGB",
    "HeWD": "HeWD",
    "COWD": "COWD",
    "ONeWD": "ONeWD",
    "NS ECSN": "NS",
    "NS CCSN": "NS",
    "BH": "BH",
    "Isolated": "Isolated"
}
markers={'HeWD':'s','COWD':'x','ONeWD':'D','NS':'d','BH':'1','MS':'^','Giant':'2','AGB':'3','Isolated':'*'}
colors = {
    'HeWD': 'green',
    'COWD': 'purple',
    'ONeWD': 'orange',
    'NS': 'black',
    'BH': 'cyan',
    'MS': 'red',
    'AGB': 'brown',
    'Giant': 'pink',
    'Isolated': 'darkgoldenrod'
}
x_min,x_max=3e-3,5e0
y_min,y_max=1e-2,10
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or 
             (companion_group == 'Isolated' and fallback_group == group)) and 
            not np.isnan(p) and
            (x_min <= m <= x_max) and
            (y_min <= p <= y_max)
            for companion_group, fallback_group, m , p 
            in zip(grouped_remnant_type, stellar_phase_str,Mc_f_spinup , P_orb)]
    if any(mask):
        plt.scatter(
            [m for m, flag in zip(Mc_f_spinup, mask) if flag],  # Filtered mass
            [p for p, flag in zip(P_orb, mask) if flag],        # Orbital period filtered
            label=f'{group}', marker=marker, s=5, color=colors[group]
        )

grouped_remnant_type2 = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str2, stellar_phase_str2)
]

plt.axvline(x=0.1, linestyle='--',color='black',label='Black widow/Redback boundary')
plt.text(1.5e0, 2e-2, "Redback side",color='red')
plt.text(4e-2, 2e-2, "Black widow side", color='black')
plt.yscale('log')
plt.xscale('log')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.ylabel(r'$P_{orb}$ (days)')
plt.xlabel(r'Companion mass ($M_{\odot}$)')
plt.legend(fontsize='small')
plt.savefig('Porb_M2_plot.pdf',dpi=300)
plt.close()

#P-Pdot with legend of companions
plt.figure(35)
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group))
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [p for p, flag in zip(P_spinup, mask) if flag],  # Filtered mass
        [pdot for pdot, flag in zip(P_dot_spinup, mask) if flag],        # Orbital period filtered
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[2]*1.0,P_dot_Edot1e23W[2]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e28W[0]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e25W[0]*1.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[700]*0.6,P_dot_B5e6[640]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(1e-24,5e-14)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(fontsize='x-small')
plt.savefig('P_Pdot_plot_sim_with_companion.pdf',dpi=300)
plt.close()

#P=f(P_orb)
for comp, p in zip(grouped_remnant_type, P_orb):
    if comp == "NS":
        print("NS orbital period =", p)
plt.figure(36)
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_orb)]
    plt.scatter(
        [p for p, flag in zip(P_orb, mask) if flag],  # Filtered mass
        [pbis for pbis, flag in zip(P_spinup, mask) if flag],        # Orbital period filtered
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-2,4000)
plt.ylim(1e-3,5e0)
plt.xlabel('Orbital period of the binary in days')
plt.ylabel(r'Spin period $P$ (s)')
plt.legend(fontsize='x-small')
plt.savefig('P_Porb_plot_sim.pdf',dpi=300)
plt.close()

#P-Pdot observed MSP
companion_styles = {
    'UL_group': {'color': 'blue', 'marker': 'o', 'label': 'Ultra light'},
    'He_group': {'color': 'green', 'marker': 's', 'label': 'He WD'},
    'MS': {'color': 'red', 'marker': '^', 'label': 'MS star'},
    'NS': {'color': 'black', 'marker': 'd', 'label': 'NS'},
    'CO': {'color': 'purple', 'marker': 'x', 'label': 'CO WD'},
    'NULL': {'color': 'yellow', 'marker': '*', 'label': 'Isolated'},
    'UNKN': {'color': 'grey', 'marker': 'h', 'label': 'Unknown companion'}
}
processed_companions = [
    'UL_group' if c in ['UL', 'ULT'] else
    'He_group' if c in ['He', 'HeT'] else
    c
    for c in type_comp2
]
plt.figure(37)
for P, P_dot, compa in zip(P2, P_dot2, processed_companions):
    style = companion_styles[compa]
    plt.scatter(P, P_dot, color=style['color'], marker=style['marker'], label=style['label'])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(),fontsize='x-small')
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[2]*1.0,P_dot_Edot1e23W[2]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e28W[0]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e25W[0]*1.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[700]*0.6,P_dot_B5e6[640]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(1e-24,1e-13)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_observation_only_with_comptype.pdf',dpi=300)
plt.close()

#e=f(P_orb) obs
plt.figure(38)
P_orb_forplot,Ecc_for_plot,compa_for_plot=[],[],[]
for porb, ecc, compa in zip(P_orb2, Ecc2, processed_companions):
    if porb!=-1 and ecc!=-1:
        P_orb_forplot.append(porb)
        Ecc_for_plot.append(ecc)
        compa_for_plot.append(compa)
for e, porb, compa in zip(Ecc_for_plot, P_orb_forplot, compa_for_plot):
    style = companion_styles[compa]
    plt.scatter(porb, e, color=style['color'], marker=style['marker'], label=style['label'])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
plt.xscale('log')
plt.yscale('log')
plt.xlim(4e-2,4000)
plt.ylim(1e-8,1)
plt.xlabel('Orbital period in days')
plt.ylabel('Eccentricity')
plt.savefig('porb_e_obs.pdf',dpi=300)
plt.close()

#e=f(P_orb) sim
plt.figure(39)
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_orb)]
    plt.scatter(
        [p for p, flag in zip(P_orb, mask) if flag],  # Filtered mass
        [e for e, flag in zip(ecc_bin_spinup, mask) if flag],        # Orbital period filtered
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
plt.legend(fontsize='small')
#plt.xlim(1e-2,4000)
#plt.ylim(0,0.01)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Orbital period in days')
plt.ylabel('Eccentricity')
plt.savefig('porb_e_sim.pdf',dpi=300)
plt.close()

#P=f(P_orb) obs
plt.figure(40)
Porb_for_plot2,P2_for_plot,compa_for_plot2=[],[],[]
for porb, p, compa in zip(P_orb2, P2, processed_companions):
    if porb!=-1:
        Porb_for_plot2.append(porb)
        P2_for_plot.append(p)
        compa_for_plot2.append(compa)
for p, porb, compa in zip(P2_for_plot, Porb_for_plot2, compa_for_plot2):
    style = companion_styles[compa]
    plt.scatter(porb, p, color=style['color'], marker=style['marker'], label=style['label'])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(),fontsize='x-small')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-2,4000)
plt.ylim(1e-3,5e0)
plt.xlabel('Orbital period of the binary in days')
plt.ylabel(r'Spin period $P$ (s)')
plt.savefig('P_Porb_plot_obs.pdf',dpi=300)
plt.close()

#M_WD comparison between real mass and relationship of Tauris & Savonije (1999) 2d <= P_orb <= 100d
plt.figure(41)
plt.hist(M_WD_sim_g2l100,bins=20,range=(0.17,0.49),edgecolor='red',color='red',alpha=0.5,label='Real masses from the simulation',histtype='step')
plt.hist(M_WD_rel1,bins=20,range=(0.17,0.49),edgecolor='blue',color='blue',alpha=0.5,label=r'Mass computed with Tauris & Savonije (1999) relationship',histtype='step') #Whole pop
plt.legend()
plt.xlabel(r'Mass of the white dwarf in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_MWD_porb_l100.pdf',dpi=300)
plt.close()

#M_WD comparison between real mass and relationship of Tauris & Savonije (1999) P_orb > 100d
plt.figure(42)
plt.hist(M_WD_sim_g100,bins=20,range=(0.17,0.49),edgecolor='red',color='red',alpha=0.5,label='Real masses from the simulation',histtype='step')
plt.hist(M_WD_rel2,bins=20,range=(0.17,0.49),edgecolor='blue',color='blue',alpha=0.5,label=r'Mass computed with Tauris & Savonije (1999) relationship',histtype='step') #Whole pop
plt.legend()
plt.xlabel(r'Mass of the white dwarf in $M_{\odot}$')
plt.ylabel('Frequency')
plt.savefig('histo_MWD_porb_g100.pdf',dpi=300)
plt.close()

#M_WD comparison between real mass and relationship of Tauris & Savonije (1999) P_orb > 100d relative difference histogram
plt.figure(43)
plt.hist(ecc_rel2,bins=20,range=(0,1),edgecolor='red',color='red',alpha=0.5,histtype='step')
plt.legend()
plt.xlabel(r'Relative difference')
plt.ylabel('Frequency')
plt.savefig('histo_reldiff_MWD_porb_g100.pdf',dpi=300)
plt.close()

#M_WD comparison between real mass and relationship of Tauris & Savonije (1999) 2 <= P_orb <= 100d relative difference histogram
plt.figure(44)
plt.hist(ecc_rel1,bins=20,range=(0,1),edgecolor='red',color='red',alpha=0.5,histtype='step')
plt.legend()
plt.xlabel(r'Relative difference')
plt.ylabel('Frequency')
plt.savefig('histo_reldiff_MWD_porb_l100.pdf',dpi=300)
plt.close()

#B histogram
plt.figure(45)
plt.hist(np.log10(Bf_spinup),bins=20,range=(3,8),edgecolor='red',color='red',alpha=0.5,label='Simulation')
plt.legend()
#plt.xscale('log')
plt.xlabel('Magnetic field (T)')
plt.ylabel('Frequency')
plt.savefig('histo_Bsim.pdf')
plt.close()

#B=f(age)
plt.figure(46)
plt.scatter(log_age,np.log10(Bf_spinup),label='Simulation')
plt.legend()
plt.ylabel('Magnetic field in logscale (B in T)')
plt.xlabel('Age in logscale (age in yr)')
plt.savefig('plot_age_Bsim.pdf')
plt.close()

#P-Pdot possible transients 
plt.figure(47)
plt.scatter(P_trans,Pdot_trans,c='red',marker='o',s=5,label='Possible transients MSPs',zorder=2)
#plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data') #whole pop
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[2]*1.0,P_dot_Edot1e23W[2]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e28W[0]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e25W[0]*1.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[700]*0.6,P_dot_B5e6[640]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(9e-25,5e-14)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(fontsize='small')
plt.savefig('P_Pdot_plot_transients.pdf',dpi=300)

#NS ellipticities
plt.figure(48)
plt.hist(np.log10(NS_ell_spinup),bins=20,range=(-12,-5),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Simulation')
plt.legend()
plt.xlabel(r'Ellipticity of NS in logscale')
plt.ylabel('Frequency')
plt.savefig('histo_ellipticityNS.pdf')
plt.close()

#fit in order to find distribution of M_WD
#M_WD_only = np.array(M_WD_only)
#counts, bin_edges = np.histogram(M_WD_only, bins=20, density=True)
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

def exponential(x, A, lambda_):
    return A * np.exp(-lambda_ * x)

def lognormal(x, A, mu, sigma):
    return (A / x) * np.exp(-((np.log(x) - mu) ** 2) / (2 * sigma ** 2))

#popt, pcov = curve_fit(lognormal, bin_centers, counts, p0=[max(counts),np.log(np.mean(M_WD_only)),np.std(np.log(M_WD_only))])
#popt2, pcov2 = curve_fit(exponential, bin_centers, counts, p0=[max(counts), 1/np.mean(M_WD_only)])
#x_fit = np.linspace(min(M_WD_only), max(M_WD_only), 1000)

#M_WD companion
#plt.figure(49)
#plt.hist(M_WD_only,bins=20,range=(0,2),edgecolor='red',color='red',density=True,alpha=0.5,histtype='step')
#plt.plot(x_fit, lognormal(x_fit, *popt), 'b-', label=f"Fit LogNorm: A={popt[0]:.2f}, μ={popt[1]:.2f}, σ={popt[2]:.2f}")
#plt.plot(x_fit, exponential(x_fit, *popt2), 'g-', label=f"Fit Exp: A={popt2[0]:.2f}, λ={popt2[1]:.2f}")
#print(f"Paramètres estimés : A = {popt[0]:.3f}, λ = {popt[1]:.3f}")
#plt.legend()
#plt.xlabel(r'M$_{WD}$ companion in $M_{\odot}$')
#plt.ylabel('Frequency')
#plt.savefig('histo_MWDonly.pdf',dpi=300)
#plt.close()

#Real age detected MSPs
plt.figure(50)
plt.hist(log_age,bins=20,range=(8,10.5),edgecolor='red',color='red',alpha=0.5,histtype='step',label='Simulation')
plt.legend()
plt.xlabel(r'Log(age) (age in yr)')
plt.ylabel('Frequency')
plt.savefig('histo_real_age.pdf',dpi=300)
plt.close()

#age_real=f(age_sp)
plt.figure(51)
age_norm=[]
for i in range(len(age_spinup)):
    age_norm.append((age_spinup[i]/(365*24*3600))/(10**(log_charac_age[i])))
log_charac_age=np.array(log_charac_age)
scatter=plt.scatter(10**(log_charac_age),age_norm,marker='o',s=5,label='Simulation')
plt.axvline(x=0.05e6,color='red',linestyle='--',label=r'$\tau_d=0.05$ Myr')
plt.axvline(x=0.6e6,color='blue',linestyle='--',label=r'$\tau_d=0.6$ Myr')
plt.axvline(x=1.3e6,color='green',linestyle='--',label=r'$\tau_d=1.3$ Myr')
plt.xlabel('Characteristic age (yr)')
plt.ylabel('Real age/characteristic age')
plt.yscale('log')
plt.xscale('log')
#plt.xlim(5e8,11e9)
#plt.ylim(4e-2,1.2)
plt.legend()
plt.grid(alpha=0.5, linestyle='-')
plt.savefig(f'plot_age_spindownage.pdf',dpi=300)
plt.close()

#P-M2 plot
plt.figure(52)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_orb)]
    plt.scatter(
        [m for m, flag in zip(Mc_f_spinup, mask) if flag],     # Filtered mass
        [p for p, flag in zip(P_spinup, mask) if flag],        # Orbital period filtered
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
plt.yscale('log')
plt.xscale('log')
#plt.xlim(1e-2,1)
#plt.ylim(1e-2,1000)
plt.ylabel('P (s)')
plt.xlabel(r'Mass of companion in $M_{\odot}$')
plt.legend(fontsize='small')
plt.savefig('Ppulsar_M2_plot.pdf',dpi=300)
plt.close()

#print(taud_spinup)
#P-Pdot plot all pulsars
plt.figure(53)
#plt.scatter(period_young,pdot_young,c='red',marker='o',s=5,label='Young < 1e9 yr',zorder=2)
#plt.scatter(period_old,pdot_old,c='blue',marker='o',s=5,label='Old > 1e9 yr',zorder=2)
sc=plt.scatter(P_spinup,P_dot_spinup,c=taud_spinup,cmap='plasma',marker='o',s=5,label='Simulation data',zorder=2)
cbar=plt.colorbar(sc)
cbar.set_label(r'$\tau_d$ (Myr)', fontsize=8)
#plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data') #whole pop
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[2]*1.0,P_dot_Edot1e23W[2]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e28W[0]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[0]*1.0,P_dot_Edot1e25W[0]*1.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[700]*0.6,P_dot_B5e6[640]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)

plt.xlim(1e-3,1e0)
plt.ylim(1e-23,1e-16)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(fontsize='x-small')
plt.savefig('P_Pdot_taud_plot.pdf',dpi=300)
plt.close()

#Final mass of the companion
plt.figure(54)
plt.hist(Mc_bigP,bins=20,range=(0,150),edgecolor='red',color='red',density=True,alpha=1,histtype='step',label=r'$P\geq0.02$ s')
plt.hist(Mc_lowP,bins=20,range=(0,150),edgecolor='blue',color='blue',density=True,alpha=1,histtype='step',label=r'$P<0.02$ s')
plt.xlabel(r'Initial mass of the companion in $M_{\odot}$')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('histo_initcompmass_P01.pdf',dpi=300)
plt.close()

for i in range(len(wr_spinup)):
    if wr_spinup[i]>360:
        wr_spinup[i]=wr_spinup[i]-int(wr_spinup[i]/360.0)*360

#For w_r with ISM+instrument effects
a_wr,b_wr,r_value,p_value,std_err=linregress(np.log10(P_r_or_rg),np.log10(wr_r_or_rg))
reglinx = np.logspace(np.log10(min(P_r_or_rg)), np.log10(max(P_r_or_rg)), num=len(P_r_or_rg))
regliny = 10**(a_wr * np.log10(reglinx) + b_wr)
residuals_y= np.log10(wr_r_or_rg) - np.log10(regliny)
sigma_y = np.sqrt(np.sum(residuals_y**2) / (len(P_r_or_rg) - 2))
Sxx=np.sum((np.log10(P_r_or_rg) - np.mean(np.log10(P_r_or_rg)))**2)
std_err_b=sigma_y*np.sqrt(1.0/len(P_r_or_rg)+((np.mean(np.log10(P_r_or_rg)))**2/Sxx))
a_wrplus=a_wr+std_err
b_wrplus=b_wr+std_err_b
a_wrminus=a_wr-std_err
b_wrminus=b_wr-std_err_b
reglinyplus = 10**(a_wr * np.log10(reglinx) + b_wrplus)
reglinyminus = 10**(a_wr * np.log10(reglinx) + b_wrminus)

#for w_r without ISM+instrument effects
indices=[i for i, val in enumerate(P_r_or_rg) if val<1]
P_r_or_rg2=[P_r_or_rg[i] for i in indices]
w_geometry_r_or_rg2=[w_geometry_r_or_rg[i] for i in indices]
a_wr2,b_wr2,r_value2,p_value2,std_err2=linregress(np.log10(P_r_or_rg2),np.log10(w_geometry_r_or_rg2))
reglinx2 = np.logspace(np.log10(min(P_r_or_rg2)), np.log10(max(P_r_or_rg2)), num=len(P_r_or_rg2))
regliny2 = 10**(a_wr2 * np.log10(reglinx2) + b_wr2)
residuals_y2= np.log10(w_geometry_r_or_rg2) - np.log10(regliny2)
sigma_y2 = np.sqrt(np.sum(residuals_y2**2) / (len(P_r_or_rg2) - 2))
Sxx2=np.sum((np.log10(P_r_or_rg2) - np.mean(np.log10(P_r_or_rg2)))**2)
std_err_b2=sigma_y2*np.sqrt(1.0/len(P_r_or_rg2)+((np.mean(np.log10(P_r_or_rg2)))**2/Sxx2))
a_wrplus2=a_wr2+std_err2
b_wrplus2=b_wr2+std_err_b2
a_wrminus2=a_wr2-std_err2
b_wrminus2=b_wr2-std_err_b2
reglinyplus2 = 10**(a_wr2 * np.log10(reglinx2) + b_wrplus2)
reglinyminus2 = 10**(a_wr2 * np.log10(reglinx2) + b_wrminus2)

#W10 Obs
P_obs_combined=P3+P5
w10_combined=w10_r+w10_rg
indices = [i for i, val in enumerate(w10_combined) if val != 0]
w10_combined_filtered = [w10_combined[i] for i in indices]
P_obs_combined_filtered = [P_obs_combined[i] for i in indices]
a_wr3,b_wr3,r_value3,p_value3,std_err3=linregress(np.log10(P_obs_combined_filtered),np.log10(w10_combined_filtered))
reglinx3 = np.logspace(np.log10(min(P_obs_combined_filtered)), np.log10(max(P_obs_combined_filtered)), num=len(P_obs_combined_filtered))
regliny3 = 10**(a_wr3 * np.log10(reglinx3) + b_wr3)
residuals_y3= np.log10(w10_combined_filtered) - np.log10(regliny3)
sigma_y3 = np.sqrt(np.sum(residuals_y3**2) / (len(P_obs_combined_filtered) - 2))
Sxx3=np.sum((np.log10(P_obs_combined_filtered) - np.mean(np.log10(P_obs_combined_filtered)))**2)
std_err_b3=sigma_y3*np.sqrt(1.0/len(P_obs_combined_filtered)+((np.mean(np.log10(P_obs_combined_filtered)))**2/Sxx3))
a_wrplus3=a_wr3+std_err3
b_wrplus3=b_wr3+std_err_b3
a_wrminus3=a_wr3-std_err3
b_wrminus3=b_wr3-std_err_b3
reglinyplus3 = 10**(a_wr3 * np.log10(reglinx3) + b_wrplus3)
reglinyminus3 = 10**(a_wr3 * np.log10(reglinx3) + b_wrminus3)

print(len(w10_combined_filtered))
print(len(w_geometry_r_or_rg))

#f(log(wr))=logP
plt.figure(55)
#plt.scatter(P_r_or_rg,wr_r_or_rg,c='green',marker='o',s=5,label='Simulation with ISM and instrumental effect')
plt.scatter(P_r_or_rg2,w_geometry_r_or_rg2,c='red',marker='o',s=5,label='Simulation')
#plt.plot(reglinx,regliny,linestyle='-',label=r'$\log$($w_r$) = %.2f $\log(P)$ + %.2f' % (a_wr, b_wr),c='green')
plt.plot(reglinx2,regliny2,linestyle='-',label=r'$\log$($w_r$) = (%.2f$\pm$%.2f) $\log(P)$ + (%.2f$\pm$%.2f) ' % (a_wr2,std_err2,b_wr2,std_err_b2),c='red')
plt.scatter(P3,w10_r,c='blue',marker='o',s=5)
plt.scatter(P5,w10_rg,c='blue',marker='o',s=5,label='ATNF data')
plt.plot(reglinx3,regliny3,linestyle='-',label=r'$\log$($w^{ATNF}_{10}$) = (%.2f$\pm$%.2f) $\log(P)$ + (%.2f$\pm$%.2f)' % (a_wr3,std_err3,b_wr3,std_err_b3),c='blue')
plt.fill_between(reglinx2, reglinyminus2, reglinyplus2, color='red', alpha=0.1)#, label=r'$\pm 1\sigma$ interval')
plt.fill_between(reglinx3, reglinyminus3, reglinyplus3, color='blue', alpha=0.1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$w_r$ (°)')
plt.xlabel(r'$P$ (s)')
plt.legend(fontsize='small')
plt.savefig('wr_P_plot.pdf',dpi=300)
plt.close()

#Comparison gamma-ray peak separation 3PC and simulation
nb_del=0
for i in range(len(g_peak_sep_sim)):
    if g_peak_sep_sim[i-nb_del]>1.0:
        del g_peak_sep_sim[i-nb_del]
        nb_del+=1

#Examine the Pdot of highg gpeak separation
#g_peak_sep_sim_nor=[i for i in g_peak_sep_sim if i!=1e55]
#P_dot_gpeakhigh=[]
#for i in range(len(g_peak_sep_sim)):
#    if g_peak_sep_sim[i]>=0.45:
#        P_dot_gpeakhigh.append(P_dot_gamma[i])

print(f'\nNumber of pulsars gamma : {len(P_dot_gamma)} and radio-gamma: {len(P_dot_radio_gamma)} and Number of gpeap available : {len(g_peak_sep_sim)}\n')

weights_sim=np.ones_like(g_peak_sep_sim) / len(g_peak_sep_sim)
weights_obs=np.ones_like(data_3PC_filtered['PKSEP']) / len(data_3PC_filtered['PKSEP'])
plt.figure(56)
plt.hist(g_peak_sep_sim,bins=10,range=(0,0.5),edgecolor='red',color='red',alpha=1,label='Simulation',density=True,zorder=3,histtype='step')
plt.hist(data_3PC_filtered['PKSEP'],bins=10,range=(0,0.5),edgecolor='blue',color='blue',alpha=1,density=True,label='3PC data',zorder=2,histtype='step') #Recycled pop
plt.legend()
#plt.yscale('log')
plt.xlabel(r'Gamma-ray peak separation')
plt.ylabel('p.d.f')
plt.savefig('histo_gpeaksep.pdf',dpi=300)
plt.close()

#P=f(t_acc)
for i in range(len(t_acc_no_zero)):
    t_acc_no_zero[i]=10**(t_acc_no_zero[i])

for i in range(len(t_acc_no_zero2)):
    t_acc_no_zero2[i]=10**(t_acc_no_zero2[i])

plt.figure(57)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [t for t, flag in zip(t_acc_no_zero, mask) if flag],     # Duration of accretion
        [p for p, flag in zip(P_spinup, mask) if flag],        # spin period pulsar
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel('Duration of accretion (yr)')
plt.ylabel(r'$P$ (s)')
plt.savefig('P_tacc_plot.pdf',dpi=300)
plt.close()

#P=f(Mc_i)
plt.figure(58)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [m for m, flag in zip(Mc_i_spinup, mask) if flag],     # Initial mass of the companion
        [p for p, flag in zip(P_spinup, mask) if flag],        # pulsar spin period 
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'Initial mass of the companion star ($M_{\odot}$)')
plt.ylabel(r'$P$ (s)')
plt.savefig('initmc_P_plot.pdf',dpi=300)
plt.close()

#Mc_i=f(t_acc)
print(f'Number of pulsars counted : {len(P_spinup2)}\n')
plt.figure(59)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    if any(mask):
        plt.scatter(
                [t for t, flag in zip(t_acc_no_zero, mask) if flag],     # Duration of accretion
                [m for m, flag in zip(Mc_i_spinup, mask) if flag],        # Initial mass of the companion
                label=f'{group}', marker=marker, s=5, color=colors[group]
                )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.ylabel(r'Initial mass of the companion star ($M_{\odot}$)')
plt.xlabel(r'Duration of accretion (yr)')
plt.savefig('mci_tacc_plot.pdf',dpi=300)
plt.close()

#B=f(t_acc)
plt.figure(60)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [t for t, flag in zip(t_acc_no_zero, mask) if flag],     # Duration of accretion
        [b for b, flag in zip(Bf_spinup, mask) if flag],        # Magnetic field
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'Duration of accretion (yr)')
plt.ylabel(r'$B$ (T)')
plt.savefig('B_tacc_plot.pdf',dpi=300)
plt.close()

#Age=f(t_acc)
age_yr=[]
for i in range(len(age_spinup)):
    age_yr.append(age_spinup[i]/(365*24*3600))
plt.figure(61)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [t for t, flag in zip(t_acc_no_zero, mask) if flag],     # Duration of accretion
        [t2 for t2, flag in zip(age_yr, mask) if flag],          # Pulsar age
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'Duration of accretion (yr)')
plt.ylabel(r'Pulsar age (yr)')
plt.savefig('age_tacc_plot.pdf',dpi=300)
plt.close()

#a_bin=f(t_acc)
plt.figure(62)
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p)
            for companion_group, fallback_group, p
            in zip(grouped_remnant_type, stellar_phase_str, P_spinup)]
    plt.scatter(
        [t for t, flag in zip(t_acc_no_zero, mask) if flag],     # Duration of accretion
        [a for a, flag in zip(a_bin_spinup, mask) if flag],          # Semi-major axis separation of the binary
        label=f'{group}', marker=marker, s=5, color=colors[group]
    )
#plt.scatter(t_acc_no_zero,P_spinup,)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'Duration of accretion (yr)')
plt.ylabel(r'Semi-major axis of the binary ($R_{\odot}$)')
plt.savefig('a_bin_tacc_plot.pdf',dpi=300)
plt.close()

#Checking charac particular population
log_age_filtered=[]
taud_filtered=[]
for i in range(len(P_gamma)):
    if (P_gamma[i]>5e-3 or P_dot_gamma[i]>1e-19):
        log_age_filtered.append(np.log10(age_gamma[i]/(365*24*3600)))
        taud_filtered.append(taud_g[i])
print(len(taud_filtered))
plt.figure(63)
plt.scatter(log_age_filtered,taud_filtered,s=5,color='red',label='gamma-ray pulsars')
#plt.hist(log_age_filtered,bins=20,range=(8,10.5),edgecolor='red',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Age (yr in logscale)')
plt.ylabel(r'$\tau_d$ in Myr')
plt.savefig('histo_gammafiltered.pdf')
plt.close()

#Flux gamma comparison obs vs sim
for i in range(len(data_3PC_filtered['G100'])):
    try:
        data_3PC_filtered['G100'][i]=np.log10(float(data_3PC_filtered['G100'][i])*1e-3)
    except ValueError:
        data_3PC_filtered['G100'][i] = np.nan

plt.figure(64)
plt.hist(np.log10(flux_gamma_only),bins=20,range=(-16,-11),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step',density=True)
plt.hist((data_3PC_filtered['G100']),bins=20,range=(-16,-11),edgecolor='blue',color='blue',alpha=1,label='3PC data',zorder=2,histtype='step',density=True) 
plt.legend()
#plt.yscale('log')
plt.xlabel(r'Gamma-ray flux in log space (W.m$^{-2}$)')
plt.ylabel('p.d.f')
plt.savefig('histo_fgamma.png',dpi=300)
plt.close()

#Beta histogram
plt.figure(65)
plt.hist(spin_orb_angle_all,bins=36,range=(0,180),edgecolor='green',color='green',alpha=0.5,density=True,label=r'Whole simulated population',histtype='step')
plt.hist(spin_orb_angle_spinup,bins=36,range=(0,180),edgecolor='red',color='red',alpha=0.5,density=True,label=r'Simulated detected population',histtype='step')
plt.hist(spin_orb_angle_init_spinup,bins=36,range=(0,180),edgecolor='blue',color='blue',alpha=0.5,density=True,label=f'Initial angle for the detected\nsimulated population',histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel(r'$\alpha$ in °')
plt.ylabel('p.d.f')
plt.savefig('histo_spinorbit.pdf')
plt.close()

#2D histogram t_acc beta
for i in range(len(t_acc_no_zero)):
    t_acc_no_zero[i]=np.log10(t_acc_no_zero[i])

spin_orb_angle_spinup_noiso,t_acc_no_zero_noiso,spin_orb_angle_spinup_noiso_cos=[],[],[]
for i in range(len(spin_orb_angle_spinup)):
    if comp_rem_type_spinup[i]!=-1:
        spin_orb_angle_spinup_noiso.append(90 - np.abs((spin_orb_angle_spinup[i]%180) - 90))
        spin_orb_angle_spinup_noiso_cos.append(np.cos(spin_orb_angle_spinup[i]))
        t_acc_no_zero_noiso.append(t_acc_no_zero[i])

plt.figure(66)
#weights=np.ones(len(spin_orb_angle_spinup))/len(spin_orb_angle_spinup)
hist=plt.hist2d(spin_orb_angle_spinup_noiso_cos,t_acc_no_zero_noiso,bins=(40,10),range=((-1,1),(5,10)),cmap='plasma',norm=LogNorm())#,weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Number of pulsars')
#plt.yscale('log')
plt.ylabel('Duration of accretion in log scale of yr')
plt.xlabel(r'cos($\alpha$)')
plt.savefig(f'2D_map_tacc_beta.pdf',dpi=300)
plt.close()

#t_acc=f(beta)
plt.figure(67)
plt.scatter(spin_orb_angle_spinup,t_acc_no_zero,c='red',marker='o',s=5,label='Simulation')
#plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'Duration of accretion (yr)')
plt.xlabel(r'$\alpha$ in °')
plt.savefig('spin_orbit_tacc_plot.pdf',dpi=300)
plt.close()

#Zeta=f(chi) radio
plt.figure(68)
plt.scatter(alpha_r,xi_r,c='red',marker='o',s=5,label='Simulation')
plt.ylabel(r'$\zeta$ in °')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.legend(fontsize='small')
plt.savefig('zeta_chi_radio.pdf',dpi=300)
plt.close()

#Zeta=f(chi) gamma
plt.figure(69)
plt.scatter(alpha_g,xi_g,c='red',marker='o',s=5,label='Simulation')
plt.ylabel(r'$\zeta$ in °')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.legend(fontsize='small')
plt.savefig('zeta_chi_gamma.pdf',dpi=300)
plt.close()

#Zeta=f(chi) radio gamma
plt.figure(70)
plt.scatter(alpha_rg,xi_rg,c='red',marker='o',s=5,label='Simulation')
plt.ylabel(r'$\zeta$ in °')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.legend(fontsize='small')
plt.savefig('zeta_chi_rg.pdf',dpi=300)
plt.close()

#Zeta=f(chi) all
plt.figure(71)
plt.scatter(alpha_rg,xi_rg,c='blue',marker='o',s=5,label='Pulsars in both radio and gamma-ray surveys')
plt.scatter(alpha_g,xi_g,c='green',marker='x',s=5,label='Pulsars in gamma-ray surveys only')
plt.scatter(alpha_r,xi_r,c='red',marker='s',s=5,label='Pulsars in radio surveys only')
plt.ylabel(r'$\zeta$ in °')
plt.legend(fontsize='small')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.savefig('zeta_chi_all.pdf',dpi=300)
plt.close()

#Count nb points in the detectable area 
#a1,b1=70/42,20
#a2,b2=42/59,-31*42/59
#alpha_g=np.array(alpha_g)
#xi_g=np.array(xi_g)
#y1=a1*alpha_g+b1
#y2=a2*alpha_g+b2
#mask_detectable = ((xi_g > np.minimum(y1, y2)) & (xi_g < np.maximum(y1, y2)))
#n_points_detectable=np.sum(mask_detectable)
#print(f"Ratio of gamma-ray pulsars in the detectable radio area : {n_points_detectable/count_gam}")

#Longitude
plt.figure(72)
plt.hist(longitude_spinup,bins=20,range=(0,360),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
plt.hist(longitude_ATNF,bins=20,range=(0,360),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step') 
plt.legend()
plt.yscale('log')
plt.xlabel('Longitude in °')
plt.ylabel('Number of pulsars')
plt.savefig('histo_longitude.pdf',dpi=300)
plt.close()

#Galactic coordinates map
longitude_spinup_rad=np.radians(longitude_spinup)
latitude_spinup_rad=np.radians(latitude_spinup)
longitude_rad=np.radians(longitude_ATNF)
latitude_rad=np.radians(latitude_ATNF)
longitude_spinup_rad = np.remainder(longitude_spinup_rad + 2*np.pi, 2*np.pi)
longitude_spinup_rad[longitude_spinup_rad > np.pi] -= 2*np.pi
#longitude_spinup_rad = -longitude_spinup_rad
longitude_rad = np.remainder(longitude_rad + 2*np.pi, 2*np.pi)
longitude_rad[longitude_rad > np.pi] -= 2*np.pi
#longitude_rad = -longitude_rad
plt.figure(73)
plt.axes(projection="mollweide")
plt.scatter(longitude_spinup_rad,latitude_spinup_rad,s=10,c='red',alpha=1,label='Simulation data',zorder=2)
plt.scatter(longitude_rad,latitude_rad,s=10,c='blue',alpha=1,label='ATNF data',zorder=2)
plt.grid(alpha=0.3)
plt.legend(loc="upper right", bbox_to_anchor=(1.05, 1.15), borderaxespad=0.)
plt.xlabel("Galactic longitude in °")
plt.ylabel("Galactic latitude in °")
plt.savefig("plot_lat_long.pdf",dpi=300)
plt.close()

#Positions plot (all MSP born)
plt.figure(74)
plt.scatter(x_MSP,y_MSP,s=2,c='red',alpha=0.4,label='Simulation data',zorder=2)
#plt.scatter(x_radio,y_radio,s=2,c='green',alpha=0.4,label='r-Simulation data',zorder=2)
#plt.scatter(x_radio_gamma,y_radio_gamma,s=2,c='purple',alpha=0.4,label='rg-Simulation data',zorder=2)
#plt.scatter(y_dirson22,x_dirson22,s=2,c='green',alpha=0.7,label='Dirson et al.(2022)')
#plt.scatter(x2,y2,s=2,c='blue',alpha=0.5,label='ATNF data',zorder=1) #Whole pop
#plt.scatter(x_init,y_init,s=2,c='blue',alpha=0.6,label='Initial simulation data position',zorder=3)
plt.scatter([0],[8.5],c='yellow',marker='o',s=10,label='The Sun',zorder=3) #position of the sun
plt.xlim(-20,20)
plt.ylim(-20,20)
#plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend()
plt.savefig('Positions_all_MSPs_born.pdf',dpi=300)
plt.close()

#2D histogram positions
plt.figure(75)
weights=np.ones(len(x_MSP))/len(x_MSP)
hist=plt.hist2d(x_MSP,y_MSP,bins=(32,32),range=((-20,20),(-20,20)),cmap='viridis',weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Fraction of occurences of MSPs')
plt.xlabel(r'x (kpc)')
plt.ylabel(r'y (kpc)')
plt.savefig(f'2D_positions_galaxy_allMSP_born.pdf',dpi=300)
plt.close()

#Theta disk distribution
plt.figure(76)
plt.hist(theta_disc,bins=20,range=(0,360),edgecolor='red',color='red',alpha=0.5,label='Simulation',zorder=3,histtype='step')
plt.legend()
plt.xlabel('Theta in degrees')
plt.ylabel('Number of MSPs')
plt.savefig('theta_disk.pdf',dpi=300)
plt.close()

#Distance histogram
plt.figure(77)
plt.hist(dist_allMSP,bins=20,range=(0,25),edgecolor='red',color='red',alpha=0.5,label='Simulated population',zorder=3,histtype='step',density=False)
#plt.hist(d2,bins=20,range=(0,25),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step') #Whole pop
#plt.hist(dist_dirson22,bins=20,range=(0,25),edgecolor='black',color='green',alpha=0.8,label='Dirson et al. (2022)')
plt.legend()
plt.yscale('log')
plt.xlabel('R (kpc)')
plt.ylabel('Number of pulsars')
#plt.title('Distance to earth of the pulsars')
plt.savefig('histo_radial_dist.pdf',dpi=300)
plt.close()

#Ratio r_kg and R_L histogram
plt.figure(78)
plt.hist(ratio_rkG_rL_rg,bins=20,range=(0,2),edgecolor='blue',color='blue',alpha=0.5,label='Pulsars in both radio and gamma-ray surveys',zorder=3,histtype='step')
plt.hist(ratio_rkG_rL_radio,bins=20,range=(0,2),edgecolor='red',color='red',alpha=0.5,label='Pulsars in radio surveys only',zorder=3,histtype='step')
plt.hist(ratio_rkG_rL_gamma,bins=20,range=(0,2),edgecolor='green',color='green',alpha=0.5,label='Pulsars in gamma-ray surveys only',zorder=3,histtype='step')
plt.legend()
plt.xlabel(r'$r_{\rm KG}/r_{\rm lc}$')
plt.ylabel('Number of pulsars')
plt.savefig('ratio_rkg_rl.pdf',dpi=300)
plt.close()

#Rho histogram
plt.figure(79)
plt.hist(rho_rg,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=0.5,label='Pulsars in both radio and gamma-ray surveys',zorder=3,histtype='step')
plt.hist(rho_r,bins=20,range=(0,180),edgecolor='red',color='red',alpha=0.5,label='Pulsars in radio surveys only',zorder=3,histtype='step')
plt.hist(rho_g,bins=20,range=(0,180),edgecolor='green',color='green',alpha=0.5,label='Pulsars in gamma-ray surveys only',zorder=3,histtype='step')
plt.legend()
plt.xlabel(r'$\rho$ (°)')
plt.ylabel('Number of pulsars')
plt.savefig('rho_sim.pdf',dpi=300)
plt.close()

#Ratio r_kg and R_L histogram
plt.figure(80)
plt.hist(r_kg_rg,bins=20,range=(0,31),edgecolor='blue',color='blue',alpha=0.5,label='Pulsars in both radio and gamma-ray surveys',zorder=3,histtype='step')
plt.hist(r_kg_r,bins=20,range=(0,31),edgecolor='red',color='red',alpha=0.5,label='Pulsars in radio surveys only',zorder=3,histtype='step')
plt.hist(r_kg_g,bins=20,range=(0,31),edgecolor='green',color='green',alpha=0.5,label='Pulsars in gamma-ray surveys only',zorder=3,histtype='step')
plt.legend()
plt.xlabel(r'Altitude of emission ($R_{\rm NS}$)')
plt.ylabel('Number of pulsars')
plt.savefig('Alt_em.pdf',dpi=300)
plt.close()

#Looking into the spin distribution 
P_spinup_filtered=[i for i in P_spinup2 if i <= 6.7e-3]
f_spinup_filtered=[1.0/i for i in P_spinup_filtered]
P_obsATNF_filtered=[i for i in P2 if i <= 6.7e-3]
f_obsATNF_filtered=[1.0/i for i in P_obsATNF_filtered]
plt.figure(81)
plt.hist(f_spinup_filtered,bins=10,range=(100,850),edgecolor='red',color='red',alpha=0.5,label='Simulation',histtype='step',density=True)
plt.hist(f_obsATNF_filtered,bins=10,range=(100,850),edgecolor='blue',color='blue',alpha=0.5,label='ATNF data',histtype='step',density=True)
plt.legend()
plt.xlabel(r'$f$ (Hz)')
plt.ylabel('p.d.f')
plt.savefig('histo_freq_pulsars.pdf',dpi=300)
plt.close()

#P_orb=f(Mc) obs general
plt.figure(82)
P_orb_forplot_mc,mc_for_plot_mc,compa_for_plot_mc=[],[],[]
for porb, mc, compa in zip(P_orb2, comp_mass_obs, processed_companions):
    if porb!=-1 and mc!=0:
        P_orb_forplot_mc.append(porb)
        mc_for_plot_mc.append(mc)
        compa_for_plot_mc.append(compa)
for mc, porb, compa in zip(mc_for_plot_mc, P_orb_forplot_mc, compa_for_plot_mc):
    style = companion_styles[compa]
    plt.scatter(mc, porb, color=style['color'], s=5, marker=style['marker'], label=style['label'])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(),fontsize='small')
#plt.xscale('log')
plt.yscale('log')
plt.xlim(0,1.8)
plt.ylim(4e-2,2000)
plt.xlabel(r'Companion mass ($M_{\odot}$)')
plt.ylabel(r'$P_{\rm orb}$ (days)')
plt.savefig('porb_mc_obs.pdf',dpi=300)
plt.close()

#P_orb=f(M2) sim more general
plt.figure(83)
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p) 
            for companion_group, fallback_group, m , p
            in zip(grouped_remnant_type, stellar_phase_str,Mc_f_spinup , P_orb)]
    if any(mask):
        plt.scatter(
            [m for m, flag in zip(Mc_f_spinup, mask) if flag],  # Filtered mass
            [p for p, flag in zip(P_orb, mask) if flag],        # Orbital period filtered
            label=f'{group}', marker=marker, s=5, color=colors[group]
        )
plt.yscale('log')
#plt.xscale('log')
plt.xlim(0,1.8)
plt.ylim(4e-2,2000)
plt.ylabel(r'$P_{orb}$ (days)')
plt.xlabel(r'Companion mass ($M_{\odot}$)')
plt.legend(fontsize='small')
plt.savefig('Porb_M2_plot_general.pdf',dpi=300)
plt.close()

#P_orb=f(Mc) obs spiders
plt.figure(84)
P_orb_obs_hr=[P_orb2[i]*24 for i in range(len(P_orb2))]
P_orb_forplot_mc2,mc_for_plot_mc2,compa_for_plot_mc2=[],[],[]
for porb, mc, compa in zip(P_orb_obs_hr, comp_mass_obs, processed_companions):
    if porb!=-1 and porb<=30 and mc!=0 and compa!="NS":
        P_orb_forplot_mc2.append(porb)
        mc_for_plot_mc2.append(mc)
        compa_for_plot_mc2.append(compa)
used_labels = set()
companion_styles = {
    'UL_group': {'color': 'blue', 'marker': 'o', 'label': 'Ultra light'},
    'He_group': {'color': 'green', 'marker': 's', 'label': 'He WD'},
    'MS': {'color': 'red', 'marker': '^', 'label': 'MS star'},
    'NS': {'color': 'black', 'marker': 'd', 'label': 'NS'},
    'CO': {'color': 'purple', 'marker': 'x', 'label': 'CO WD'},
    'NULL': {'color': 'yellow', 'marker': '*', 'label': 'Isolated'},
    'UNKN': {'color': 'grey', 'marker': 'h', 'label': 'Unknown companion'}
}
for mc, porb, compa in zip(mc_for_plot_mc2, P_orb_forplot_mc2, compa_for_plot_mc2):
    style = companion_styles[compa]
    label = style['label'] if style['label'] not in used_labels else None
    plt.scatter(mc, porb, color=style['color'], s=5, marker=style['marker'], label=label)
    if label is not None:
        used_labels.add(label)
#handles, labels = plt.gca().get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys(),fontsize='small')
plt.axvline(x=0.1, linestyle='--',color='black',label='Black widow/Redback boundary')
plt.text(1.5e-1, 1e0, "Redback side",color='red')
plt.text(1.5e-3, 1e0, "Black widow side", color='black')
plt.legend(fontsize='x-small')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-3,1e0)
plt.ylim(0,30)
plt.xlabel(r'Companion mass ($M_{\odot}$)')
plt.ylabel(r'$P_{\rm orb}$ (hr)')
plt.savefig('porb_mc_obs_spiders.pdf',dpi=300)
plt.close()

#Corbet diagram to check fig 20 of Tauris et al. (2013)
#P=f(P_orb)
plt.figure(85)
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p) and
            (compa=='COWD' or compa=='HeWD')
            for companion_group, fallback_group, p, compa
            in zip(grouped_remnant_type, stellar_phase_str, P_orb,remnant_type_str)]
    if any(mask):
        plt.scatter(
            [p for p, flag in zip(P_orb, mask) if flag],  # Filtered mass
            [pbis for pbis, flag in zip(P_spinup, mask) if flag],        # Orbital period filtered
            label=f'{group}', marker=marker, s=5, color=colors[group]
        )
plt.scatter(P_orb_AIC_saved,P_AIC_saved,facecolors='red',edgecolors='black',marker='o',s=20,linewidths=0.5,label='Recycled pulsars formed via AIC')
plt.axvline(x=200, linestyle='--',color='black')
plt.axhline(y=1.5e-2, linestyle='--',color='black')
plt.text(2e-2, 2e0, "II",color='black')
plt.text(2e-2, 2e-3, "I",color='black')
plt.text(6e2, 1e-1, "III",color='black')
plt.text(6e2, 2e-3, "IV",color='black')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-2,4000)
plt.ylim(1e-3,5e0)
plt.xlabel(r'Orbital period ($P_{\rm orb}$) of the binary in days')
plt.ylabel(r'Spin period $P$ (s)')
plt.legend(loc='upper center')
plt.savefig('P_Porb_plot_sim_AIC_check.pdf',dpi=300)
plt.close()

#Spider plot with obs/sim together
plt.figure(86)
x_min,x_max=1e-3,1e0
y_min,y_max=7e-1,30
#Simulations
category_to_group = {
    "Main Sequence": "Sim MS",
    "Terminal Main Sequence": "Sim MS",
    "Shell H Burning": "Sim Giant",
    "Core He Burning": "Sim Giant",
    "Terminal Core HeBurning": "Sim Giant",
    "Shell He Burning": "Sim AGB",
    "HeWD": "Sim HeWD",
    "COWD": "Sim COWD",
    "ONeWD": "Sim ONeWD",
    "NS ECSN": "Sim NS",
    "NS CCSN": "Sim NS",
    "BH": "Sim BH",
    "Isolated": "Sim Isolated"
}
markers={'Sim HeWD':'s','Sim COWD':'s','Sim ONeWD':'s','Sim NS':'s','Sim BH':'s','Sim MS':'s','Sim Giant':'s','Sim AGB':'s','Sim Isolated':'s'}
colors = {
    'Sim HeWD': 'green',
    'Sim COWD': 'purple',
    'Sim ONeWD': 'orange',
    'Sim NS': 'black',
    'Sim BH': 'cyan',
    'Sim MS': 'red',
    'Sim AGB': 'brown',
    'Sim Giant': 'royalblue',
    'Sim Isolated': 'darkgoldenrod'
}
grouped_remnant_type = [
    category_to_group[companion] if companion in category_to_group else
    category_to_group[fallback] if fallback in category_to_group else 'Unknown'
    for companion, fallback in zip(remnant_type_str, stellar_phase_str)
]
for group, marker in markers.items():
    mask = [(companion_group == group or
             (companion_group == 'Isolated' and fallback_group == group)) and
            not np.isnan(p) and
            (x_min <= m <= x_max) and
            (y_min <= p <= y_max)
            for companion_group, fallback_group, m , p
            in zip(grouped_remnant_type, stellar_phase_str,Mc_f_spinup , P_orb_hr)]
    if any(mask):
        plt.scatter(
            [m for m, flag in zip(Mc_f_spinup, mask) if flag],  # Filtered mass
            [p for p, flag in zip(P_orb_hr, mask) if flag],        # Orbital period filtered
            label=f'{group}', marker=marker, s=5, color=colors[group]
        )

#Observations
P_orb_obs_hr=[P_orb2[i]*24 for i in range(len(P_orb2))]
P_orb_forplot_mc2,mc_for_plot_mc2,compa_for_plot_mc2=[],[],[]
for porb, mc, compa in zip(P_orb_obs_hr, comp_mass_obs, processed_companions):
    if porb!=-1 and porb<=30 and mc!=0 and compa!="NS":
        P_orb_forplot_mc2.append(porb)
        mc_for_plot_mc2.append(mc)
        compa_for_plot_mc2.append(compa)
used_labels = set()
companion_styles = {
    'UL_group': {'color': 'blue', 'marker': 'o', 'label': 'Obs ultra light'},
    'He_group': {'color': 'green', 'marker': 'o', 'label': 'Obs He WD'},
    'MS': {'color': 'red', 'marker': 'o', 'label': 'Obs MS star'},
    'NS': {'color': 'black', 'marker': 'o', 'label': 'Obs NS'},
    'CO': {'color': 'purple', 'marker': 'o', 'label': 'Obs CO WD'},
    'NULL': {'color': 'yellow', 'marker': 'o', 'label': 'Obs Isolated'},
    'UNKN': {'color': 'grey', 'marker': 'o', 'label': 'Unknown (obs)'}
}
for mc, porb, compa in zip(mc_for_plot_mc2, P_orb_forplot_mc2, compa_for_plot_mc2):
    style = companion_styles[compa]
    label = style['label'] if style['label'] not in used_labels else None
    plt.scatter(mc, porb, facecolors='none', edgecolors=style['color'], s=5, marker=style['marker'], label=label)
    if label is not None:
        used_labels.add(label)

plt.axvline(x=0.1, linestyle='--',color='black',label='BW/RB boundary')
plt.text(1.5e-1, 1e0, "RB side",color='red')
plt.text(1.5e-3, 1e0, "BW side", color='black')
plt.yscale('log')
plt.xscale('log')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.ylabel(r'$P_{orb}$ (hr)')
plt.xlabel(r'Companion mass ($M_{\odot}$)')
plt.legend(fontsize='small')#,loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('Porb_M2_plot_sim_obs.pdf',dpi=300, bbox_inches="tight")
plt.close()
