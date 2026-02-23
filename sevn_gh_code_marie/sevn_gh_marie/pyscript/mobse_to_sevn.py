#!/usr/bin/env python
"""
SCRIPT: input_to_sevn1
PURPOSE: Take an input file for MOBSE and transform it to SEVN1 or SEVN2 style
USAGE;
Note: If the outputname is not given the default output will be input_sevn1.dat
REQUIRED PACKAGE:
    -pandas

V 1.1: Changed function p2a to avoid a small offset in the semimajor axis  // Giuliano Iorio
V 1.0: First version // Giuliano Iorio
"""

SEVN_V=2
input_columns={
    "id" : 0,
    "m1" : 1,
    "m2" : 2,
    "P"  : 3,
    "ecc" : 4,
    "Z"   : 5,
    "Tint" : 6
}
skiprows=1 #In MOBSE input first row usually contains the total number of objects
MIN_MASS=2.2
MAX_MASS=600.
N_RESAMPLE=1001000
Tini="zams"
SN="delayed"
Dtout="all"
Zplaceholder=True



















Header_SEVN2=("Mass_0", "Z_0", "spin_0","SN_0", "Tstart_0","Mass_1","Z_1","spin_1","SN_1","Tstart_1","a","e","Tend","Dtout")
Header_SEVN1=("Mass_0", "Mass_1", "Z_0","Z_1", "spin_0","spin_1","a","e","Tend","Tstart","dt","SN_0","SN_1","Dtout")

MOBSE_SEVN_coupling={
    "Mass_0" : "m1",
    "Mass_1" : "m2",
    "Z_0" : "Z",
    "Z_1" : "Z",
}




import sys
import numpy as np
import pandas as pd

def p2a(p,m1,m2):
    G = 3.925125598496094e8 #RSUN^3 YR^-2 MSUN^-1 (astropy: constant.G.to("Rsun^3/yr^2/Msun"))
    yeardy=365.2425
    p = p/yeardy
    a = ((p*p*G*(m1+m2))/(4*np.pi*np.pi))**(1./3.)
    return a



if __name__=="__main__":

    PLACEHOLDER="xxx"

    argv = sys.argv

    if len(argv)==1:
        raise IOError("Please write the PATH to the input file")
    elif len(argv)==2:
        input_file=sys.argv[1]
        output_name=sys.argv[1].split(".")[0]+"_sevn%i.dat"%SEVN_V
    elif len(argv)==3:
        input_file, output_name=sys.argv[1:]
    else:
        raise IOError("Maximum number of script parameters is 2, you give %i parameters"%(len(argv)-1))


    if (SEVN_V==1):
        Header=Header_SEVN1
    elif (SEVN_V==2):
        Header=Header_SEVN2
    else:
        raise IOError("SEVN_V %s unkwnown"%(str(SEVN_V)))


    if Tini is None:  Tini=PLACEHOLDER
    if SN is None:  SN=PLACEHOLDER
    if Dtout is None: Dtout=PLACEHOLDER

    mobse_input = np.loadtxt(input_file,skiprows=1)
    idx = (mobse_input[:,input_columns["m1"]]>=MIN_MASS) & (mobse_input[:,input_columns["m1"]]<=MAX_MASS) & (mobse_input[:,input_columns["m2"]]>=MIN_MASS) & (mobse_input[:,input_columns["m2"]]<=MAX_MASS)

    out_dict={}

    out_dict["Mass_0"]=mobse_input[idx,input_columns["m1"]]
    out_dict["Mass_1"] = mobse_input[idx,input_columns["m2"]]
    out_dict["spin_0"] = out_dict["spin_1"] = np.zeros_like(out_dict["Mass_0"])
    out_dict["SN_0"] = out_dict["SN_1"] =  np.array([SN,]*len(out_dict["Mass_0"]))
    out_dict["Tstart_0"] = out_dict["Tstart_1"] = np.array([Tini,]*len(out_dict["Mass_0"]))
    out_dict["Dtout"] = np.array([Dtout,]*len(out_dict["Mass_0"]))
    out_dict["a"] = p2a(mobse_input[idx,input_columns["P"]],out_dict["Mass_0"],out_dict["Mass_1"])
    out_dict["e"] = mobse_input[idx,input_columns["ecc"]]
    out_dict["Tend"] = mobse_input[idx,input_columns["Tint"]]
    if Zplaceholder:
        out_dict["Z_0"] = out_dict["Z_1"] = np.array([PLACEHOLDER,]*len(out_dict["Mass_0"]))
    else:
        out_dict["Z_0"] = out_dict["Z_1"] = mobse_input[idx,input_columns["Z"]]

    df=pd.DataFrame.from_dict(out_dict)

    if N_RESAMPLE is not None and N_RESAMPLE<len(df):
        df=df.sample(n=N_RESAMPLE)

    df[list(Header)].to_csv(output_name, sep="\t",header=False,index=False)



