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

reg_int = re.compile(r"-?\d+")
#reg_5=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
#reg_sci = re.compile(r"-?\d+\.\d+e[+-]?\d+")
#reg_14 = re.compile(r"^((?:-?\d+\.\d+e[+-]?\d+\s+){13}-?\d+\.\d+e[+-]?\d+)")
reg_5 = re.compile(r"-?\d\.\d{6}[eE][+-]\d{2}")

with open("data_subms_all.txt","r") as f:
    data=re.findall(reg_int,f.read())

with open("data_pop_ms_all.txt","r") as f:
    data_msp=re.findall(reg_5,f.read())

P,Pdot=[],[]
for i in range(int(len(data_msp)/14)):
    if data_msp[14*i]=='inf':
        P.append(1e55)
    else:
        P.append(float(data_msp[14*i]))
    if data_msp[14*i+1]=='inf':
        Pdot.append(1e55)
    else:
        Pdot.append(float(data_msp[14*i+1]))

print(len(P))
print(P[-1])

index_saved,P_saved,Pdot_saved=[],[],[]
for i in range(len(data)):
    data[i]=float(data[i])
    if data[i]==1:
        index_saved.append(i)
        P_saved.append(P[i])
        Pdot_saved.append(Pdot[i])

sum_allmsp=sum(data)

print(f'Number of pulsars that went under the millisecond : {sum_allmsp}')
print(f'Characteristics of the pulsars that went under the millisecond:')
#print(P_saved)
#print(Pdot_saved)
