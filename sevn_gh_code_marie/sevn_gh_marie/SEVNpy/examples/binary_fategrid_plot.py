"""
Simple script to produce a grid plot similar to Fig. 2 in Gallegos-Garcia+21 paper
(https://ui.adsabs.harvard.edu/abs/2021ApJ...922..110G/abstract)
"""
from sevnpy.sevn import SEVNmanager,sevnwrap
from sevnpy.io import readlogstring
import matplotlib.pyplot as plt
import numpy as np

# Mass and metallicity of the donor star for which we want to generate the gris
Mdonor = 25
Zdonor = 0.002

# SEVN option for mass transfer stability, use qcrit_startrack to be consistent with Gallegos-Garcia+21
qcrit_option = "qcrit_startrack"
# SEVN alpha option for the CE,  use qcrit_startrack to be consistent with Gallegos-Garcia+21
alpha_CE = 1

# mass-ration and Period grid
qgrid = np.linspace(0.05,1,50)
Pgrid = np.logspace(0.0,4,50)


def find_channel(log):
    """
    Find the formation channel based on the SEVN log
    """
    alog=readlogstring(log,events=["CE","RLO_BEGIN","MERGER"],capturing_cols=("time","event"))
    dfCE=alog["CE"]
    dfRL=alog["RLO_BEGIN"]
    dfM=alog["MERGER"]

    if len(dfRL)==0 and len(dfCE)==0 and len(dfM)==0:
        return 0
    elif len(dfRL)>0 and len(dfCE)==0 and len(dfM)==0:
        return 1
    elif len(dfCE)>0 and len(dfM)==0:
        return 2
    elif len(dfCE)>0 and len(dfM)>0:
        return 3
    elif len(dfM)>0:
        return 4

def find_status(df):
    """
    Find the final status of the  binary  based on the results of the SEVN evolution
    (-1 Destroyed binary, 1 Merge within an Hubble time, 0 Not merge within an Hubble time)
    """
    if np.isnan(df["GWtime"][-1]):
        return -1
    else:
        tdel= df["Worldtime"][-1] + df["GWtime"][-1]
        if tdel<=14000:return 1
        else: return 0

def find_label(df,log):
    """
    Find the label for the grid plot
    """
    ch = find_channel(log)
    st = find_status(df)

    if ch>=3:
        return "CE"
    elif ch==2 and st==1:
        return "CEm"
    elif ch==1 and st==1:
        return "MT"
    else:
        return "wide"


def calc_a(per,m1,m2):
    """
    From period to semi-major axss
    :param per:  Periods in  year
    :param m1: Mass in Msun
    :param m2:  Mass in Msun
    :return: semi-major axis in Rsun
    """
    G = 3.925125598496094e8
    return (per*per*(G*(m1+m2))/(np.pi*np.pi*4))**(1./3.)

if __name__=="__main__":

    SEVNmanager.init({"ce_alpha":alpha_CE, "rlo_stability":qcrit_option})

    ql=[] # To store the mass-ratios
    pl=[] # To store the periods
    ll=[] # To store the label
    for p in Pgrid:
        for q in qgrid:
            m1=Mdonor
            m2=m1*q
            a=calc_a(p/365.25,m1,m2)
            #Evolve a binary initialising the second star as a BH
            df,ldic=sevnwrap.evolve_binary(a,0.,m1,Zdonor,m2,Zdonor,star_flag_1="BH")
            ll.append(find_label(df,ldic["Log"]))
            ql.append(q)
            pl.append(p)
    SEVNmanager.close()

    ql =  np.array(ql)
    pl =  np.array(pl)
    ll =  np.array(ll)


    fig,ax=plt.subplots(1,1,figsize=(5.5,6))
    idx=ll=="CE"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="gold",marker="s",label="merger during CE")
    idx=ll=="wide"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="lightskyblue",marker="s",label="wide binary")
    idx=ll=="MT"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="blue",marker="s",label="stable MT to BBH merger")
    idx=ll=="CEm"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="orange",marker="s",label="BBH merger following CE")
    plt.legend(ncols=2,bbox_to_anchor=(1.0,1.11),frameon=False)
    plt.xlabel("$q$ [Macc/Mdonor]",fontsize=15)
    plt.ylabel("$\log_{10} P_\mathrm{orb,i}$ [days]",fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    fig.tight_layout()
    plt.show()


