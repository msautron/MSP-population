# Examine or display the 3PC sensitivity map. 
# Reproduce 3PC Fig 25.
# Obtain sensitivity limit for a given point on the sky.
# Toby Burnett, July 2023.

import matplotlib.pyplot as plt
import healpy as hp
import re
import numpy as np
import os
filename='3PC_SensitivityMap_20230629.fits'

# Generate plot similar to 3PC Figure 25.
sensitivity_map = hp.read_map(filename)
#hp.mollview(sensitivity_map, norm='log', title=filename, 
            #unit=r'$\mathrm{G_{100}\ limit\ (erg\  cm^{-2}\ s^{-1})}$', 
            #xsize=1200, cmap='plasma', badcolor='white',min=4e-13, max=2e-12);
#plt.show()

#Remove old data
os.system(f'rm fermi_fmin.txt')

# Load the coordinates data
reg_4=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
with open("l_b_coord_sim.txt","r") as f:
    data_lb=re.findall(reg_4,f.read())

number_of_sources=int(len(data_lb)/2)

l,b=[],[]
for i in range(number_of_sources):
    l.append(float(data_lb[2*i]))
    b.append(float(data_lb[2*i+1]))

def get_lb(name):
    """ funcion which uses astropy to look up source names 
    returns Galactic longitude, l, and latitude, b, of known source, in degrees.
    """
    from astropy.coordinates import SkyCoord
    sc = SkyCoord.from_name(name).galactic
    return sc.l.deg, sc.b.deg

def sensitivity(l,b):
    """
    function that returns the sensitivity map value at (l,b):
    Sensitivity is the phase-integrated (that is, point source) 95% Confidence Level integral energy flux, integrated above 100 MeV,
    in units of erg/cm2/s.
    """
    return sensitivity_map[hp.ang2pix(hp.get_nside(sensitivity_map),l,b, lonlat=True)]

# demonstrations
#print( f'Sensitvity at  (270,0.3) is {sensitivity (-90,0.3):.2e} erg s-1 cm-2')
for i in range(number_of_sources):
    s=sensitivity(l[i],b[i])*1e-3 #Conversion in W.m^-2
    if (np.isnan(s)):
        s=1.0e55
    with open('fermi_fmin.txt','a') as f:
        f.write(f'{s}\n')

