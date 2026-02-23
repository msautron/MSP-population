from astropy.table import Table 
import math
import os 
import random
import argparse

parser = argparse.ArgumentParser(description="Id sim")
parser.add_argument("Id_sim", type=int)
args = parser.parse_args()
Id_sim_int=args.Id_sim
Id_sim=str(Id_sim_int)

a=-1.0
e=-1.0
i=0

data_binary=Table.read(f'/home/matteo.sautron/sevn_gh/SEVN_run/sevn_output{Id_sim}/output_0.csv',format='csv')

M_companion_init=data_binary['Mass_1'][0]
if (math.isnan(data_binary['Mass_1'][-1])):
    M_companion_final=-1.0
else :
    M_companion_final=data_binary['Mass_1'][-1]
if (math.isnan(data_binary['Semimajor'][-1]) and data_binary['RemnantType_1'][-1]<4):
    a=-1.0
elif (data_binary['RemnantType_1'][-1]>=4 and math.isnan(data_binary['Semimajor'][-1])):
    while (a==-1.0 or i<10):
        a_test=data_binary['Semimajor'][-1-i]
        if (math.isnan(a_test)):
            a=-1.0
            i+=1
        else : 
            a=a_test
else :
    a=data_binary['Semimajor'][-1]
i=0
if (math.isnan(data_binary['Eccentricity'][-1]) and data_binary['RemnantType_1'][-1]<4):
    e=-1.0
elif (data_binary['RemnantType_1'][-1]>=4 and math.isnan(data_binary['Eccentricity'][-1])):
    while (e==-1.0 or i<10):
        e_test=data_binary['Eccentricity'][-1-i]
        if (math.isnan(e_test)):
            e=-1.0
            i+=1
        else : 
            e=e_test
else :
    e=data_binary['Eccentricity'][-1]
i=0
star_type=float(data_binary['Phase_1'][-1])
remnant_type=float(data_binary['RemnantType_1'][-1])

dict_star_stage = {
        0:"PreMainSequence",
        1:"MainSequence",
        2:"TerminalMainSequence",
        3:"ShellHBurning",
        4:"CoreHeBurning",
        5:"TerminalCoreHeBurning",
        6:"ShellHeBurning",
        7:"Remnant"
}

dict_remnant = {
        0:"NotARemnant",
        1:"HeWD",
        2:"COWD",
        3:"ONeWD",
        4:"NS_ECSN",
        5:"NS_CCSN",
        6:"BH",
        -1:"Empty"
}

with open(f'/home/matteo.sautron/sevn_gh/SEVN_run/data_pop_ms{Id_sim}.txt', 'a') as f:
    content_to_add = f"{M_companion_init} {M_companion_final} {a} {e} {star_type} {remnant_type}\n"
    f.write(content_to_add)

#try:
#    with open(f'/home/matteo.sautron/sevn_gh/SEVN_run/data_pop_ms{Id_sim}.txt', "r") as f:
#        lines = f.readlines()
#except FileNotFoundError:
#    print("The file you are looking for does not exist")
#    exit(1)

#last_line = lines[-1].rstrip()

#with open(f'/home/matteo.sautron/sevn_gh/SEVN_run/data_pop_ms{Id_sim}.txt','w') as f:
#    f.writelines(lines[:-1])
#    content_to_add=f"{M_companion_init} {M_companion_final} {a} {e} {star_type} {remnant_type}\n"
#    # Rewrite last line with the added content
#    f.write(last_line +" "+ content_to_add)
