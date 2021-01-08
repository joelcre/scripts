######################################################
# Script to extract  information from electric dipole#    
# polarizability  in dirac and construct a           #
# spectrum written to file                           #
######################################################

import matplotlib.pyplot as plt
from matplotlib import collections as matcoll
import math
import numpy as np

#---------------------------------------------------------------#
#---remove element from list that are empty--#
def clean_list(old_list):
    new_list=[]
    for element in old_list:
        if len(element)>0:
            new_list.append(element)
    

    return new_list


#----------------Function to extract polarizabilites from file and construct a spectra-------------#
def construct_spectrum(f):
    PI=math.pi
    C= 137.036 #Speed of light in atomic units
    
    file=open(f,"r")
    xx=[]
    yy=[]
    zz=[]
    freq=[]
    line_ptr=False

    #extract trace of the electric dipole polarizability tensor
    for line in file:
        if '<< 1, 1>>:' in line:
            xx.append((clean_list(line.split(' ')))[3])
            freq.append((clean_list(line.split(' ')))[8])
        if '<< 2, 2>>:' in line:
            yy.append((clean_list(line.split(' ')))[3])
        if '<< 3, 3>>:' in line:
            zz.append((clean_list(line.split(' ')))[3])

    # take out the imaginary/real part of the traces of the electric dipole polarizability tensor      
    xx=xx[1:len(xx):2]  # real if from 0, from 1 imaginary 
    yy=yy[1:len(yy):2]
    zz=zz[1:len(zz):2]
    freq=freq[0:len(freq):2]
    file.close

    #
    spectrum=np.zeros((len(freq),2))
    for i in range(0,len(freq)):
        spectrum[i,0]=float(freq[i])
        spectrum[i,1]=-((4*PI*spectrum[i,0])/C)*((1/3)*(float(xx[i])+float(yy[i])+float(zz[i]))) #convert to cross-section
        
    return spectrum


#----function that writes the spectrum to file---#
def write_to_file(in_file,out_file):
    
    spectrum=construct_spectrum(in_file)
    
    for i in range(0,len(spectrum[:,0])):
        out_file.write(str(spectrum[i,0])+ ' ')
        out_file.write(str(spectrum[i,1])+ ' ')
        out_file.write('\n')



files=['camb3lyp_cpp_1_5AA_HF_3c_opt_5AA_HF_3c_opt.out','camb3lyp_cpp_2_5AA_HF_3c_opt_5AA_HF_3c_opt.out','camb3lyp_cpp_3_5AA_HF_3c_opt_5AA_HF_3c_opt.out','camb3lyp_cpp_4_5AA_HF_3c_opt_5AA_HF_3c_opt.out','camb3lyp_cpp_5_5AA_HF_3c_opt_5AA_HF_3c_opt.out','camb3lyp_cpp_6_5AA_HF_3c_opt_5AA_HF_3c_opt.out']

f = open('spect.txt', 'w')

for file in files:
    write_to_file(file,f)

f.close()

