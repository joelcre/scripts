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
def construct_spectrum_dalton(f):
    PI=math.pi
    C= 137.036 #Speed of light in atomic units
    file=open(f,"r")
    lines=[]
    X_pol=[]
    Y_pol=[]
    Z_pol=[]

    for line in file:
        if  'XDIPLEN   XDIPLEN' in line:  
            X_pol.append(line.split(' '))
        if  'YDIPLEN   YDIPLEN' in line:  
            Y_pol.append(line.split(' '))
        if  'ZDIPLEN   ZDIPLEN' in line:  
            Z_pol.append(line.split(' '))

    for i in range(0,len(X_pol)):
        X_pol[i]=clean_list(X_pol[i])
        Y_pol[i]=clean_list(Y_pol[i])
        Z_pol[i]=clean_list(Z_pol[i])

    spectrum=np.zeros((len(X_pol),2))
    for j in range(0,len(X_pol)):
        spectrum[j,0]=float(X_pol[j][3]) #energy
        spectrum[j,1]= ((4*PI*spectrum[j,0])/C)*(1.0/3.0)*(float(X_pol[j][5])+float(Y_pol[j][5])+float(Z_pol[j][5])) #
     
    file.close
    return spectrum
    
def construct_spectrum_dirac(f,nr_freq):
    PI=math.pi
    C= 137.036 #Speed of light in atomic units
    
    file=open(f,"r")
    T_im=np.zeros((nr_freq,3)) #imaginary part of polarizability tensor
    freq=np.zeros((3))
    line_ptr=0
    lines=file.readlines()
    #extract trace of the electric dipole polarizability tensor
    for line in lines:
        if '<<A( 1),B( 1)>>' in line:
            temp_lines=lines[line_ptr+7:line_ptr+7+nr_freq]  
            for i in range(0,len(temp_lines)):
                temp_lines[i]=clean_list(temp_lines[i].split(' '))
                freq[i]=float(temp_lines[i][0])
                T_im[i,0]=float(temp_lines[i][4])
                
        if '<<A( 2),B( 2)>>' in line:
            temp_lines=lines[line_ptr+7:line_ptr+7+nr_freq]  
            for i in range(0,len(temp_lines)):
                temp_lines[i]=clean_list(temp_lines[i].split(' '))
                T_im[i,1]=float(temp_lines[i][4])
                
        if '<<A( 3),B( 3)>>' in line:
            temp_lines=lines[line_ptr+7:line_ptr+7+nr_freq]  
            for i in range(0,len(temp_lines)):
                temp_lines[i]=clean_list(temp_lines[i].split(' '))
                T_im[i,2]=float(temp_lines[i][4])
           
        line_ptr=line_ptr+1

    

    file.close

    
    spectrum=np.zeros((len(freq),2))
    for i in range(0,len(freq)):
        spectrum[i,0]=freq[i]
        spectrum[i,1]=-((4*PI*spectrum[i,0])/C)*((1/3)*(T_im[i,0]+T_im[i,1]+T_im[i,2])) #convert to cross-section
    return spectrum
    

#----function that writes the spectrum to file---#
def write_to_file(in_file,out_file,nr_freq,program):

    if program=='dirac':
        spectrum=construct_spectrum_dirac(in_file,nr_freq)
    elif program=='dalton':
        spectrum=construct_spectrum_dalton(in_file)
          
    for i in range(0,len(spectrum[:,0])):
        out_file.write(str(spectrum[i,0])+ ' ')
        out_file.write(str(spectrum[i,1])+ ' ')
        out_file.write('\n')



files_dalton=['dal_test.out']

f = open('spect.txt', 'w')

for file in files_dalton:
    write_to_file(file,f,3,'dalton')

f.close()

