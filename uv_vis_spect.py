###########################################################   
# Script to extract energies and oscillator strengths     #
# from ORCA output file and compute a uv-vis spectrum     #                                                         
###########################################################     

import matplotlib.pyplot as plt
from matplotlib import collections as matcoll
import math
import numpy as np

import matplotlib.ticker
#----------------Functions for convolution-------------------------------------#
def grid(energies,size):

    max_en=max(energies)+60 #extending grid by adding 60
    min_en=min(energies)-60
    dE=float((max_en-min_en)/size)
    E_grid=np.zeros(size)

    for i in range(0,size):
        E_grid[i]=min_en +float(i*dE)

    
    return E_grid


def gaussian(x,e,sigma):
      
    const=1.3062974e+8 # constant obtained from gaussian homepage
    G=(const/sigma)*np.exp(-(np.power(1/x-1/e,2))/(np.power(sigma,2)))
    
    return G


def gaussian_convolution(X,e,T,sigma):
 
  
    G_spect=np.zeros(len(X))
    G_spect= (T/(1e+7))*(gaussian(X, e,sigma))
   # print(np.trapz(gaussian(X,e,sigma)))
    return G_spect
          

#---------------------------------------------------------------#

#--------------help-functions-----------------------------------#
            

def strip_blankspace(lines):

#---stripping lines of blank-spaces----#
    c_lines=[]
    c_list=[]

    for line in lines:
        for c in line:
            if c!='':
                if '\n' in c:
                    c.strip('\n')
                    c_list.append(c)
                else:
                    c_list.append(c)
    
        c_lines.append(c_list)
        c_list=[]
    return c_lines
    return energies


def line_counter(file):
    #---counting the total number of lines in a file#
    fil=open(file,"r")
    f=fil.readlines()
    max_lines=len(f)
    fil.close()
    return max_lines


#----------------Functions to extract information from file-------------#
def extract_spectrum(f):

    file=open(f,"r")
    lines=[]
    h=-1
 
    for line in file:
    
        if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
            h=1
        elif 'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS' in line:
            h=-1

        if h==1:
            lines.append(line.split(' '))

            

    line_s= strip_blankspace(lines)


    spectrum=np.zeros((len(line_s)-7,2))
    for i in range(5,len(line_s)-2):
        spectrum[i-5,1]=float(line_s[i][3]) #oscillator strengths
        spectrum[i-5,0]=float(line_s[i][2]) # wavenumbers in nanometer


    file.close
    

    return spectrum

#--------------------------------------------------------------------#


def spectrum_gaussian(wav_len,intensities,sigma,size,color,scale_fact,labl):
        #-----constructing spectrum-----#
        const=1.3062974e+8/1e+7  # constant obtained from the gaussian hompage
         
        X=grid(wav_len,size) # constructs grid 
        min_en=min(wav_len)
        G_spect=np.zeros(size)
            
        for i in range(0,len(wav_len)):  
            G=gaussian_convolution(X,wav_len[i],intensities[i],sigma)
            G_spect=G_spect+G
        ax1.plot(X,G_spect,color,label=labl)

                
        for i in range(0,len(wav_len)): # constructs stick-spectra for wavelengths
            en=[wav_len[i],wav_len[i]]
            intens=[0,intensities[i]]
            ax3.plot(en,intens,color)
             
        ax3.scatter(wav_len,intensities,color=color)

      

    


def plot_spectra(file,color,label,broadening):
    
    spect=extract_spectrum(file)               
    wavenumbers=spect[:,0]
    oscillator_strengths=spect[:,1]
       
   
    spectrum_gaussian(wavenumbers,oscillator_strengths,broadening,1000,color,1.0,label)



def nm_t_ev(arr): #converts array from nanometers to electronvolts

    h=6.62607004e-34 # plancks constant                                                                                                                                                        
    c=299792458 # speed of light                                                                                                                                                               
    j_to_ev=1/1.60218e-19
    ev_to_j=1.60218e-19 # conversion of joule to electronvolts                                                                                                                                 
    hartree_to_ev=27.2116
    energies=((1.0/((arr[:])*(1e-9)))*(h*c))*j_to_ev

    return energies

    

#---files---#
fil1="ex1.out"
#fil2="tddft_triplet.out"
#-----------#

#---Setting axes--------#
fig, ax1 = plt.subplots()
ax1.set_ylabel('$\epsilon$ (L mol$^{-1}$cm$^{-1}$)', color='k',size=14)
ax1.set_xlabel('Wavelength(nm)',size=14)
ax1.tick_params('y', colors='k')


ax2=ax1.twiny()
ax2.set_xlabel('Energy(eV)',size=14)
ax2.tick_params('y', colors='k')
ax2.invert_xaxis()


ax3 = ax1.twinx()
ax3.set_ylabel('Osc.str', color='k',size=14)
ax3.tick_params('y', colors='k')


fig.tight_layout()



#-------main--------------------#

#---constants-----#
h=6.62607004e-34 # plancks constant                                                                                                                                                                   
c=299792458 # speed of light                                                                                                                                                                    
ev_to_j=1.60218e-19 # conversion from electronvolts to joule                                                                                                                                           
#-------------------#

sigma=0.1 # broadening in electronvolts 
broadening=1/((h*c/(ev_to_j*sigma))/1e-9)
   
plot_spectra(fil1,'b','',broadening) #add label
#plot_spectra(fil2,'r','triplet_state',broadening) #YUNHAO

ax1.legend( loc='upper right',bbox_to_anchor=(1, 1),prop={'size': 12})



# get the primary axis x tick locations in plot units
xtickloc = ax1.get_xticks()
# avoid conversion to an infinite energy value
if xtickloc[0]==0:
    xtickloc=xtickloc[1:len(xtickloc)]

# set the second axis ticks to the same locations
ax2.set_xticks(xtickloc)
# calculate new values for the second axis tick labels, format them, and set them    
x2labels = ['{:.3g}'.format(x) for x in nm_t_ev(xtickloc)]
ax2.set_xticklabels(x2labels)
# force the bounds to be the same
ax2.set_xlim(ax1.get_xlim())


plt.show()


