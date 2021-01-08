########################################################################
# Script to plot CPP spectra from dirac                                #
#                                                                      #                                                                      
########################################################################

import matplotlib.pyplot as plt
from matplotlib import collections as matcoll
import math
import numpy as np
#---------------------------------------------------------------#

def plot_spect(file,legend,plot_style,hline):

    f=open(file,"r")
    lines=[]
    for line in f:
        lines.append(line.split(' '))

        spect=np.zeros((2,len(lines)))
    for i in range(0,len(lines)):
        spect[1,i]=float(lines[i][1])
        spect[0,i]=float(lines[i][0])

    
        
    plt.plot(27.211*spect[0,:],spect[1,:],plot_style,label=legend)
    if hline==True:
        plt.plot(27.211*spect[0,:],np.zeros((len(lines))),'--k',linewidth=0.5)


fil1='spect.txt'
fil2='aq_5AA.txt'
#fil2='ll_ecd_spect.txt'

plot_spect(fil1,'4C','--.b',False)
plot_spect(fil2,'4C','--.r',False)
#plot_spect(fil2,'LL','--.r',True)
#plt.title('3R-iodo-1-butyne',size=14)
#plt.legend((),('label1', 'label2'),loc='upper right',bbox_to_anchor=(1, 1))                                                                                                                               
plt.legend(loc='upper right',bbox_to_anchor=(1, 1),prop={'size': 12})
plt.xlabel('Photon energy (eV)',size=14)
#plt.ylabel('degrees/(dm g/cm$^3$)',size=14) #optical rotation(No London)                                                                                                                                  
plt.ylabel('$\sigma(\omega)$',size=14) #cross section                                                                                                                       
#plt.grid(True)                                                                                                                                                                                            

#plt.ylim([-70,70])
#plt.xlim([120,220])
plt.show()

