####################################################################################
# Program that computes the exponents (alpha) for rydberg and continuum functions  #
####################################################################################

import math
import numpy as np

a=[0.584342,0.452615,0.382362,0.337027,0.304679] # taken from article  K Kaufmann et al 1989 J. Phys. B: At. Mol. Opt. Phys. 22 2223
b=[0.424483,0.309805,0.251333,0.251013,0.189944] # taken from article  K Kaufmann et al 1989 J. Phys. B: At. Mol. Opt. Phys. 22 2223


def alpha_c(n,l,z):
    '''
    Continuum functions 
    Generate exponents
    n - principal quantum number
    l - orbital angular momentum number
    z - net charge of ionic core
    '''
    a_l=a[l]
    b_l=b[l]
    alpha_nl=(1/(4*np.power(a_l*n+b_l,2))) 
    return  alpha_nl

def alpha(n,l,z):
    '''
    Rydberg functions
    Generate exponents
    n - principal quantum number
    l - orbital angular momentum number
    z - net charge of ionic core
    '''
    a_l=a[l]
    b_l=b[l]
    alpha_nl=np.power(z/(2*n),2)*(1/np.power(a_l*n+b_l,2)) 
    return  alpha_nl

n=[1.5,2,2.5,3,3.5]
l=2
z=1
alphas=[]


for j in n:
    alphas.append(alpha(j,l,z))

round_to_eights = [round(num, 8) for num in alphas]

for i in range(0,len(round_to_eights)):
    print(round_to_eights[i])
