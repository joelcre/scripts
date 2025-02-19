#!/usr/bin/python
#########################################################################
# Program that calculates the spatial overlap of two molecular orbitals #
#########################################################################

import sys
import math
import numpy as np


def remove_blankspaces(lists,option):

    if option==True:
        new_lists=[]
        new_list=[]

        for list in lists:
            for el in list:
                if el!='':
                    if '\n' in el:
                        el.strip('\n')
                        new_list.append(el)
                    else:
                        new_list.append(el)
    
            new_lists.append(new_list)
            new_list=[]
        return new_lists
    elif option==False:
        
        new_list=[]
        
        for el in lists:
            if el!='':
                if '\n' in el:
                    el.strip('\n')
                    new_list.append(el)
                else:
                    new_list.append(el)
    
        return new_list
        
   

'''
Function that returns the index of the molecular orbital and number of atoms directly from file
'''
def mo_info(f):

    file=open(f,"r")
    line_index=0

    for line in file:
        if line_index==1:
            mo_index=line.split(' ')[2]

        if line_index==2:
            nr_atoms=abs(int(remove_blankspaces(line.split(' '),False)[0]))
            break
        
        line_index=line_index+1

    return [mo_index,nr_atoms]

'''
 Read in a cube file and convert it into densities in
 the format of a numpy array
'''
def extract_cube_file(f,ngrid,mo_index,nr_atoms):

    file=open(f,"r")
    lines=[]
    append_condition=False

    line_index=0
    for line in file:

        if append_condition==True:
            lines.append(line.split(' '))


        if line_index>(nr_atoms+5):
            append_condition=True

        line_index=line_index+1
        
    line_s= remove_blankspaces(lines,True)

    
    list_cube_points=[]
    for i in range(0,len(line_s)):
        for j in range(0,len(line_s[i])):
            if line_s[i][j]!='\n':
                list_cube_points.append(line_s[i][j])


    cube_points=np.zeros((ngrid,ngrid,ngrid))
    for i in range(0,ngrid):
        for j in range(0,ngrid):
            for k in range(0,ngrid):
                cube_points[i,j,k]=np.abs(float(list_cube_points[(i*ngrid*ngrid)+(j*ngrid+k)]))
                
    
    file.close
    return cube_points


'''
Extract spacing in each dimension for grid in cube file
- input a cube file
- output 
'''
def construct_grid(f):
    file=open(f,"r")
    lines=[]
    append_condition=False

    counter=0
    
    for line in file:

        if counter>2 and counter<6:
            lines.append(line.split(' '))


        counter=counter+1

    dim=[]    
    line_s= remove_blankspaces(lines,True)
    
    for i in range(0,len(line_s)):
        dim.append(float(line_s[i][i+1]))
        if i==0:
            ngrid=int(line_s[0][0])
    
 
    return [ngrid,dim[0],dim[1],dim[2]]
                

'''
Function for test of integrator
'''
def func(x,y,z):
    return(x*x+y+4*z)
#    return (x+y+z)


'''
Array to store function values for test of integrator
'''
def func_arr(hx,hy,hz,nx,ny,nz):

    func_arr=np.zeros((nx,ny,nz))

    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                func_arr[i,j,k]=func(hx*i,hy*j,hz*k)

    return func_arr


'''
3D integrator using the  trapezoidal method
- func is the function to be integrated
- hx,hy,hz are the stepsizes
- nx,ny,nz are the number of elements for each dimension
'''
def integrate_trapz_3d(func,hx,hy,hz,nx,ny,nz):
    no_sum_term=func[0,0,0]+func[0,0,nz-1]+func[0,ny-1,0]+func[0,ny-1,nz-1]+func[nx-1,0,0]+func[nx-1,0,nz-1]+func[nx-1,ny-1,0]+func[nx-1,ny-1,nz-1]

    one_sum_term=(np.sum(func[0,0,1::])+np.sum(func[0,ny-1,1::])+np.sum(func[nx-1,0,1::])+np.sum(func[nx-1,ny-1,1::]))+(np.sum(func[0,1::,0])+np.sum(func[0,1::,nz-1])+np.sum(func[nx-1,1::,0])+np.sum(func[nx-1,1::,nz-1]))+(np.sum(func[1::,0,0])+np.sum(func[1::,0,nz-1])+np.sum(func[1::,ny-1,0])+np.sum(func[1::,ny-1,nz-1]))

    two_sum_term=(np.sum(func[0,:,1::])+np.sum(func[nx-1,:,1::]))+(np.sum(func[1::,0,1::])+np.sum(func[1::,ny-1,1::]))+(np.sum(func[1::,:,0])+np.sum(func[1::,:,nz-1]))

    three_sum_term=np.sum(func[1::,:,1::])
    
    integral=(((hx*hy*hz)/(8.0))*(no_sum_term+2*one_sum_term+4*two_sum_term+8*three_sum_term))
    return integral


def main():
    file_1=sys.argv[1]
    file_2=sys.argv[2]


    [ngrid,hx,hy,hz]=construct_grid(file_1)
    [ngrid2,hx2,hy2,hz2]=construct_grid(file_2)

    if ngrid!=ngrid2:
        print('Number of grid points of both files must be equal for current implementation')
        quit()
    
    [mo_1,nr_atoms_1]=mo_info(file_1)
    [mo_2,nr_atoms_2]=mo_info(file_2)

    cube_points_1=extract_cube_file(file_1,ngrid,mo_1,nr_atoms_1)
    cube_points_2=extract_cube_file(file_2,ngrid,mo_2,nr_atoms_1)

    integrand=np.multiply(cube_points_1,cube_points_2)
    integral=integrate_trapz_3d(integrand,hx,hy,hz,ngrid,ngrid,ngrid)
    print('Spatial overlap of molecular orbitals ',mo_1, ' and ', mo_2,': ',integral)



#main()


#Test of integrator
hx=3.0/500.0
hy=2.0/500.0
hz=4.0/500.0

func=func_arr(hx,hy,hz,500,500,500)

integral=integrate_trapz_3d(func,hx,hy,hz,500,500,500)
print(integral)

