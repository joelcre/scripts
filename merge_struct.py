##########################################################################
# Script to merge water molecules from optimization with the remainder
##########################################################################

'''
1. copy optimized structures back to directory
2. read in xyz file
3. convert structure that the optimized structure was cut out from into an xyz file 
3. merge xyz files into one
4. convert merged xyz file into a pdb file
5. run through the order_residues function
6. cut out 15 Å sphere
7. remove platinum complex
8. copy structures to separate directories
'''
import math
import numpy as np
import shutil
import os
import subprocess


'''
define type, index and atoms for each residue
'''
class Residue:
    def __init__(self,residue_type,index,atoms):
        self.residue_type=residue_type
        self.index=index
        self.atoms=atoms

    def displ_atoms(self):
        [print(a.nr) for a in self.atoms]

'''
define element, coordinates and atom number for an atom
'''
class Atom:
    def __init__(self,element,x,y,z,nr):
        self.element=element
        self.x=x
        self.y=y
        self.z=z
        self.nr=nr
        
'''
Calculate the distance between two atoms
'''
def distance(atom_1,atom_2):
    x_2=np.power(atom_1.x-atom_2.x,2)
    y_2=np.power(atom_1.y-atom_2.y,2)
    z_2=np.power(atom_1.z-atom_2.z,2)
    r=np.sqrt(x_2+y_2+z_2)
    return r

'''
read in pdb-file and sort the residues into a list 
'''
def readin(file):
    f=open(file,'r')

    lines=[]

    for line in f:
        if 'HETATM' in line:
            line_split=line.split(' ')
            line_strip=[]
            [line_strip.append(el) for el in line_split if el!='']
            lines.append(line_strip)
        
    residue_ind=[] #list of the residue indices
    for line in lines:
        if int(line[4]) not in residue_ind:
            residue_ind.append(int(line[4]))


    atom_nr=0
    residues=[]
    nr_res=len(residue_ind)
    for i in range(0,nr_res):
        atoms=[]
        for line in lines:
            if int(line[4])==residue_ind[i]:
                res_type=line[3]
                x=float(line[5])
                y=float(line[6])
                z=float(line[7])
                element=line[10]
                atoms.append(Atom(element,x,y,z,atom_nr))
            
                atom_nr=atom_nr+1


        
        residues.append(Residue(res_type,residue_ind[i],atoms))
    
    return residues

'''
write cutout of the input structure (in_file) to a new file (out_file)
'''
def write_pdb_file(in_file,out_file,residue_list):

    fin=open(in_file,'r')
    #cut_file=in_file+'_15_cut'
    #fcut=open(cut_file,'w')
    fout=open(out_file,'w')
    for line in fin:
        if 'HETATM' in line:
            line_split=line.split(' ')
            line_strip=[]
            [line_strip.append(el) for el in line_split if el!='']
            if int(line_strip[4]) in residue_list: 
    
                fout.write(line) 

           # if int(line_strip[4]) not in residue_list:
            #    fcut.write(line)
                

    #fcut.close()
    fin.close()
    fout.close()

'''
write lines containing structural information of the platinum complex into separate file, then convert to xyz
'''    
def write_xyz_file(in_file,out_file,xyz_file,residue_list):
    fin=open(in_file,'r')
    fout=open(out_file,'w')
    
    for line in fin:
        if 'HETATM' in line:
            line_split=line.split(' ')
            line_strip=[]
            [line_strip.append(el) for el in line_split if el!='']
            if int(line_strip[4]) in residue_list:
    
                fout.write(line)

    fin.close()
    fout.close()
    babel_runstr="babel -i pdb " +out_file+ " -o xyz " +xyz_file 
    subprocess.run(babel_runstr, shell=True)

    
'''
calculate the maximum distance to an atomic center within the same residue 
'''
def radius_sphere(residues,center,residue_list):
    distances=[]
    for residue in residues:
        if residue.residue_type in residue_list:
            atoms=residue.atoms
            for atom in atoms:
                distances.append(distance(center,atom))

    return max(distances)

'''
identify all  atoms within a cutoff radius from center
and return a list of atoms either outside (GT) or inside (LT) the cutoff
'''
def identify_within_radius(residues,center,cutoff,operation):
    residue_list=[]
    for residue in residues:
        atoms_inside=0
        if residue.residue_type=='HOH':
            for atom in residue.atoms:
                r=distance(center,atom)
                if operation=='GT': #greater then
                    if r>cutoff:
                        atoms_inside=atoms_inside+1

                elif operation=='LT':#less then
                         if r<cutoff:
                             atoms_inside=atoms_inside+1

                if atoms_inside==3:
                    residue_list.append(residue.index)

                    
    return residue_list


def order_residue(in_file,out_file):
    fin = open(in_file, "rt")
    fout = open(out_file, "wt")

    i=1
    error=False
    for line in fin:
        
        if 'O   HOH'  in line:
            if error==False:
                i=i+1
            elif error==True:
                i=i+1
                
            shift=len(str(i))
            if shift==2:
                line=line.replace( str(i+1)+' ',str(i)+' ')
            elif shift==3:
                line=line.replace(str(i+1),str(i))
            else:
                line=line.replace(' '+str(i+1)+' ',' '+str(i)+' ' )
                
            fout.write(line)
        elif  'H   HOH'  in line:
            shift=len(str(i))
            if shift==2:
                fout.write(line.replace(' 0 ', str(i)+' '))
            elif shift==3:
                fout.write(line.replace('  0 ', str(i)))
            else:
                fout.write(line.replace(' 0 ', ' '+str(i)+' '))
        

        elif 'H   UNK' in line:
            #fout.write(line.replace('ATOM', 'HETATM'))
            
            line=line.replace('ATOM', 'HETATM')
            line=line.replace('UNK','HOH')
            
            shift=len(str(i+1))
            if shift==2:
                line=line.replace( str(i+1)+' ',str(i)+' ')
            elif shift==3:
                line=line.replace(str(i+1),str(i))
            else:
                line=line.replace(' '+str(i+1)+' ',' '+str(i)+' ' )
            fout.write(line)
            error=True

        else:
            fout.write(line)


    fin.close()
    fout.close()



def merge_xyz(opt_file,old_file,new_file):
    fopt=open(opt_file,'r')
    fold=open(old_file,'r')
    fnew=open(new_file,'w')

    lines_1=fopt.readlines()
    lines_2=fold.readlines()
    nr_atoms=int(lines_1[0])+int(lines_2[0])
    merge_str=opt_file+old_file+'\n'

    fnew.write(str(nr_atoms)+'\n')
    fnew.write(merge_str)
    lines=lines_1[2:]+lines_2[2:]
    
    for i in range(0,len(lines)):
        fnew.write(lines[i])

    fopt.close()
    fold.close()
    fnew.close()

    
'''
'''
def make_cutouts(in_file,out_file,radius):
    residues=readin(in_file)
    
    for f in residues:
        for atom in f.atoms:
            if atom.element=='Pt':
                center=atom


    residue_list=['LIG']
    dist=radius_sphere(residues,center,residue_list)
    cutoff=dist+float(radius)
    res_ind=identify_within_radius(residues,center,cutoff,'LT')
    res_ind.append(1)
    #[res_ind.append(i) for i in range(1,8)]
    
    #out_file='test_'+in_file
    write_pdb_file(in_file,out_file,res_ind)
    

def cutout_pt(opt_file,pt_file):
    fin=open(opt_file,'r')
    fout=open(pt_file,'w')

    lines=fin.readlines()
    for i in range(0,21):
        line=lines[i]
        fout.write(line)

    fin.close()
    fout.close()


'''
convert unaltered strucutures from which the cutouts were taken into xyz format
'''

struct_indexes=[117,127,141,151,161,171,181,189]
file_prefix='opt_'
for i in struct_indexes:
    babel_runstr="babel -i xyz opt_" +str(i)+'.xyz' " -o pdb " +"opt_" +str(i)+'.pdb' #convert to xyz-file
    subprocess.run(babel_runstr, shell=True)
    order_residue("opt_" +str(i)+'.pdb',"opt_" +str(i)+'ordered.pdb')


for j in struct_indexes:
    opt_file='opt_'+str(j)+'.xyz'
    pt_file='Pt_'+str(j)+'.xyz'
    cutout_pt(opt_file,pt_file)
