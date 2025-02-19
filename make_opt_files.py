################################################################################################
# Script to cutout a sphere of molecules within a certain distance from a defined atomic center
# and then construct orca geometry input files,run-files and directories
################################################################################################
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
        if 'ATOM' in line:
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
    cut_file=in_file+'_cut'
    fcut=open(cut_file,'w')
    fout=open(out_file,'w')
    for line in fin:
        if 'ATOM' in line:
            line_split=line.split(' ')
            line_strip=[]
            [line_strip.append(el) for el in line_split if el!='']
            if int(line_strip[4]) in residue_list:
    
                fout.write(line)

            if int(line_strip[4]) not in residue_list:
                fcut.write(line)
                

    fcut.close()
    fin.close()
    fout.close()

'''
Makes an orca input file for geometry optimization  
'''
def write_input_file(input_file,atom_list,xyzfile):
    
    fout=open(input_file,'w')
    struct_str='*xyzfile 0 1 '+xyzfile +'\n'

    fout.write('! PBEh-3c Opt Freq PAL8\n') #BP86 RI def2-SV(P) def2/J
    fout.write(struct_str)
    fout.write('%maxcore 5000\n')
    fout.write('%geom\n')
    fout.write('Maxiter 600\n')
    fout.write('Constraints\n')

    for a in atom_list:
        c_line='{C '+str(a) + ' C}\n'
        fout.write(c_line)

    fout.write('end\n')
    fout.write('end\n')
    fout.close()
    
'''
write batch submit script of geometry optimization in orca
'''
def write_run_file(run_file,input_file,output_file):
    fout=open(run_file,'w')
    fout.write('#!/bin/sh\n')
    fout.write('#SBATCH -N 1\n')
    fout.write('#SBATCH --tasks-per-node=8\n')
    fout.write('#SBATCH -t 168:00:00\n')   
    fout.write('#SBATCH --mem-per-cpu=5000\n')
    fout.write('')
    fout.write('module purge\n')
    fout.write('ml GCC/8.3.0  OpenMPI/3.1.4\n')
    fout.write('ml ORCA/4.2.1\n')
    run_str='/sw/easybuild/software/ORCA/4.2.1-gompi-2019b/orca '+input_file+' > '+output_file 
    fout.write(run_str)
    fout.close()

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
        if residue.residue_type=='WAT':
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

'''
return list of atoms to be frozen in a geometry optimization 
'''
def freeze_atom_nr(residues,res_ind):

    centers_freeze=[] #lists the atom number that are to be kept frozen 
    
    for ind in res_ind:
        for residue in residues:
            if residue.index==ind:
                for atom in residue.atoms:
                    centers_freeze.append(atom.nr)

    return centers_freeze

'''
Make files and directories for geomtry optimizations
'''
def make_opt(pdb_file,input_file,xyz_file,run_file,output_file,run_dir,radius):

    residues=readin(pdb_file)
    for f in residues:
        if f.residue_type=='PT1':
            center=f.atoms[0]
    residue_list=['AM1','AM2','AZ1','AZ2','HD1','HD2']
    dist=radius_sphere(residues,center,residue_list)
    cutoff=dist+radius
    res_ind=identify_within_radius(residues,center,cutoff,'GT')
    atom_list=freeze_atom_nr(residues,res_ind)
    write_input_file(input_file,atom_list,xyz_file)
    write_run_file(run_file,input_file,output_file)
    #create directory for calculation
    #move files to directory
    currDir = os.getcwd()
    inp_src=currDir+'/'+input_file
    inp_dest=currDir +'/'+run_dir+'/'+input_file
    shutil.move(inp_src,inp_dest)
    shutil.move(currDir+'/'+xyz_file,currDir+'/'+run_dir+'/'+xyz_file)
    shutil.move(currDir+'/'+run_file,currDir+'/'+run_dir+'/'+run_file)


'''
'''
def make_cutouts(in_file,radius):
    residues=readin(in_file)
    
    for f in residues:
        if f.residue_type=='PT1':
            center=f.atoms[0]
            
    residue_list=['AM1','AM2','AZ1','AZ2','HD1','HD2']
    dist=radius_sphere(residues,center,residue_list)
    cutoff=dist+radius
    res_ind=identify_within_radius(residues,center,cutoff,'LT')
    [res_ind.append(i) for i in range(1,8)]
    
    out_file='cutout_'+in_file
    write_pdb_file(in_file,out_file,res_ind)
   


'''
make cutouts of pdb-files from md snapshots
'''

file_prefix='mask.pdb.'
files=[]
md_snapshots=range(1,200,2)
for i in range(1,200,2):
    temp_file=file_prefix+str(i)
    files.append(temp_file)
    
for f in files:
    print(f)
    make_cutouts(f,6.0)

    
'''
Construct filenames for xyz, input pdb files and orca run-files
'''    

pdb_files=[]
pdbfile_prefix='cutout_mask.pdb.'
xyz_files=[]
xyzfile_prefix='cutout_mask.xyz.'
input_files=[]
run_files=[]
output_files=[]
for i in range(1,200,2):
    pdb_file_str=pdbfile_prefix+str(i)
    xyz_file_str=xyzfile_prefix+str(i)
    input_file_str='opt_'+str(i)+'.inp'
    output_file_str='opt_'+str(i)+'.out'
    run_file_str='run_'+str(i)+'.sh'
    babel_runstr="babel -i pdb cutout_mask.pdb." +str(i)+ " -o xyz " +"cutout_mask.xyz."+str(i) #convert to xyz-file
    subprocess.run(babel_runstr, shell=True)
    pdb_files.append( pdb_file_str)
    xyz_files.append(xyz_file_str)
    input_files.append(input_file_str)
    run_files.append(run_file_str)
    output_files.append(output_file_str)


'''
Make input files,run files and directories for all cutouts 
'''

for i in range(0,len(pdb_files)):
    os.mkdir('opt_'+str(i+1))
    make_opt(pdb_files[i],input_files[i],xyz_files[i],run_files[i],output_files[i],'opt_'+str(i+1),4.0)



