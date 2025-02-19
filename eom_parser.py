
#######################################################################################################################
# Script that extracts a state with a specific character (amplitude transitions) from a series of qchem output files  #
#######################################################################################################################

import pandas as pd

'''
A state is a class that contains the attributes energy (real and imaginary), amplitudes and the angle
'''
class state:
    def __init__(self,re_ener_,im_ener_,amplitudes_,irrep_,angle_):
        self.re_ener=re_ener_
        self.im_ener=im_ener_
        self.amplitudes=amplitudes_
        self.irrep=irrep_
        self.angle=angle_

    def print_state(self):
        print('Irrep: ', self.irrep)
        print('Im(E): ', self.re_ener)
        print('Re(E): ', self.im_ener)
        print('Angle: ', self.angle)
        print('Amplitudes: ',self.amplitudes)

    def get_im(self):
        return self.im_ener

    def get_re(self):
        return self.re_ener

    def get_angle(self):
        return self.angle

'''
A dictionary of the most common point groups and their irreducible representations
'''
point_groups={
    'D2h':['Ag','B1g','B2g','B3g','Au','B1u','B2u','B3u'],
    'C2v':['A1','A2','B1','B2']
    }


'''
returns a dictionary with the number of states per irrep
'''
def create_states_per_sym(pg,states_list):
    states_per_sym={}
    irreps=point_groups[pg]
  
    for i in range(0,len(states_list)):
        states_per_sym.update({irreps[i]:states_list[i]})

    return states_per_sym



'''
Search file for states within an irrep and return a list with state objects containing
information about energy, amplitudes and irrep
'''            
def extract_states_file(qchem_file,states_per_sym):

    states=[]
    angle=get_angle(qchem_file)

    for irrep in states_per_sym.keys():
        if states_per_sym.get(irrep)!=0:
            nr_states=states_per_sym.get(irrep)

            for i in range(1,nr_states+1):
                search_str='CS/CAP-EOMEE-CCSD transition '+ str(i)+'/'+irrep
                state=get_state(search_str,qchem_file,irrep,angle) # returns a state object containing energy,amplitudes and angle
                states.append(state)


    return states


def get_angle(qfile):
    
    for line in qfile:
        if 'complex_theta =' in line:
            angle=int(int(line.split(' ')[2].strip('\n'))/10)

    return angle
    

'''
extract amplitudes for a state and returns the information into a dataframe
'''
def extract_amplitude(qchem_file,line_nr):
    
    for i in range(line_nr,len(qchem_file)):
        if 'Amplitude    Transitions between orbitals' in qchem_file[i]:

            init=i

        if 'Summary of significant orbitals:' in qchem_file[i]:
            end=i
            break

    amplitudes_unformatted=qchem_file[init+1:end-1]
    amplitudes_formatted=[]
    for a in amplitudes_unformatted:
        amplitudes_formatted.append(strip_blankspaces(a.split('  ')))

    df=pd.DataFrame(amplitudes_formatted)
    df2 = df[[0,1,3]]
 
    
    return df2

'''
extracts the energy (real and imaginary) and the amplitudes for a state, and returns the instance of a state object
'''
def get_state(state_str,qfile,irrep,angle):

    for i in range(0,len(qfile)):
        if state_str in qfile[i]:
            temp=[]
            re_e=qfile[i+1].split(' ')[13].strip('(')
            im_e=qfile[i+1].split(' ')[15].strip(')')
            line_nr=i

    amplitudes=extract_amplitude(qfile,line_nr)
    state_=state(re_e,im_e,amplitudes,irrep,angle)
    return state_

'''
Find the point group in qchem file
'''
def find_pg(file_lines):
    irrep=0
    for pg in point_groups.keys():
        for line in file_lines:
            if pg in line:
                irrep=pg

    return irrep

'''
Removes all the empty elements from a list
'''
def strip_blankspaces(in_list):
    out_list=[]

    for el in in_list:
        if el!='':
            out_list.append(el)

    return out_list

'''
find the input lines that specifies the number of states in the input of a qchem file
'''
def find_states(qfile):
 state_list=[]
 for line in qfile:
     if 'ee_states =' in line:
         temp_list=strip_blankspaces(line.split(' '))[2].strip('[').strip(']\n').split(',')

         
 [state_list.append(int(el)) for el in temp_list]
 return state_list
    


'''
 Reads in a qchem output file containing series of jobs and divides it into a dictionary with the keys being the number of file and 
the values being the file output
'''
def divide_to_files(in_file):
    f=open(in_file,'r')

    file_lines=f.readlines()
    nr_files=0
    start=[]
    end=[]

    for i,line in enumerate(file_lines):
        if 'Welcome to Q-Chem' in line:
            start.append(i)
           
            nr_files=nr_files+1

        if 'Thank you very much for using Q-Chem.  Have a nice day.' in line:
           
            end.append(i)

        
    dict_files={}

    for j in range(0,len(start)-1):
        dict_files.update({j:file_lines[start[j]:end[j]]})

    return dict_files


'''
searches list with states to find the sought after amplitudes
'''
def search_amplitudes(amplitude,states):

    return_state=None
    for s in states:

        if s.amplitudes.iloc[0,1:3].to_list()==amplitude[0] and s.amplitudes.iloc[2,1:3].to_list()==amplitude[1]:
            return_state=s
            return return_state
        elif s.amplitudes.iloc[0,1:3].to_list()==amplitude[1] and s.amplitudes.iloc[2,1:3].to_list()==amplitude[0]:
            return_state=s
            return return_state


'''
makes a dataframe with the energies and angles 
'''        
def make_dataframe(states):

    re_e=[]
    im_e=[]
    angle=[]

    for s in states:
        re_e.append(s.get_re())
        im_e.append(s.get_im())
        angle.append(s.get_angle())


    data=[re_e,im_e,angle]
    labels = ['Re(E) (eV)','Im(E) (eV)','Angle (deg.)']
    df=pd.DataFrame(data,index=labels).transpose()
    print(df)


'''
search through qchem file for states that match the one with the character of the amplitudes that are given
'''
def states_per_file(dict_files,i,amplitudes):

    pg=find_pg(dict_files[i]) # find point group
    states_list=find_states(dict_files[i]) 
    states_per_sym=create_states_per_sym(pg,states_list)

    states=extract_states_file(dict_files[i],states_per_sym)
    sought_states=search_amplitudes(amplitudes,states)
   
    return sought_states

in_file='co_ee_cleaned.out'
dict_1=divide_to_files(in_file)


amplitude_1=[' 1 (B1) A','6 (A1) A']
amplitude_2=[' 1 (B1) A','8 (A1) A']
amplitudes=[amplitude_1,amplitude_2]
search_states=[]

search_state=[]
for i in range(0,len(dict_1.keys())):
    
    search_state=states_per_file(dict_1,i,amplitudes)
    #search_state.print_state()
    if search_state!=None:
        search_states.append(search_state)



make_dataframe(search_states)

