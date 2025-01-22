import numpy as np
import pandas as pd

class state:
    def __init__(self,re_ener_,im_ener_,amplitudes_,irrep_):
        self.re_ener=re_ener_
        self.im_ener=im_ener_
        self.amplitudes=amplitudes_
        self.irrep=irrep_

    def print_state(self):
        print('Irrep: ', self.irrep)
        print('Im(E): ', self.re_ener)
        print('Re(E): ', self.im_ener)
        print('Amplitudes: ',self.amplitudes)
        

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
    #print(irreps)
    for i in range(0,len(states_list)):
        states_per_sym.update({irreps[i]:states_list[i]})

    return states_per_sym



def find_state(qchem_file,state_str):

    for line in qchem_file:
        if state_str in line:
            print(line)


'''
Search file for states within an irrep and return a list with state objects containing
information about energy, amplitudes and irrep
'''            
def extract_states_file(qchem_file,states_per_sym):

    states=[]
    for irrep in states_per_sym.keys():
        if states_per_sym.get(irrep)!=0:
            nr_states=states_per_sym.get(irrep)

            for i in range(1,nr_states+1):
                search_str='CS/CAP-EOMEE-CCSD transition '+ str(i)+'/'+irrep
                #print(search_str)
                state=get_energy(search_str,qchem_file,irrep)
                states.append(state)

    #[state.print_state() for state in states]

    return states

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

    amplitudes_unform=qchem_file[init+1:end-1]
    amplitudes_form=[]
    for a in amplitudes_unform:
        amplitudes_form.append(strip_blankspaces(a.split('  ')))

    df=pd.DataFrame(amplitudes_form)
    #print(df)
    df2 = df[[0,1,3]]
    #print(df2)
    
    return df2

'''
'''
def get_energy(state_str,qfile,irrep):

    for i in range(0,len(qfile)):
        if state_str in qfile[i]:
            temp=[]
            #print(qfile[i+1].split(' ')[13].strip('('))
            #print(qfile[i+1].split(' ')[15].strip(')'))
            re_e=qfile[i+1].split(' ')[13].strip('(')
            im_e=qfile[i+1].split(' ')[15].strip(')')
            line_nr=i

    amplitudes=extract_amplitude(qfile,line_nr)
    state_=state(re_e,im_e,amplitudes,irrep)
    #state_.print()
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
    


def divide_to_files(in_file):
    f=open(in_file,'r')

    file_lines=f.readlines()
    nr_files=0
    start=[]
    end=[]

    for i,line in enumerate(file_lines):
        if 'Welcome to Q-Chem' in line:
            #print(i)
            start.append(i)
            #print(line)
            nr_files=nr_files+1

        if 'Thank you very much for using Q-Chem.  Have a nice day.' in line:
            #print(i)
            end.append(i)
            #print(line)
        
    dict_1={}

    for j in range(0,len(start)-1):
        dict_1.update({j:file_lines[start[j]:end[j]]})

    return dict_1


def search_amplitudes(amplitude,states):
    for s in states:
        #print(s.amplitudes.iloc[0,2])
        #print(amplitude[1])
        if s.amplitudes.iloc[0,1]==amplitude[0] and s.amplitudes.iloc[0,2]==amplitude[1]:
            s.print_state()

            


def states_per_file(dict_1,i):
    #print(dict_1[0])
    pg=find_pg(dict_1[i])
    #print(pg)
    states_list=find_states(dict_1[i])
    #print(states_list)
    states_per_sym=create_states_per_sym(pg,states_list)
    #print(states_per_sym)
    states=extract_states_file(dict_1[i],states_per_sym)

    amplitude=[' 1 (B3u) A','6 (B3u) A']
    search_amplitudes(amplitude,states)

    #print(type(states[0].amplitudes.iloc[0,1]))
    #print(states[0].amplitudes.iloc[0,2])
    #[state.print_state() for state in states]

in_file='ne_ee.out'
dict_1=divide_to_files(in_file)
print(len(dict_1))

for i in range(0,30):
    states_per_file(dict_1,i)
