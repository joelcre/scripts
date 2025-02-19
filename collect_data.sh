#!/bin/bash

#########################################################################################
# script that goes thorugh a series of structure catalogues with data and transforms it
# into files that can be plotted as spectra
###########################################################################################

Structures=('opt_103' 'opt_109' 'opt_113' 'opt_119' 'opt_121' 'opt_125' 'opt_129' 'opt_131' 'opt_133' 'opt_135' 'opt_145' 'opt_157' 'opt_167' 'opt_169' 'opt_173' 'opt_177' 'opt_179' 'opt_183' 'opt_187')



functionals=("CAMB3LYP")
Hamiltonians=("x2c")

for structure in ${Structures[*]};
do
    echo $structure
    
    cd $structure
    rm construct_cpp_spectra.py
    cp ../construct_cpp_spectra.py .
    for hamiltonian in ${Hamiltonians[*]};
    do
	
	cd ${hamiltonian}
	rm construct_cpp_spectra.py
	cp ../construct_cpp_spectra.py .
	for functional in ${functionals[*]};
	do

	    cd ${functional}
	    rm construct_cpp_spectra.py
	    cp ../construct_cpp_spectra.py .
	    #mkdir output
	    #cp interval_*/*.out output
	    cd output
	    #rm slurm*
	    
	    if grep -q 'NOT ' *.out; then
		echo "not converged"
	       
	    else
		statement="$(cat *.out | grep -c 'E N D   of   D I R A C  output')"

		if [ "$statement" -eq 7 ]; then
		    #if grep -q 'E N D   of   D I R A C  output' *.out; then
		    rm construct_cpp_spectra.py
		    cp ../construct_cpp_spectra.py .
		    str_1='_'
		    str_2='.out'
		    str_3='spect_'
		    str_4='.txt'
		    file_suffix=$str_1$structure$str_1$structure$str_2
		    spect_file=$str_3$structure$str_4
		    sed -i 's/_opt_13_opt_13.out/'$file_suffix'/g' construct_cpp_spectra.py
	  	    sed -i 's/spect_13.txt/'$spect_file'/g' construct_cpp_spectra.py

		    python3 construct_cpp_spectra.py
		    cp $spect_file ../../../../
		    cd ..
		fi
	    fi
	    
	    cd EEF
	    
	    #mkdir output
	    #cp interval_*/*.out output
	    cd output
	    #rm slurm*
	    
		if grep -q 'NOT ' *.out; then

		    echo 'not converged'
	        else
		    statement="$(cat *.out | grep -c 'E N D   of   D I R A C  output')"

		    if [ "$statement" -eq 7 ]; then
			#if grep -q 'E N D   of   D I R A C  output' *.out; then
			rm construct_cpp_spectra.py
			cp ../../construct_cpp_spectra.py .
			str_1='_'
			str_2='.out'
			str_3='spect_'
			str_4='_EEF.txt'
			file_suffix=$str_1$structure$str_1$structure$str_2
			spect_file=$str_3$structure$str_4
			sed -i 's/_opt_13_opt_13.out/'$file_suffix'/g' construct_cpp_spectra.py
			sed -i 's/spect_13.txt/'$spect_file'/g' construct_cpp_spectra.py
			python3 construct_cpp_spectra.py
			cp $spect_file ../../../../../
			cd ..
		    fi
		fi
		
		cd ..
		

	    cd ..
	    done
	    
	cd ..
    done
    cd ..
    
    done
		     
