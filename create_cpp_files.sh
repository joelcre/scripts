#!/bin/bash

####################################################################
#Script that creates input files and run files for cpp-calculations
####################################################################

input_temp=$1 #name of input template
end=6.5
start=2.0
increment=0.0225 

python3 create_intervals.py $start $end $increment  #python script to create the smaller energy intervals 

file="intervals.txt"  #python script prints to temporary file
lines=`cat $file`

i=0  #counter
for line in $lines; do   #loop copies the template file and inserts energy values copied from file
    if [ $i -gt "0" ]; then
	cp ${input_temp}.inp  ${input_temp}_$i.inp    
	sed -i 's/end/'$line'/g' ${input_temp}_$i.inp
	sed -i 's/start/'$temp_line'/g' ${input_temp}_$i.inp
    fi

   temp_line=$line
   i=$((i+1)) 
   
done

nr_files=$((i-1))
echo $nr_files

# making run files

for j in $(seq 1 $nr_files)
 do  
	cp runtemp.sh  run_$j.sh    
	sed -i 's/index/'$j'/g' run_$j.sh   
	mkdir interval_$j
	mv run_$j.sh  ${input_temp}_$j.inp interval_$j
	cp  $2  interval_$j
	if [ -f $3 ]; then #apparently does not work, but makes no difference
	    echo $3
	    cp $3 interval_$j
	fi
done


