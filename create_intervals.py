#########################################################################################
#Script that cuts an energy interval into smaller pieces, and puts them into a text-file#
#########################################################################################

import math
import sys


start=float(sys.argv[1])
end= float(sys.argv[2])
increment=float(sys.argv[3])

ev=27.2114
start=start/ev
end=end/ev     
val=start


nr_files=round((end-start)/increment)
inter='intervals.txt'
i= open(inter,"w+")

for j in range(0,nr_files+1):
    val= float("{0:.4f}".format(val))
    i.write(str(val))
    i.write('\n')
    val=val+increment



   
