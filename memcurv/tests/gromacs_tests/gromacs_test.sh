#!/bin/bash

# generate a system and try to grompp it with topology and index 
folder="gromacs_test_files"

rm $folder/*.tpr $folder/*.gro $folder/test*

python3 ../../main.py -f asymm.pdb -s sphere -z 10.5 -o $folder/test_start.pdb -p $folder/test_start.top -n $folder/test_start.ndx -g r_sphere:100

echo '#include "toppar/martini_v2.0_ions.itp"'              | cat - $folder/test_start.top > temp && mv temp $folder/test_start.top
echo '#include "toppar/martini_v2.0_lipids_all_201506.itp"' | cat - $folder/test_start.top > temp && mv temp $folder/test_start.top
echo '#include "toppar/martini_v2.2.itp"'                   | cat - $folder/test_start.top > temp && mv temp $folder/test_start.top


gmx grompp -f $folder/step6.0_equilibration.mdp -p $folder/test_start.top -n $folder/test_start.ndx -c $folder/test_start.pdb -o $folder/test_minimization.tpr 
gmx mdrun -deffnm $folder/test_minimization -nsteps 40 -v
gmx grompp -f $folder/step7_production.mdp -p $folder/test_start.top -n $folder/test_start.ndx -c $folder/test_minimization.gro -o $folder/test_grompp_production.tpr



if [ $? -eq 0 ]; then
	echo 'Success'
else
	echo 'Failure :('
fi

rm $folder/#*# mdout.mdp ./#*#

exit 
