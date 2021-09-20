#!/bin/sh
EM='../../../pack/em.mdp'
GRO='../../npt/equil.gro'
TOP='../../../pack/topol.top'
XTC='../output.xtc'
UXTC='unwrapped.xtc'
TPR='em.tpr'

gmx_mpi grompp -f $EM -c $GRO -p $TOP -o $TPR
echo System | gmx_mpi trjconv -s $TPR -f $XTC -pbc nojump -o $UXTC
echo 2 | gmx_mpi dipoles -f $UXTC -s $EM > dipoles.log
gmx_mpi rdf -f $UXTC -s $EM -ref 'name "S*"' -sel 'name "S*"' -o s-s-rdf.xvg
dielectric="$(tail -1 dipoles.log | getcol -1)"
sep="$(grep -v '^[#@]' s-s-rdf.xvg | sort -k2n | tail -1 | awk '{print $1*10}')"
calc_bjerrum_length.py -d $dielectric -s $sep > bjerrum.txt
bj="$(getcol -2 bjerrum.txt)"

traj_convert.py -c $GRO -t $UXTC --sel "name CL and around $bj name S*" -o chlorides.pdb
traj_convert.py -c $GRO -t $UXTC --sel "name NA and around $bj name S*" -o sodiums.pdb

travis -p sodiums.pdb -i $thisdir/na_input.txt && mv travis.log sodium.log
travis -p chlorides.pdb -i $thisdir/cl_input.txt && mv travis.log chloride.log

# cube root the product of box dimensions - accounts for non-standard boxes
boxsize="$(tail -1 $gro | awk '{print ($1*$2*$3)**(1/3)*10}')"

# convert to cm^2/s
na_msd=$(grep "m^2/s" sodium.log | awk '{print $(NF-1)*10**4}')
cl_msd=$(grep "m^2/s" chloride.log | awk '{print $(NF-1)*10**4}')

echo "# Sodium" >> corrected_diff_coeffs.txt
corrected_diff_coeff.py -b $boxsize -d $na_msd >> corrected_diff_coeffs.txt
echo "# Chloride" >> corrected_diff_coeffs.txt
corrected_diff_coeff.py -b $boxsize -d $cl_msd >> corrected_diff_coeffs.txt
