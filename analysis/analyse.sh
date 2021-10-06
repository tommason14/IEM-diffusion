#!/bin/sh
MDP='../../../pack/em.mdp' # any options can be included - just need to create a tpr for gmx tools
GRO='../../npt/equil.gro'
TOP='../../../pack/topol.top'
XTC='../output.xtc'
UXTC='unwrapped.xtc'
TPR='em.tpr'

# Remember to check Na and Cl atom names otherwise traj_convert will fail!

gmx_mpi grompp -f $MDP -c $GRO -p $TOP -o $TPR
echo System | gmx_mpi trjconv -s $TPR -f $XTC -pbc nojump -o $UXTC
echo 2 | gmx_mpi dipoles -f $UXTC -s $TPR > dipoles.log
gmx_mpi rdf -f $UXTC -s $TPR -ref 'name "S*"' -sel 'name "S*"' -o s-s-rdf.xvg
dielectric="$(tail -1 dipoles.log | getcol -1)"
sep="$(grep -v '^[#@]' s-s-rdf.xvg | sort -k2n | tail -1 | awk '{print $1*10}')"
calc_bjerrum_length.py -d $dielectric -s $sep > bjerrum.txt
bj="$(getcol -2 bjerrum.txt)"

traj_convert.py -c $GRO -t $UXTC --sel "name CL and around $bj name S*" -o chlorides.pdb
traj_convert.py -c $GRO -t $UXTC --sel "name NA and around $bj name S*" -o sodiums.pdb

curl -O https://raw.githubusercontent.com/tommason14/IEM-diffusion/main/analysis/cl_input.txt
curl -O https://raw.githubusercontent.com/tommason14/IEM-diffusion/main/analysis/na_input.txt

travis -p sodiums.pdb -i na_input.txt && mv travis.log sodium.log
travis -p chlorides.pdb -i cl_input.txt && mv travis.log chloride.log

# cube root the product of box dimensions - accounts for non-standard boxes
boxsize="$(tail -1 $GRO | awk '{print ($1*$2*$3)**(1/3)*10}')"

# convert to cm^2/s
na_diff=$(grep "m^2/s" sodium.log | awk '{print $(NF-1)*10**4}')
cl_diff=$(grep "m^2/s" chloride.log | awk '{print $(NF-1)*10**4}')

echo "# Sodium" >> corrected_diff_coeffs.txt
corrected_diff_coeff.py -b $boxsize -d $na_diff >> corrected_diff_coeffs.txt
echo "# Chloride" >> corrected_diff_coeffs.txt
corrected_diff_coeff.py -b $boxsize -d $cl_diff >> corrected_diff_coeffs.txt

### Report results ###
# sodium ran before chloride
na_corrected=$(grep "PBCs" corrected_diff_coeffs.txt | head -1 | getcol -2)
cl_corrected=$(grep "PBCs" corrected_diff_coeffs.txt | tail -1 | getcol -2)

# CR61 values
exp_na="2e-6"
exp_cl="1.8e-6"

perc_diff_na=$(awk -v ref="$exp_na" -v calc="$na_diff" 'BEGIN{print (ref-calc)/ref * 100}')
perc_diff_na_corr=$(awk -v ref="$exp_na" -v calc="$na_corrected" 'BEGIN{print (ref-calc)/ref * 100}')
perc_diff_cl=$(awk -v ref="$exp_cl" -v calc="$cl_diff" 'BEGIN{print (ref-calc)/ref * 100}')
perc_diff_cl_corr=$(awk -v ref="$exp_cl" -v calc="$cl_corrected" 'BEGIN{print (ref-calc)/ref * 100}')

cat << EOF | column -t -s,
Dielectric,S-S dist (Å),Bjerrum length (Å),Boxsize (Å),D_Na (cm²/s),D_Cl (cm²/s),% diff Na,% diff Cl,Corrected D_Na (cm²/s),Corrected D_Cl (cm²/s),% diff Na,% diff Cl
$dielectric,$sep,$bj,$boxsize,$na_diff,$cl_diff,$perc_diff_na,$perc_diff_cl,$na_corrected,$cl_corrected,$perc_diff_na_corr,$perc_diff_cl_corr
EOF
