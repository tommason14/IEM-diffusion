#!/bin/sh
# MDP='em.mdp' # any options can be included - just need to create a tpr for gmx tools
GRO='../../../pack/pack.gro'
TOP='../../../pack/topol-scaled.top'
XTC='../traj.xtc'
UXTC='unwrapped.xtc'
TPR='../../../pack/run-scaled.tpr'
module load gromacs/2020.1
# curl -O https://raw.githubusercontent.com/tommason14/IEM-diffusion/main/setup/pack/em.mdp
# Remember to check Na and Cl atom names otherwise traj_convert will fail!

# gmx grompp -f $MDP -c $GRO -p $TOP -o $TPR
echo System | gmx trjconv -s $TPR -f $XTC -pbc nojump -o $UXTC
echo 2 | gmx dipoles -f $UXTC -s $TPR > dipoles.log
dielectric="$(tail -1 dipoles.log | getcol -1)"
gmx rdf -f $UXTC -s $TPR -ref 'name "S*"' -sel 'name "S*"' -o s-s-rdf.xvg
sep="$(grep -v '^[#@]' s-s-rdf.xvg | sort -k2n | tail -1 | awk '{print $1*10}')"
# gmx rdf -f $UXTC -s $TPR -ref 'name "N[0-9]*"' -sel 'name "N[0-9]*"' -o n-n-rdf.xvg
# sep="$(grep -v '^[#@]' n-n-rdf.xvg | sort -k2n | tail -1 | awk '{print $1*10}')"
calc_bjerrum_length.py -d $dielectric -s $sep > bjerrum.txt
bj="$(getcol -2 bjerrum.txt)"

traj_convert.py -c $GRO -t $UXTC --sel "name Cl" -o chlorides.pdb
traj_convert.py -c $GRO -t $UXTC --sel "name Na" -o sodiums.pdb

curl -O https://raw.githubusercontent.com/tommason14/IEM-diffusion/main/analysis/cl_input.txt
curl -O https://raw.githubusercontent.com/tommason14/IEM-diffusion/main/analysis/na_input.txt

# trajectories saved every 20 ps, not 10 ps as with a 1 fs timestep
sed -i 's/10000/20000/' cl_input.txt
sed -i 's/10000/20000/' na_input.txt

travis -p sodiums.pdb -i na_input.txt && mv travis.log sodium.log
travis -p chlorides.pdb -i cl_input.txt && mv travis.log chloride.log

boxsize="$(sed -n '/Step/,//p' ../nvt.log | csvcut -c "Box Volume (nm^3)" | tail -1 | awk '{print ($1**(1/3))*10}')"

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

# AR103 values
# exp_na="9e-7"
# exp_cl="1.8e-6"

perc_diff_na=$(awk -v ref="$exp_na" -v calc="$na_diff" 'BEGIN{print (calc-ref)/ref * 100}')
perc_diff_na_corr=$(awk -v ref="$exp_na" -v calc="$na_corrected" 'BEGIN{print (calc-ref)/ref * 100}')
perc_diff_cl=$(awk -v ref="$exp_cl" -v calc="$cl_diff" 'BEGIN{print (calc-ref)/ref * 100}')
perc_diff_cl_corr=$(awk -v ref="$exp_cl" -v calc="$cl_corrected" 'BEGIN{print (calc-ref)/ref * 100}')

cat << EOF | column -t -s,
Dielectric,S-S dist (Å),Bjerrum length (Å),Boxsize (Å),D_Na (cm²/s),D_Cl (cm²/s),% diff Na,% diff Cl,Corrected D_Na (cm²/s),Corrected D_Cl (cm²/s),% diff Na,% diff Cl
$dielectric,$sep,$bj,$boxsize,$na_diff,$cl_diff,$perc_diff_na,$perc_diff_cl,$na_corrected,$cl_corrected,$perc_diff_na_corr,$perc_diff_cl_corr
EOF
