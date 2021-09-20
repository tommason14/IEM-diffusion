#!/usr/bin/env bash

# balance charge of polymer
gmx_mpi insert-molecules -f polymer.gro -ci sodium.gro -nmol 24 -o neutral.gro

# position in box - CC polymer isn't uniform, so this is trial and error
gmx_mpi editconf -f neutral.gro -c -center 2.5 2.5 2 -box 5 5 4 -o edge.gro

# stack
gmx_mpi genconf -f edge.gro -nbox 1 1 5 -o stacked.gro
echo "Edit topol.top to include 5 polymer and associated sodiums, 1 by 1. So 1 polymer, 24 sodiums, 1 polymer 24 sodiums.
Desired output:
[ molecules ]
pol                 1
NA                 24
pol                 1
NA                 24
...
"
read -p "Press any key to open the file..."
vim topol.top

gmx_mpi solvate -cp stacked.gro -cs tip4p.gro -p topol.top -o stacked_solv.gro
# now add 0.2 M ions 
gmx_mpi grompp -f ions.mdp -c stacked_solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx_mpi genion -s ions.tpr -p topol.top -o pack.gro -pname NA -nname CL -conc 0.2

find . -maxdepth 1 -name "#*" -delete
