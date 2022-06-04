####################################################################
# A program that does the following:
# 	1. run <gmx mdrun> for $nsteps steps
#	2. run <gmx energy> to calc. P-ij.
#	3. run <python multiple-tau.py> to calculate the correlator
# 	   for each component of the pressure tensor P-ij.
#	4. Copy the energy file with a different frequency 
#          with <gmx eneconv>
#	6. Remove the large original energy file.
#	7. repeat from step 1.
####################################################################
nsteps=1000000 # 10'000 steps = 2MB .edr file
echo Beginning run 01 Initialization run.
gmx mdrun -deffnm $1 -v -nsteps ${nsteps}
echo 'pres-' | gmx energy -f $1.edr
python3 multiple-tau.py energy.xvg
gmx eneconv -f $1.edr -dt 200 -o $1.part0001-conv.edr
rm energy.xvg
rm $1.edr


for i in {0002..0240}
do
	echo Beginning run $i.
	gmx mdrun -deffnm $1 -cpi $1.cpt -nsteps ${nsteps} -noappend -v
	echo 'pres-' | gmx energy -f $1.part$i.edr
	python3 multiple-tau.py energy.xvg
	rm energy.xvg ; rm $1.part$i.edr
	gmx eneconv -f $1.part$i.edr -dt 200 -o $1.part$i-conv.edr
        rm energy.xvg
        rm $1.part$i.edr
done

