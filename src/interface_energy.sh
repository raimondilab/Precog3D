#!/bin/bash

#PBS -l select=1:ncpus=12:mem=40gb
#PBS -l walltime=48:00:00
#PBS -q q02anacreon
#PBS -N interface_energy

module load gcc
cd $INPUT_FOLDER
echo $INPUT_NAME
gprotein_list=("P29992" "Q03113" "Q14344" "O95837" "P30679" "P63096" "P04899" "P08754" "P38405" "P50148" "P63092" "P19086" "P09471")

for gprotein in ${gprotein_list[@]}; do
    /projects/bioinformatics/SW/rosetta.source.release-371/main/source/bin/InterfaceAnalyzer.linuxgccrelease -s $INPUT_NAME\_${gprotein,,}/$INPUT_NAME\_${gprotein,,}_filtered_0001.pdb -interface A_B -out:file:score_only $OUTPUT_FOLDER/$INPUT_NAME\_${gprotein,,}/$INPUT_NAME\_${gprotein,,}_interface.tsv -pack_input -pack_separated -ignore_unrecognized_res &
done
wait
