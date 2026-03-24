#!/bin/bash
#PBS -l select=1:ncpus=2:mem=32gb:ngpus=1
#PBS -l walltime=48:00:00
#PBS -q q02gaia
#PBS -N model

cd $HOME_FOLDER
source activate TabPFN
python3 src/model.py \
    --input_file $OUTPUT_FOLDER/$INPUT_NAME\_table.tsv \
    --output_file $OUTPUT_FOLDER/$INPUT_NAME\_predictions.tsv