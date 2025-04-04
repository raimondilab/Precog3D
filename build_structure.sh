#!/bin/bash
#PBS -l select=1:ncpus=16:mem=64gb:ngpus=1
#PBS -l walltime=48:00:00
#PBS -q q02gaia
#PBS -N af3

cd $HOME_FOLDER
singularity exec -H $HOME_FOLDER --nv --bind $DATABASE_FOLDER:/root/public_databases --bind $INPUT_FOLDER:/root/af_input --bind $OUTPUT_FOLDER:/root/af_output --bind $HOME_FOLDER:/root/models alphafold3.sif python run_alphafold.py --db_dir=/root/public_databases --json_path=/root/af_input/$INPUT_FILE --output_dir=/root/af_output --model_dir=/root/models