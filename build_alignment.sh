#!/bin/bash
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -l walltime=48:00:00
#PBS -q q02anacreon
#PBS -N af_alignment

cd $HOME_FOLDER
singularity exec -H $HOME_FOLDER --bind $DATABASE_FOLDER:/root/public_databases --bind $INPUT_FOLDER:/root/af_input --bind $OUTPUT_FOLDER:/root/af_output alphafold3.sif python run_alphafold.py --db_dir=/root/public_databases --json_path=/root/af_input/$INPUT_FILE --output_dir=/root/af_output --norun_inference