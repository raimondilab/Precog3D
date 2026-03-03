#!/bin/bash

# Ensure the script stops if a command fails
set -e

# Check if the correct number of arguments is provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 uniprot_accessions.txt"
    echo "uniprot_accessions.txt: File containing UniProt accessions (one per line)."
    exit 1
fi

# Input parameters
UNIPROT_LIST=$1

# Split the input UniProt accessions file into chunks of 20 accessions per file
echo "Splitting input file into smaller files, each containing 20 accessions..."
split -l 20 "$UNIPROT_LIST" "uniprot_chunk_"

# Process each chunk and generate a single PBS script
for CHUNK in uniprot_chunk_*; do
    echo "Processing chunk: $CHUNK"

    # Create PBS script for this chunk
    PBS_SCRIPT="pbs_${CHUNK}_job.sh"
    echo "#!/bin/bash" > $PBS_SCRIPT
    echo "#PBS -l select=1:ncpus=16:mem=32gb" >> $PBS_SCRIPT
    echo "#PBS -l walltime=48:00:00" >> $PBS_SCRIPT
    echo "#PBS -q q02anacreon" >> $PBS_SCRIPT
    echo "#PBS -N ${CHUNK}_alignment" >> $PBS_SCRIPT
    echo "" >> $PBS_SCRIPT
    echo "cd /projects/bioinformatics/SW/alphafold3" >> $PBS_SCRIPT

    # Download FASTA sequences, generate JSONs, and add Singularity runs to PBS script
    while read -r ACCESSION; do
        if [[ -n "$ACCESSION" ]]; then
            # Download FASTA
            wget "https://rest.uniprot.org/uniprotkb/${ACCESSION}.fasta"
            echo "Downloaded sequence for UniProt accession: $ACCESSION"

            # Generate JSON file using Python script
            python3 fasta_to_json_AF3.py \
                --fasta_file ${ACCESSION}.fasta \
                --job_name ${ACCESSION} \
                --output_file ${ACCESSION}.json

            # Add Singularity run command to PBS script for this accession
            echo "singularity exec -H /projects/bioinformatics/ --bind /projects/bioinformatics/DB/AF3_DB/:/root/public_databases --bind /home/pmiglionico/GPCR_AF3/:/root/af_input --bind /home/pmiglionico/GPCR_AF3:/root/af_output alphafold3.sif python run_alphafold_short.py --db_dir=/root/public_databases --json_path=/root/af_input/${ACCESSION}.json --output_dir=/root/af_output --norun_inference" >> $PBS_SCRIPT
            echo "rm ${ACCESSION}.json" >> $PBS_SCRIPT

            # Increment counter
            #rm ${ACCESSION}.fasta
        fi
    done < "$CHUNK"

    # Cleanup FASTA files for the current chunk

    echo "Generated PBS script for ${CHUNK}: ${PBS_SCRIPT}"
    qsub $PBS_SCRIPT
done

# Cleanup split chunks
rm uniprot_chunk_*
echo "All PBS scripts generated and temporary files cleaned up."