#!/bin/bash

usage() {
    echo "Usage: $0 -h HOME_FOLDER -d DATABASE_FOLDER -i INPUT -o OUTPUT_FOLDER -a ALIGNMENT_FOLDER -g GPROTEIN_LIST" 
    echo "HOME_FOLDER: Home folder for alphafold3 singularity container."
    echo "DATABASE_FOLDER: Folder containing the databases."
    echo "INPUT: Input file in FASTA format."
    echo "OUTPUT_FOLDER: Folder where the output files will be stored."
    echo "ALIGNMENT_FOLDER: Folder where the alignments of the G-proteins are stored."
    exit 1
}

while getopts ":h:d:i:o:a:" opt; do
    case ${opt} in
        h )
            HOME_FOLDER=$OPTARG
            ;;
        d )
            DATABASE_FOLDER=$OPTARG
            ;;
        i )
            INPUT=$OPTARG
            ;;
        o )
            OUTPUT_FOLDER=$OPTARG
            ;;
        a )
            ALIGNMENT_FOLDER=$OPTARG
            ;;
        \? )
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Split the input file into file and folder
INPUT_FILE=$(basename $INPUT)
INPUT_FOLDER=$(dirname $INPUT)

python3 src/fasta_to_json_AF3.py \
    --fasta_file $INPUT \
    --job_name ${INPUT_FILE%%.fasta} \
    --output_file ${INPUT%%.fasta}.json

INPUT_NAME=${INPUT_FILE%%.fasta}

qsub -Wblock=true -v HOME_FOLDER=$HOME_FOLDER,DATABASE_FOLDER=$DATABASE_FOLDER,INPUT_FOLDER=$INPUT_FOLDER,INPUT_FILE=$INPUT_NAME.json,OUTPUT_FOLDER=$OUTPUT_FOLDER src/build_alignment.sh

gprotein_list=("P29992" "Q03113" "Q14344" "O95837" "P30679" "P63096" "P04899" "P08754" "P38405" "P50148" "P63092" "P19086" "P09471")

for gprotein in ${gprotein_list[@]}; do
    python3 src/generate_json.py \
        ${INPUT%%.fasta} \
        $ALIGNMENT_FOLDER/$gprotein \
        $ALIGNMENT_FOLDER/P62873 \
        $ALIGNMENT_FOLDER/P59768 \
        --output_dir $OUTPUT_FOLDER \
        --name $INPUT_NAME\_$gprotein
    qsub -Wblock=true \
        -v HOME_FOLDER=$HOME_FOLDER,DATABASE_FOLDER=$DATABASE_FOLDER,INPUT_FOLDER=$INPUT_FOLDER,INPUT_FILE=$INPUT_NAME\_$gprotein.json,OUTPUT_FOLDER=$OUTPUT_FOLDER \
        src/build_structure.sh
    python3 src/filter.py \
        --input_cif $OUTPUT_FOLDER/${INPUT_NAME,,}_${gprotein,,}/${INPUT_NAME,,}_${gprotein,,}_model.cif \
        --output_pdb $OUTPUT_FOLDER/${INPUT_NAME,,}_${gprotein,,}/${INPUT_NAME,,}_${gprotein,,}_filtered.pdb
    python3 src/pdockq.py \
        --pdbfile $OUTPUT_FOLDER/${INPUT_NAME,,}_${gprotein,,}/${INPUT_NAME,,}_${gprotein,,}_filtered.pdb \
        > $OUTPUT_FOLDER/${INPUT_NAME,,}_${gprotein,,}/${INPUT_NAME,,}_${gprotein,,}_pdockq.txt
done

qsub -Wblock=true -v INPUT_FOLDER=$OUTPUT_FOLDER,INPUT_NAME=${INPUT_NAME,,},OUTPUT_FOLDER=$OUTPUT_FOLDER src/relax.sh
qsub -Wblock=true -v INPUT_FOLDER=$OUTPUT_FOLDER,INPUT_NAME=${INPUT_NAME,,},OUTPUT_FOLDER=$OUTPUT_FOLDER src/interface_energy.sh
blastp -query $INPUT -db GPCRDB -out $OUTPUT_FOLDER/$INPUT_NAME.out -outfmt 0
python3 src/extract_plddt_json.py \
    --gpcr $INPUT_NAME \
    --input_dir $OUTPUT_FOLDER \
    --output_json $OUTPUT_FOLDER/${INPUT_NAME,,}_plddt.json
python3 src/create_table.py \
    --gpcr-id $INPUT_NAME \
    --plddt-json $OUTPUT_FOLDER/${INPUT_NAME,,}_plddt.json \
    --blastp-out $OUTPUT_FOLDER/$INPUT_NAME.out \
    --results-dir $OUTPUT_FOLDER \
    --output-file $OUTPUT_FOLDER/$INPUT_NAME\_table.tsv
qsub -Wblock=true -v OUTPUT_FOLDER=$OUTPUT_FOLDER,INPUT_NAME=$INPUT_NAME,HOME_FOLDER=$PWD src/model.sh