import json
import string
import argparse

def parse_fasta(fasta_file):
    sequences = []
    header = None
    sequence = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        
        # Add the last record
        if header:
            sequences.append((header, "".join(sequence)))
    
    return sequences

def parse_fasta_to_json(
    fasta_file, job_name="AlphaFoldJob", model_seeds=None, ligands=None, output_file="alphafold_input.json"
):
    """
    Converts a FASTA file with protein sequences and ligands into the AlphaFold 3 JSON format.
    """
    if model_seeds is None:
        model_seeds = [42]
    if ligands is None:
        ligands = []

    # Initialize chain IDs (progressive letters: A, B, C, ...)
    chain_ids = iter(string.ascii_uppercase)

    # Parse FASTA file
    parsed_sequences = parse_fasta(fasta_file)
    sequences = []
    for _, seq in parsed_sequences:
        chain_id = next(chain_ids)  # Get the next chain ID
        sequences.append({
            "protein": {
                "id": chain_id,
                "sequence": seq
            }
        })

    # Add ligands to sequences
    for ligand in ligands:
        sequences.append({
            "ligand": ligand
        })

    # Build the JSON structure
    alphafold_input = {
        "name": job_name,
        "modelSeeds": model_seeds,
        "sequences": sequences,
        "bondedAtomPairs": [],
        "dialect": "alphafold3",
        "version": 1 
    }

    # Save to file
    with open(output_file, "w") as outfile:
        json.dump(alphafold_input, outfile, indent=4)
    print(f"AlphaFold input JSON saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Convert a FASTA file to AlphaFold JSON input format.")
    parser.add_argument("--fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("--job_name", type=str, help="Name of the AlphaFold job.")
    parser.add_argument("--output_file", type=str, default="alphafold_input.json", help="Path to save the output JSON.")
    parser.add_argument("--model_seeds", type=int, nargs="+", default=[42], help="List of random seeds for AlphaFold.")
    parser.add_argument("--ligands", type=json.loads, default="[]", help="JSON-formatted string of ligand definitions.")
    
    args = parser.parse_args()

    # Convert ligands argument from JSON string to list if provided
    ligands = json.loads(args.ligands) if args.ligands else []

    # Run the conversion
    parse_fasta_to_json(
        fasta_file=args.fasta_file,
        job_name=args.job_name,
        model_seeds=args.model_seeds,
        ligands=ligands,
        output_file=args.output_file
    )

if __name__ == "__main__":
    main()