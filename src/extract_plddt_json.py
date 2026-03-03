import sys
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import json
import argparse
import os

def find_plddt(cif_file):
    """
    Extracts pLDDT of CA atoms.
    """
    plddts = []
    mmcif_dict = MMCIF2Dict(cif_file)

    b_factor_key = "_atom_site.B_iso_or_equiv"

    # Iterate through atoms to find CA atoms and their B-factors
    for i, atom_label in enumerate(mmcif_dict["_atom_site.label_atom_id"]):
        if atom_label == "CA":
            b_factor = float(mmcif_dict[b_factor_key][i])
            plddts.append(b_factor)

    return plddts

def main():
    parser = argparse.ArgumentParser(description="Extract CA pLDDT scores for GPCR-Gprotein pairs from MMCIF files.")
    parser.add_argument("--gpcr", required=True, help="Identifier for the GPCR.")
    parser.add_argument("--gprotein_list", default="gproteins.txt", help="Path to the file containing Gprotein names.")
    parser.add_argument("--input_dir", required=True, help="Input directory")
    parser.add_argument("--output_json", required=True,help="Path to the output JSON file where results will be saved.")

    args = parser.parse_args()

    gpcr = args.gpcr
    gpcr_lower = gpcr.lower()

    gproteins_df = pd.read_csv(args.gprotein_list, header=None)
    gproteins = gproteins_df[0].tolist()

    # Dictionary to store the results
    combinations = {}

    for gprotein in gproteins:
        gprotein_lower = gprotein.lower()
        cif_file_path = f"{args.input_dir}/{gpcr_lower}_{gprotein_lower}/{gpcr_lower}_{gprotein_lower}_model.cif"
        plddt_scores = find_plddt(cif_file_path)
        combination_key = f"{gpcr}_{gprotein}"
        combinations[combination_key] = plddt_scores

    # Save the results to a JSON file
    with open(args.output_json, 'w') as f_out:
        json.dump(combinations, f_out, indent=4) # Use indent for readability
    print(f"Successfully saved results to {args.output_json}")

if __name__ == "__main__":
    main()