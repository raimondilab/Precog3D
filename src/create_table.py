import pandas as pd
import numpy as np
import json
import pickle
from tqdm import tqdm
import argparse

def normalize(df, column):
    norm_col_name = f"{column}_norm"
    if norm_col_name not in df.columns:
        df[norm_col_name] = pd.NA
    # Group by GPCR and normalize
    group_means = df.groupby("GPCR")[column].transform('mean')
    df[norm_col_name] = df[column] - group_means
    return df

def read_score(GPCR, Gprotein, base_dir="."):
    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/ranking_scores.csv"
    ranking = pd.read_csv(filepath)
    return max(ranking["ranking_score"])

def read_iptm(GPCR, Gprotein, base_dir="."):
    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/{GPCR.lower()}_{Gprotein.lower()}_summary_confidences.json"
    with open(filepath, 'r') as f:
        iptm = json.load(f)
        return iptm["chain_pair_iptm"][0][1]

def parse_blastp_output(blastp_output_path):
    """
    Parses a BLASTP output file (format 0, from command line) and extracts
    the alignment correspondence between the query sequence and its best hit.
    """

    query_seq = []
    sbjct_seq = []
    uniac = ""

    with open(blastp_output_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if uniac:
                    break
                uniac = line.split()[1].split("|")[1]
            if line.startswith("Query") and uniac:
                parts = line.split()
                query_seq.append((int(parts[1]), parts[2]))
            elif line.startswith("Sbjct"):
                parts = line.split()
                sbjct_seq.append((int(parts[1]), parts[2]))

    # Build the correspondence dictionary
    correspondence_dict = {}
    for i in range(len(query_seq)):
        query_pos = query_seq[i][0]
        hit_pos = sbjct_seq[i][0]
        for j in range(len(query_seq[i][1])):
            if query_seq[i][1][j] == '-':
                hit_pos += 1
            elif sbjct_seq[i][1][j] == '-':
                query_pos += 1
            else:
                correspondence_dict[query_pos] = hit_pos
                query_pos += 1
                hit_pos += 1

    return correspondence_dict, uniac

def BW_mapping(mapping_dict, uniac, pos, BW_data):
    """Looks up BW mapping using the BLAST alignment."""
    if pos in mapping_dict:
        return BW_data[(uniac, mapping_dict[pos])][1]
    return "-"

def get_interaction_and_contact_probs(GPCR, Gprotein, base_dir="."):
    """
    Calculates both the interaction probability and contact probabilities between a GPCR and G-protein.
    """

    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/{GPCR.lower()}_{Gprotein.lower()}_confidences.json"
    with open(filepath, "r") as f:
        data = json.load(f)
    A_end = 0
    B_end = 0
    for i, chain in enumerate(data["token_chain_ids"]):
        if chain == "A":
            A_end = i
        if chain == "C":
            B_end = i - 1
            break

    # Calculate interaction probability
    interaction_prob_list = []
    for i in range(A_end + 1):
        interaction_prob_list.append(max(data["contact_probs"][i][A_end + 1 : B_end + 1]))
    interaction_prob = max(interaction_prob_list)

    # Calculate contact probabilities for the last 7 residues of the G-protein
    contact_probs = []
    for i in range(7):
        contact_probs.append(max(data["contact_probs"][B_end - i][: A_end + 1]))

    return [interaction_prob] + contact_probs

def find_interface_energy(GPCR, Gprotein, base_dir="."):
    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/{GPCR.lower()}_{Gprotein.lower()}_interface.tsv"
    with open(filepath) as f:
        lines = f.readlines()
        line = lines[2].strip().split()
        # Filter out empty strings resulting from multiple spaces
        line = [val for val in line if val]
        dG = float(line[5])
        dG_dSASA = float(line[6])
        dSASA = float(line[8])
        return [dG, dG_dSASA, dSASA]

def find_pdockq(GPCR, Gprotein, base_dir="."):
    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/{GPCR.lower()}_{Gprotein.lower()}_pdockq.txt"
    with open(filepath) as f:
        row = f.readline()
        parts = row.split()
        return float(parts[2])

def get_interaction_prob_specific(GPCR, Gprotein, Gprotein_name, pos1_BW, pos2_CGN, BW2pos_map, G2pos_map, base_dir=".", plddts_data=None, mode="contact"):
    """
    Gets specific contact probabilities or pLDDT values between GPCR and G-protein residues.
    """

    # Map BW/CGN positions to indices
    positions1_indices = [] # 0-based indices for GPCR
    for bw in pos1_BW:
        seq_pos_1based = BW2pos_map.get(bw)
        if seq_pos_1based is not None:
            positions1_indices.append(seq_pos_1based - 1) # Convert to 0-based
        else:
            print(f"Warning: BW position '{bw}' not found in BW2pos map for {GPCR}.")
            positions1_indices.append(None) # Mark as unmappable

    positions2_indices = [] # 1-based indices for G-protein
    for cgn in pos2_CGN:
        positions2_indices.append(G2pos_map[(Gprotein_name, cgn)])

    # Load data
    filepath = f"{base_dir}/{GPCR.lower()}_{Gprotein.lower()}/{GPCR.lower()}_{Gprotein.lower()}_confidences.json"
    with open(filepath, "r") as f:
        conf_data = json.load(f)
    # Find chain boundaries
    A_end = -1
    for i, chain in enumerate(conf_data["token_chain_ids"]):
        if chain == "A":
            A_end = i
        else:
            break

    # Calculate results based on mode
    results = []

    if mode == 'contact':
        for i in range(len(pos1_BW)):
            idx1 = positions1_indices[i]
            idx2 = positions2_indices[i]
            if idx1 is None or idx2 is None:
                results.append(np.nan)
                continue
            results.append(conf_data["contact_probs"][idx1][idx2+A_end])

    elif mode == 'plddt':
        pair_key = f"{GPCR}_{Gprotein}"
        plddt_values = plddts_data.get(pair_key)
        num_tokens = len(plddt_values)
        # Get pLDDT for GPCR positions
        gpcr_plddts = []
        for idx1 in positions1_indices:
            if idx1 is not None:
                gpcr_plddts.append(plddt_values[idx1])
            else:
                gpcr_plddts.append(np.nan)
        # Get pLDDT for G-protein positions
        gprot_plddts = []
        for idx2 in positions2_indices:
            if idx2 is not None:
                gprot_plddts.append(plddt_values[idx2 + A_end])
            else:
                gprot_plddts.append(np.nan)
        results = gpcr_plddts + gprot_plddts

    return results

def main(args):

    gpcr = args.gpcr_id
    print(f"Processing GPCR: {gpcr}")

    # Load G-protein list
    gproteins_df = pd.read_csv(args.gproteins_list, header=None)
    gprotein_list = gproteins_df[0].tolist()

    # Load BW mapping data
    with open(args.bw_pickle, "rb") as f:
        BW_data = pickle.load(f)

    # Load pLDDT data
    with open(args.plddt_json, 'r') as f:
        plddts = json.load(f)

    # Load possible GPCR-Gprotein contacts
    contacts = pd.read_csv(args.contacts_tsv, sep='\t')

    # Load conserved GPCR positions (BW format)
    with open(args.valid_gpcr_pos, "r") as f:
        valid_gpcr_positions_bw = [line.strip() for line in f if line.strip()]
    valid_gpcr_positions_bw.sort()

    # Load consensus G-protein positions (CGN format)
    with open(args.consensus_gprot_pos, "r") as f:
        consensus_gprotein_pos_cgn = [line.strip() for line in f if line.strip()]
    consensus_gprotein_pos_cgn.sort()

    # Load G-protein CGN to position mapping
    Gprot_positions_df = pd.read_csv(args.gprot_pos_map)

    # Create G2pos dictionary (Gprotein_name, CGN) -> 0-based structure index
    G2pos = {}
    for _, row in Gprot_positions_df.iterrows():
        key = (row["Protein"], row["CGN"])
        G2pos[key] = row["Position"] # Assuming 'Position' is the 0-based index needed

    # DataFrame Setup
    combinations = [[gpcr, gprot] for gprot in gprotein_list]
    pairs = pd.DataFrame(combinations, columns=['GPCR', 'Gprotein'])

    # Get G-protein Names from DB
    names = json.load(open("data/gprotein_names.json", "r"))

    pairs["Gprotein_name"] = pairs["Gprotein"].map(lambda x: names[x])

    # Add Basic Metrics
    print("Calculating ranking scores, iPTM, interaction/contact probs, interface energy, pDockQ...")
    tqdm.pandas(desc="Basic Metrics")
    pairs["ranking_score"] = pairs.progress_apply(lambda x: read_score(x["GPCR"], x["Gprotein"], args.results_dir), axis=1)
    pairs['iptm'] = pairs.progress_apply(lambda x: read_iptm(x["GPCR"], x["Gprotein"], args.results_dir), axis=1)

    # Calculate interaction and last 7 contact probabilities
    prob_results = pairs.progress_apply(lambda x: get_interaction_and_contact_probs(x["GPCR"], x["Gprotein"], args.results_dir), axis=1, result_type='expand')
    prob_cols = ["interaction_prob", "contact_prob26", "contact_prob25", "contact_prob24", "contact_prob23", "contact_prob22", "contact_prob21", "contact_prob20"]
    pairs[prob_cols] = prob_results

    # Calculate interface energy
    energy_results = pairs.progress_apply(lambda x: pd.Series(find_interface_energy(x["GPCR"], x["Gprotein"], args.results_dir)), axis=1)
    pairs[['dG', 'dG/dSASA', 'dSASA']] = energy_results

    # Calculate pDockQ
    pairs['pDockQ'] = pairs.progress_apply(lambda x: find_pdockq(x["GPCR"], x["Gprotein"], args.results_dir), axis=1)

    # Normalize Basic Metrics
    print("Normalizing basic metrics...")
    cols_to_normalize = pairs.select_dtypes(include=np.number).columns.tolist()

    for feature in tqdm(cols_to_normalize, desc="Normalizing"):
        pairs = normalize(pairs, feature)

    # Prepare for Position-Specific Metrics
    # Filter contacts involving GPCR conserved positions (BW)
    contacts_filtered = contacts[contacts['BW'].isin(valid_gpcr_positions_bw)]
    GPCR_cont_bw = list(contacts_filtered['BW'])
    Gprot_cont_cgn = list(contacts_filtered['Gprot_pos'])

    # Parse BLAST output to create BW to sequence position mapping for the specific GPCR
    print(f"Parsing BLAST output for {gpcr}...")
    mapping_dict, uniac = parse_blastp_output(args.blastp_out)

    # Create BW2pos map using the alignment and BW_data
    BW2pos = {}
    # Find the maximum possible query position from the alignment dict keys
    max_query_pos = max(mapping_dict.keys())
    print(f"Building BW2pos map using alignment (up to query position {max_query_pos})...")
    for i in tqdm(range(1, max_query_pos), desc="BW Mapping"):
        BWposition = BW_mapping(mapping_dict, uniac, i, BW_data)
        if BWposition != "-":
            BW2pos[BWposition] = i

    # Calculate Specific Contact Probabilities
    print("Calculating specific contact probabilities...")
    num_contacts = len(GPCR_cont_bw)
    contact_prob_cols = [f'contact_prob_{GPCR_cont_bw[i]}_{Gprot_cont_cgn[i]}' for i in range(num_contacts)]

    # Apply function row-wise
    contact_results = []
    for index, row in tqdm(pairs.iterrows(), total=len(pairs), desc="Specific Contacts"):
        result = get_interaction_prob_specific(
            row['GPCR'], row['Gprotein'], row['Gprotein_name'],
            GPCR_cont_bw, Gprot_cont_cgn, # BW list, CGN list
            BW2pos, G2pos, # Maps
            args.results_dir,
            mode='contact'
        )
        contact_results.append(result)

    # Assign results
    if contact_results:
        contact_df = pd.DataFrame(contact_results, index=pairs.index, columns=contact_prob_cols)
        pairs[contact_prob_cols] = contact_df

    # Normalize specific contact probabilities
    print("Normalizing specific contact probabilities...")
    for col in tqdm(contact_prob_cols, desc="Normalizing Contacts"):
        pairs = normalize(pairs, col)

    # Calculate pLDDT for Consensus Positions
    print("Calculating pLDDT for consensus positions...")

    plddt_gpcr_cols = [f'plddt_{bw}' for bw in valid_gpcr_positions_bw]
    plddt_gprot_cols = [f'plddt_{cgn}' for cgn in consensus_gprotein_pos_cgn]
    plddt_cols = plddt_gpcr_cols + plddt_gprot_cols

    # Apply function row-wise
    plddt_results_list = []
    for index, row in tqdm(pairs.iterrows(), total=len(pairs), desc="Consensus pLDDT"):
        result = get_interaction_prob_specific(
            row['GPCR'], row['Gprotein'], row['Gprotein_name'],
            valid_gpcr_positions_bw, consensus_gprotein_pos_cgn, # BW list, CGN list
            BW2pos, G2pos, # Maps
            args.results_dir,
            plddts_data=plddts, # Pass preloaded plddt data
            mode='plddt'
        )
        plddt_results_list.append(result)

    # Assign results
    if plddt_results_list:
        plddt_df = pd.DataFrame(plddt_results_list, index=pairs.index, columns=plddt_cols)
        pairs = pd.concat([pairs, plddt_df], axis=1)

    # Normalize pLDDT values
    print("Normalizing pLDDT values...")
    for col in tqdm(plddt_cols, desc="Normalizing pLDDT"):
        pairs = normalize(pairs, col)

    # Save Final Results
    output_filename = args.output_file if args.output_file else f"{gpcr}_pairs.tsv"
    print(f"Saving results to {output_filename}...")
    pairs.to_csv(output_filename, sep="\t", index=False, na_rep='NaN')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze GPCR-GProtein interaction data from modeling results.")

    # --- Input Files ---
    parser.add_argument("--gpcr-id", required=True, help="Identifier for the target GPCR.")
    parser.add_argument("--gproteins-list", default="data/gproteins.txt", help="File listing G-protein Uniprot accessions.")
    parser.add_argument("--plddt-json", required=True, help="Path to the JSON file containing pLDDT scores for all pairs.")
    parser.add_argument("--blastp-out", required=True, help="Path to the BLASTP output file (format 0) for the GPCR sequence against reference.")
    parser.add_argument("--bw-pickle", default="data/GPCRDB_pos_latest.pickle", help="Path to the BW mapping pickle file (GPCRdb positions).")
    parser.add_argument("--contacts-tsv", default="data/contacts.tsv", help="Path to the TSV file defining specific contacts (BW, Gprot_pos).")
    parser.add_argument("--valid-gpcr-pos", default="data/valid_positions.txt", help="File listing valid GPCR BW positions to consider.")
    parser.add_argument("--consensus-gprot-pos", default="data/consensus_Gprotein_positions.txt", help="File listing consensus G-protein positions (CGN format) for pLDDT analysis.")
    parser.add_argument("--gprot-pos-map", default="data/Processed_Gprotein_Position.csv", help="CSV file mapping G-protein common name and CGN to structure position index.")

    # --- Directories ---
    parser.add_argument("--results-dir", default=".", help="Base directory where individual pair results (subdirectories like gpcr_gprot/) are located. Default: current directory.")

    # --- Output File ---
    parser.add_argument("--output-file", help="Optional: Name for the output TSV file. If not provided, defaults to '<gpcr-id>_pairs.tsv'.")

    args = parser.parse_args()
    main(args)