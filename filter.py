from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import argparse


def filter_residues(input_cif, output_pdb, b_factor_threshold=70.0):
    """Filters residues based on B-factor, preserving last 7 residues of chain B."""

    # Structure initialization
    parser = MMCIFParser()
    structure = parser.get_structure("protein", input_cif)
    cif_dict = MMCIF2Dict(input_cif)

    # Identify chain B residues
    chain_b_residues = []
    for residue in structure[0]['B']:
        chain_b_residues.append(residue.id)
    last_seven = chain_b_residues[-7:] if len(chain_b_residues) >=7 else chain_b_residues

    # Filtering criteria
    keep_residues = set()
    for i, (atom, b_val, chain, seq_id, ins_code) in enumerate(zip(
        cif_dict["_atom_site.label_atom_id"],
        cif_dict["_atom_site.B_iso_or_equiv"],
        cif_dict["_atom_site.auth_asym_id"],
        cif_dict["_atom_site.label_seq_id"],
        cif_dict["_atom_site.pdbx_PDB_ins_code"]
    )):
        if atom != "CA":
            continue
        
        residue_id = (' ', int(seq_id), " " if ins_code == "?" else ins_code)
        b_val = float(b_val)

        if chain == 'B' and residue_id in last_seven:
            keep_residues.add((chain, residue_id))
        elif b_val >= b_factor_threshold:
            keep_residues.add((chain, residue_id))

    # Create filtered structure
    filtered = structure.copy()
    for chain in list(filtered[0].get_chains()):
        for res in list(chain.get_residues()):
            if (chain.id, res.id) not in keep_residues:
                chain.detach_child(res.id)

    # Save as PDB
    pdb_io = PDBIO()
    pdb_io.set_structure(filtered)
    pdb_io.save(output_pdb)
    print(f"Filtered structure saved to {output_pdb}")

def main():
    parser = argparse.ArgumentParser(description="Filter GPCR-Gprotein complex structures")
    parser.add_argument("--input_cif", help="Input CIF file")
    parser.add_argument("--output_pdb", help="Output PDB file")
    parser.add_argument("--b_factor_threshold", type=float, default=70.0, help="B-factor threshold for filtering")
    args = parser.parse_args()
    
    filter_residues(args.input_cif, args.output_pdb, args.b_factor_threshold)

if __name__ == "__main__":
    main()