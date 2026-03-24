# Precog3D

Pipeline for **GPCR–G protein coupling prediction** based on AlphaFold3 modeling, structure post-processing, interface metrics extraction, and TabPFN regression described in:

Miglionico, P., Matic, M., Franchini, L., Hiroki, A., Nemati Fard L.A., Arora C., Gherghinescu, M., De Oliveira Rosa, N., Ryoji, K., Gutkind, J. S., Orlandi, C., Inoue, A., Raimondi F. Computed atlas of the human GPCR-G protein signaling complexes. BioRxiv, 2026, https://doi.org/10.64898/2026.03.07.710286 

The repository combines:
- AlphaFold3 input generation from FASTA
- batch structure modeling against a fixed G-protein panel
- structure filtering and interface scoring (pDockQ + Rosetta interface energy)
- feature-table generation
- final coupling-value prediction with TabPFN

## Requirements

## Python dependencies

Install Python packages:

```bash
pip install -r requirements.txt
```

## External tools / runtime environment

The pipeline assumes an HPC-like environment with:
- PBS scheduler (`qsub`, `-Wblock=true`)
- Singularity/Apptainer with an AlphaFold3 container (`alphafold3.sif`)
- AlphaFold3 databases directory
- BLAST+ (`blastp`) with a local GPCR database
- Rosetta binaries (`relax.default.linuxgccrelease`, `InterfaceAnalyzer.linuxgccrelease`)

## Input

- A GPCR sequence in FASTA format (single target)
- AlphaFold3 home/database/output folders
- A folder containing alignment-based JSON resources for the G-protein panel

## Quick start

From repository root:

```bash
bash pipeline.sh \
  -h /path/to/alphafold3_home \
  -d /path/to/alphafold3_databases \
  -i /path/to/TARGET.fasta \
  -o /path/to/output_dir \
  -a /path/to/alignment_dir
```

This launches the full workflow:
1. Convert FASTA to AF3 JSON (`src/fasta_to_json_AF3.py`)
2. Build alignments with AF3 (`src/build_alignment.sh`)
3. Create pairwise JSON jobs (`src/generate_json.py`)
4. Run AF3 structure inference (`src/build_structure.sh`)
5. Filter residues and compute pDockQ (`src/filter.py`, `src/pdockq.py`)
6. Relax structures and compute interface energy (`src/relax.sh`, `src/interface_energy.sh`)
7. Run BLASTP and extract pLDDT profiles (`src/extract_plddt_json.py`)
8. Build feature table (`src/create_table.py`)
9. Final model scoring (`src/model.py`)

## Main outputs

For each GPCR–G-protein pair (`<gpcr>_<gprot>/`):
- `<gpcr>_<gprot>_model.cif`
- `<gpcr>_<gprot>_filtered.pdb`
- `<gpcr>_<gprot>_pdockq.txt`
- `<gpcr>_<gprot>_interface.tsv`

Global outputs:
- `<GPCR>.out` (BLASTP output)
- `<gpcr>_plddt.json` (per-pair pLDDT traces)
- `<GPCR>_table.tsv` (feature table)
- final coupling prediction table produced by `src/model.py`

## Notes

- This project is HPC-oriented, expects scheduler/container infrastructure, shell scripts need to be adjusted to your specific system