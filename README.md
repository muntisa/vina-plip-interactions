# Vina-PLIP Interactions
**A standalone Windows executable for automated AutoDock Vina pose extraction and PLIP interaction profiling**

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Automated pipeline that extracts the top-ranked docking pose from AutoDock Vina output, converts it to PDB, builds the protein–ligand complex, runs [PLIP](https://github.com/pharmbio/plip) (Protein–Ligand Interaction Profiler), and exports all detected non-covalent interactions as a clean CSV table along with detailed report files (text and XML).

## Features

- Extracts the best pose (`MODEL 1`) from a multi-model Vina PDBQT file
- Converts PDBQT to standard PDB via OpenBabel Python API
- Combines receptor and ligand into a single complex PDB
- Runs PLIP via Python API (no subprocess calls) to detect hydrophobic contacts, hydrogen bonds, salt bridges, π-stacking, π-cation, water bridges, halogen bonds, and metal complexes
- Outputs a human-readable CSV (protein-as-donor only for hydrogen bonds) plus detailed PLIP report files (`report.txt` and `report.xml`)
- Available as a single-file Windows executable (no Python installation required)

## Pre-built Executable (Windows)

A single-file executable is available as `plip_interactions.exe`. No Python installation or virtual environment is needed — just run it from the command line:

```bash
# Run with default paths (molecules/ and results/ in current directory)
plip_interactions.exe

# Specify custom paths
plip_interactions.exe -r protein.pdbqt -v docking.pdbqt -o results --csv results/interactions.csv
```

All arguments are identical to the Python script (see [Arguments](#arguments) below).

### Building the EXE yourself

If you need to rebuild the executable (e.g., after modifying the script):

```bash
python -m venv .venv
.venv\Scripts\activate
pip install openbabel==3.2.0 plip==1.4.2 pyinstaller

# Build
pyinstaller --onefile --console --name plip_interactions \
    --add-data ".venv\Lib\site-packages\openbabel\bin\*.obf;openbabel\bin" \
    --add-data ".venv\Lib\site-packages\openbabel\share;openbabel\share" \
    plip_interactions.py
```

The resulting `dist/plip_interactions.exe` (~24 MB) contains all dependencies.

## Installation (Python)

```bash
# Clone the repository
git clone https://github.com/muntisa/vina-plip-interactions.git
cd vina-plip-interactions

# Create and activate a virtual environment
python -m venv .venv
# Windows:
.venv\Scripts\activate
# Linux/macOS:
# source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

Place your input files in the `molecules/` folder:

| File | Description |
|------|-------------|
| `molecules/receptor.pdbqt` | Protein receptor prepared for AutoDock/Vina (standard amino acids only) |
| `molecules/vina_output.pdbqt` | Multi-model Vina docking output with `MODEL 1` as the top pose |

Run with defaults:

```bash
python plip_interactions.py
```

Or specify custom paths:

```bash
python plip_interactions.py \
    --receptor protein.pdbqt \
    --vina-output docking_results.pdbqt \
    --output-dir results \
    --csv results/interactions.csv
```

### Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--receptor` | `-r` | `molecules/receptor.pdbqt` | Receptor PDBQT file (standard amino acids only) |
| `--vina-output` | `-v` | `molecules/vina_output.pdbqt` | Multi-model Vina output PDBQT |
| `--output-dir` | `-o` | `results` | Directory for PLIP output files |
| `--csv` | | `results/plip_interactions.csv` | Path for the CSV output |

## Pipeline

1. **Extract** — `MODEL 1` is split from the Vina output PDBQT into `molecules/pose_1.pdbqt`
2. **Convert** — Both ligand and receptor are converted from PDBQT → PDB via OpenBabel Python API
3. **Combine** — Receptor PDB + ligand ATOM records → `molecules/complex.pdb`
4. **Analyze** — PLIP runs via Python API on `molecules/complex.pdb`, outputs to `results/`
5. **Parse** — Interaction data is extracted directly from PLIP's in-memory results
6. **Export** — All interactions are written to `results/plip_interactions.csv`, `results/report.txt`, and `results/report.xml`

## Output

| File | Description |
|------|-------------|
| `results/report.txt` | PLIP full text report |
| `results/report.xml` | PLIP XML report with per-interaction geometry |
| `results/plipfixed.complex_*.pdb` | PLIP-fixed PDB file used for atom-level resolution |
| `results/plip_interactions.csv` | Tabular interaction summary |

### CSV columns

| Column | Example | Description |
|--------|---------|-------------|
| `interaction_type` | `hydrophobic` | Interaction category |
| `protein_chain` | `A` | Protein chain identifier |
| `amino_acid` | `GLU:40` | Residue name and sequence number |
| `protein_atom` | `OE2:297` | Atom name and PDB serial in the protein |
| `ligand_name` | `UNL` | Ligand residue name |
| `ligand_atom` | `O20` | Element symbol and heavy-atom index in the ligand |
| `distance_angstrom` | `3.15` | Distance in angstroms |

### Parsed interaction types

| Type | PLIP API attribute |
|------|---------------------|
| Hydrophobic | `hydrophobic_contacts` |
| Hydrogen bond | `hbonds_pdon` (protein as donor) |
| Salt bridge | `saltbridge_lneg` + `saltbridge_pneg` |
| π-Stacking | `pistacking` |
| π-Cation | `pication_laro` + `pication_paro` |
| Water bridge | `water_bridges` |
| Halogen bond | `halogen_bonds` |
| Metal complex | `metal_complexes` |

## Notes

- **Receptor PDBQT** must contain only standard amino acids. Non-standard residues (LIG, UNK) will be treated as ligands by PLIP.
- **Hydrogen bonds** are reported once with the protein as the donor. The reciprocal ligand-as-donor entry is omitted.
- **Intermediate files** (`molecules/pose_1.pdbqt`, `molecules/pose_1.pdb`, `molecules/receptor.pdb`, `molecules/complex.pdb`) are overwritten on each run.
- The script was developed and tested on **Windows** with Python 3.10. Paths may need adjustment for Linux/macOS.
- **OpenBabel plugin warning** — The "Unable to find OpenBabel plugins" message from the EXE is cosmetic and does not affect functionality.

## License

MIT
