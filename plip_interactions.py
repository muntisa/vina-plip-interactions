#!/usr/bin/env python
"""
plip_interactions.py - Automated PLIP analysis pipeline for protein-ligand docking results.
Extracts pose 1 from Vina output, converts to PDB, runs PLIP via Python API, and exports CSV.
"""

import argparse
import sys
import os
import csv
import re
import glob

if getattr(sys, 'frozen', False):
    _meipass = sys._MEIPASS
    os.environ['BABEL_LIBDIR'] = os.path.join(_meipass, 'openbabel', 'bin')
    os.environ['BABEL_DATADIR'] = os.path.join(_meipass, 'openbabel', 'share', 'openbabel', '3.2.0')
else:
    _meipass = None

from openbabel import openbabel as ob
from plip.modules.preparation import PDBComplex
from plip.modules.report import StructureReport

if getattr(sys, 'frozen', False):
    SCRIPT_DIR = os.getcwd()
else:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MOLECULES_DIR = os.path.join(SCRIPT_DIR, 'molecules')
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract the top pose from Vina output, run PLIP, and export interactions as CSV.')
    parser.add_argument('--receptor', '-r',
                        default=os.path.join(MOLECULES_DIR, 'receptor.pdbqt'),
                        help='Receptor PDBQT file (default: molecules/receptor.pdbqt)')
    parser.add_argument('--vina-output', '-v',
                        default=os.path.join(MOLECULES_DIR, 'vina_output.pdbqt'),
                        help='Multi-model Vina docking output PDBQT (default: molecules/vina_output.pdbqt)')
    parser.add_argument('--output-dir', '-o',
                        default=RESULTS_DIR,
                        help='PLIP output directory (default: results)')
    parser.add_argument('--csv',
                        default=os.path.join(RESULTS_DIR, 'plip_interactions.csv'),
                        help='CSV output path (default: results/plip_interactions.csv)')
    return parser.parse_args()


def build_file_paths():
    return {
        'pose_pdbqt': os.path.join(MOLECULES_DIR, 'pose_1.pdbqt'),
        'pose_pdb': os.path.join(MOLECULES_DIR, 'pose_1.pdb'),
        'receptor_pdb': os.path.join(MOLECULES_DIR, 'receptor.pdb'),
        'complex_pdb': os.path.join(MOLECULES_DIR, 'complex.pdb'),
    }


def extract_pose_1(vina_path, pose_path):
    """Extract MODEL 1 from multi-model vina output PDBQT."""
    with open(vina_path, 'r') as f:
        content = f.read()
    parts = re.split(r'(?=^MODEL\s+\d+)', content, flags=re.MULTILINE)
    for part in parts:
        if part.strip().startswith('MODEL 1'):
            with open(pose_path, 'w') as f:
                f.write(part.strip() + '\n')
            print(f"  Extracted to {pose_path}")
            return
    print("ERROR: MODEL 1 not found in vina output")
    sys.exit(1)


def convert_pdbqt_to_pdb(src, dst):
    """Convert PDBQT to PDB using OpenBabel Python API."""
    conv = ob.OBConversion()
    conv.SetInFormat("pdbqt")
    conv.SetOutFormat("pdb")
    mol = ob.OBMol()
    if not conv.ReadFile(mol, src):
        print(f"ERROR: Failed to read {src}")
        sys.exit(1)
    conv.WriteFile(mol, dst)
    print(f"  Converted {os.path.basename(src)} -> {os.path.basename(dst)}")


def combine_complex(receptor_pdb, pose_pdb, complex_pdb):
    """Combine receptor and ligand PDB files into one complex."""
    with open(receptor_pdb) as f:
        rec_lines = [l for l in f if not l.startswith(('COMPND', 'AUTHOR', 'END'))]
    with open(pose_pdb) as f:
        lig_lines = [l for l in f if l.startswith(('ATOM', 'HETATM'))]
    combined = rec_lines + lig_lines + ['END\n']
    with open(complex_pdb, 'w') as f:
        f.writelines(combined)
    print(f"  Complex written: {len(rec_lines)} receptor + {len(lig_lines)} ligand atoms")


def find_fixed_pdb(output_dir):
    """Find the PLIP-fixed PDB file in the output directory."""
    pattern = os.path.join(output_dir, 'plipfixed.complex_*.pdb')
    files = glob.glob(pattern)
    if not files:
        print(f"ERROR: PLIP-fixed PDB not found in {output_dir}")
        sys.exit(1)
    return files[0]


def build_atom_lookup(fixed_pdb):
    """Build atom lookups from the PLIP-fixed PDB file."""
    lig_atoms = {}
    prot_atoms = {}
    lig_count = 0

    with open(fixed_pdb) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                serial = int(line[6:11].strip())
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resnum = int(line[22:26].strip())
                element = line[76:78].strip() or atomname[0]
            except (ValueError, IndexError):
                continue

            if resname == 'UNL':
                if element.upper() != 'H':
                    lig_count += 1
                    lig_atoms[serial] = {'element': element, 'idx': lig_count}
            else:
                prot_atoms[serial] = {
                    'resname': resname, 'resnum': resnum,
                    'chain': chain, 'atomname': atomname,
                }

    return lig_atoms, prot_atoms


def run_plip_analysis(complex_pdb, output_dir):
    """Run PLIP analysis using the PLIP Python API."""
    if os.path.exists(output_dir):
        import shutil
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    mol = PDBComplex()
    mol.output_path = output_dir
    mol.load_pdb(complex_pdb)
    mol.analyze()

    if not mol.interaction_sets:
        print("  No interaction sets found by PLIP")
        return mol, None

    key = list(mol.interaction_sets.keys())[0]
    pli = mol.interaction_sets[key]
    print(f"  PLIP analysis complete: {key}")
    return mol, pli


def _get_prot_idx(sb, protispos):
    """Get the first protein PDB serial from a salt bridge."""
    return sb.positive.atoms_orig_idx[0] if protispos else sb.negative.atoms_orig_idx[0]


def _get_lig_idx(sb, protispos):
    """Get the first ligand PDB serial from a salt bridge."""
    return sb.negative.atoms_orig_idx[0] if protispos else sb.positive.atoms_orig_idx[0]


def extract_interactions(pli, lig_atoms, prot_atoms):
    """Extract all interaction types from PLIP results."""
    interactions = []

    for hyd in pli.hydrophobic_contacts:
        prot = prot_atoms.get(hyd.bsatom_orig_idx, {})
        lig = lig_atoms.get(hyd.ligatom_orig_idx, {})
        interactions.append({
            'type': 'hydrophobic',
            'chain': prot.get('chain', hyd.reschain),
            'aa': f"{hyd.restype}:{hyd.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{hyd.bsatom_orig_idx}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{hyd.distance:.2f}",
        })

    for hb in pli.hbonds_pdon:
        prot = prot_atoms.get(hb.d_orig_idx, {})
        lig = lig_atoms.get(hb.a_orig_idx, {})
        interactions.append({
            'type': 'hydrogen_bond',
            'chain': prot.get('chain', hb.reschain),
            'aa': f"{hb.restype}:{hb.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{hb.d_orig_idx}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{hb.distance_ah:.2f}",
        })

    for sb in pli.saltbridge_lneg + pli.saltbridge_pneg:
        prot_serial = _get_prot_idx(sb, sb.protispos)
        lig_serial = _get_lig_idx(sb, sb.protispos)
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'salt_bridge',
            'chain': prot.get('chain', sb.reschain),
            'aa': f"{sb.restype}:{sb.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{sb.distance:.2f}",
        })

    for ps in pli.pistacking:
        prot_serial = ps.proteinring.atoms_orig_idx[0]
        lig_serial = ps.ligandring.atoms_orig_idx[0]
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'pi_stack',
            'chain': prot.get('chain', ps.reschain),
            'aa': f"{ps.restype}:{ps.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{ps.distance:.2f}",
        })

    for pc in pli.pication_laro + pli.pication_paro:
        prot_serial = pc.ring.atoms_orig_idx[0]
        lig_serial = pc.charge.atoms_orig_idx[0]
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'pi_cation',
            'chain': prot.get('chain', pc.reschain),
            'aa': f"{pc.restype}:{pc.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{pc.distance:.2f}",
        })

    for wb in pli.water_bridges:
        prot_serial = wb.d_orig_idx if wb.protisdon else wb.a_orig_idx
        lig_serial = wb.a_orig_idx if wb.protisdon else wb.d_orig_idx
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'water_bridge',
            'chain': prot.get('chain', wb.reschain),
            'aa': f"{wb.restype}:{wb.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{wb.distance_aw:.2f}",
        })

    for hg in pli.halogen_bonds:
        prot_serial = hg.acc_orig_idx
        lig_serial = hg.don_orig_idx
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'halogen_bond',
            'chain': prot.get('chain', hg.reschain),
            'aa': f"{hg.restype}:{hg.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{hg.distance:.2f}",
        })

    for mc in pli.metal_complexes:
        prot_serial = mc.metal_orig_idx
        lig_serial = mc.target_orig_idx
        prot = prot_atoms.get(prot_serial, {})
        lig = lig_atoms.get(lig_serial, {})
        interactions.append({
            'type': 'metal_complex',
            'chain': prot.get('chain', mc.reschain),
            'aa': f"{mc.restype}:{mc.resnr}",
            'prot_atom': f"{prot.get('atomname', '?')}:{prot_serial}",
            'lig_name': 'UNL',
            'lig_atom': f"{lig.get('element', '?')}{lig.get('idx', '?')}",
            'distance': f"{mc.distance:.2f}",
        })

    return interactions


def write_csv(interactions, csv_path):
    """Write interactions to CSV file."""
    fieldnames = [
        'interaction_type', 'protein_chain', 'amino_acid',
        'protein_atom', 'ligand_name', 'ligand_atom', 'distance_angstrom'
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for item in interactions:
            writer.writerow({
                'interaction_type': item['type'],
                'protein_chain': item['chain'],
                'amino_acid': item['aa'],
                'protein_atom': item['prot_atom'],
                'ligand_name': item['lig_name'],
                'ligand_atom': item['lig_atom'],
                'distance_angstrom': item['distance'],
            })
    print(f"  CSV written: {len(interactions)} interactions")


def main():
    args = parse_args()
    paths = build_file_paths()
    receptor_pdbqt = args.receptor
    vina_output = args.vina_output
    plip_output = args.output_dir
    csv_path = args.csv

    if not os.path.isfile(receptor_pdbqt):
        print(f"\nERROR: Receptor file not found: {receptor_pdbqt}")
        print("Use --receptor (-r) to specify the correct path.")
        print("Run with --help for usage information.")
        sys.exit(1)

    if not os.path.isfile(vina_output):
        print(f"\nERROR: Vina output file not found: {vina_output}")
        print("Use --vina-output (-v) to specify the correct path.")
        print("Run with --help for usage information.")
        sys.exit(1)

    os.makedirs(MOLECULES_DIR, exist_ok=True)
    os.makedirs(plip_output, exist_ok=True)

    step = 0
    total = 5

    step += 1
    print(f"\n[{step}/{total}] Extracting pose 1 from {os.path.basename(vina_output)}...")
    extract_pose_1(vina_output, paths['pose_pdbqt'])

    step += 1
    print(f"\n[{step}/{total}] Converting PDBQT to PDB...")
    convert_pdbqt_to_pdb(paths['pose_pdbqt'], paths['pose_pdb'])
    convert_pdbqt_to_pdb(receptor_pdbqt, paths['receptor_pdb'])

    step += 1
    print(f"\n[{step}/{total}] Combining receptor + ligand into complex.pdb...")
    combine_complex(paths['receptor_pdb'], paths['pose_pdb'], paths['complex_pdb'])

    step += 1
    print(f"\n[{step}/{total}] Running PLIP analysis...")
    mol, pli = run_plip_analysis(paths['complex_pdb'], plip_output)

    if pli is None or pli.no_interactions:
        print("\nNo interactions found. Writing empty CSV.")
        write_csv([], csv_path)
        return

    step += 1
    print(f"\n[{step}/{total}] Parsing PLIP results and writing reports...")
    fixed_pdb = find_fixed_pdb(plip_output)
    lig_atoms, prot_atoms = build_atom_lookup(fixed_pdb)
    interactions = extract_interactions(pli, lig_atoms, prot_atoms)
    write_csv(interactions, csv_path)

    csv_dir = os.path.dirname(os.path.abspath(csv_path))
    report = StructureReport(mol)
    report.outpath = csv_dir
    report.write_xml()
    report.write_txt()
    print(f"  Reports written: report.xml, report.txt")

    print(f"\nAll done! {len(interactions)} interactions found.")
    for i in interactions:
        d = i['distance'] if i['distance'] else ''
        print(f"  {i['type']:20s} {i['aa']:8s}({i['chain']}) "
              f"{i['prot_atom']:12s} <-> {i['lig_atom']:5s}  {d}A")


if __name__ == '__main__':
    main()
