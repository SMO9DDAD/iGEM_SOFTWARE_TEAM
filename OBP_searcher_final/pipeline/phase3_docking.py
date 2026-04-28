"""
phase3_docking.py — Funcions auxiliars compartides de docking
Usades per main2.py (Fase C) i phase4_selectivity.py (Fase D).
"""

import os
from pathlib import Path


def prepare_ligand(smiles: str, out_pdbqt: str):
    """
    Converteix un SMILES a format PDBQT per a AutoDock Vina.
    Passos: SMILES 2D -> geometria 3D (ETKDGv3) -> optimització (MMFF94) -> PDBQT (Meeko)
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)                          # afegeix hidrògens explícits
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())  # genera geometria 3D
    AllChem.MMFFOptimizeMolecule(mol)              # optimitza amb MMFF94
    prep = MoleculePreparation()
    prep.prepare(mol)
    prep.write_pdbqt_file(out_pdbqt)


def prepare_receptor(pdb_path: str, pdbqt_path: str):
    """
    Converteix un fitxer PDB a format PDBQT per a AutoDock Vina.
    Usa Open Babel amb flag -xr (receptor rígid): elimina aigua,
    afegeix hidrògens polars i assigna tipus d'àtom AutoDock4.
    """
    os.system(f"obabel {pdb_path} -O {pdbqt_path} -xr 2>nul")


def get_center(pdb_path: str):
    """
    Calcula el centre geomètric de tots els àtoms pesants de la proteïna.
    Retorna (cx, cy, cz) en Angstroms per centrar la caixa de docking.
    """
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except:
                    pass
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)