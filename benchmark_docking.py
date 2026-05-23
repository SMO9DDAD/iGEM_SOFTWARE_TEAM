"""
benchmark_docking.py
====================
Fa docking de 8 VOCs seleccionats (els amb més cobertura experimental a iobpdb)
contra les seves TOP 5 OBPs i compara el score de Vina amb el Kd experimental.

Els VOCs i els seus SMILES estan hardcodejats per evitar dependre de PubChem.

Ús:
    python benchmark_docking.py

Fitxers necessaris:
    - data/Compound_OBP_binding.csv
    - vina/vina_1.2.7_win.exe
"""

import os
import re
import math
import time
import subprocess
import requests
import pandas as pd
from pathlib import Path

# ─── CONFIGURACIÓ ─────────────────────────────────────────────────────────────
IOBPDB_PATH    = "data/Compound_OBP_binding.csv"
VINA_EXE       = str(Path("vina/vina_1.2.7_win.exe"))
TOP_N          = 5
RESULTS_DIR    = "results"
DOCKING_DIR    = "docking/benchmark"
IOBPDB_PDB_URL = "https://raw.githubusercontent.com/sshuklz/iobpdb_app/master/alpha_pdbs/{}.pdb"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(DOCKING_DIR, exist_ok=True)

# VOCs SELECCIONATS 


VOC_DATA = [
    {
        "name":       "beta-ionone",
        "smiles":     "CC1=CCCC(=C1/C=C/C(=O)C)C",   # (E)-beta-ionone
        "iobpdb_key": "ionone (beta)",
    },
    {
        "name":       "cis-3-hexen-1-ol",
        "smiles":     "OCC/C=C\\CC",
        "iobpdb_key": "cis-3-Hexen-1-ol",
    },
    {
        "name":       "benzaldehyde",
        "smiles":     "O=Cc1ccccc1",
        "iobpdb_key": "benzaldehyde",
    },
    {
        "name":       "beta-myrcene",
        "smiles":     "C=CC(=C)CCC=C(C)C",
        "iobpdb_key": "beta-myrcene",
    },
    {
        "name":       "nonanal",
        "smiles":     "CCCCCCCCC=O",
        "iobpdb_key": "nonanal",
    },
    {
        "name":       "linalool",
        "smiles":     "CC(C)=CCCC(C)(O)C=C",
        "iobpdb_key": "linalool",
    },
    {
        "name":       "methyl-salicylate",
        "smiles":     "COC(=O)c1ccccc1O",
        "iobpdb_key": "methyl salicylate",
    },
    {
        "name":       "hexanal",
        "smiles":     "CCCCCC=O",
        "iobpdb_key": "hexanal",
    },
]


#  FUNCIONS REUTILITZADES

def is_mutant(obp_name: str) -> bool:
    return bool(re.search(r'\s+[A-Z]\d+[A-Z]', obp_name))


def get_top_obps(iobpdb_key: str, df: pd.DataFrame, top_n: int) -> list:
    """Retorna les TOP N OBPs naturals amb Kd mesurat per a un VOC."""
    mask = df["Compound name"].str.contains(iobpdb_key, case=False, na=False)
    row  = df[mask]
    if row.empty:
        return []

    results = []
    for obp_name in df.columns[2:]:
        if is_mutant(obp_name):
            continue
        val = row.iloc[0][obp_name]
        if pd.isna(val):
            continue
        val_str = str(val).strip()
        if val_str.startswith(">"):
            continue
        try:
            kd_nm = int(round(float(val_str) * 1000))
            results.append({"obp_name": obp_name, "kd_nm_experimental": kd_nm})
        except:
            continue

    return sorted(results, key=lambda x: x["kd_nm_experimental"])[:top_n]


def download_iobpdb_structure(obp_name: str, out_pdb: str) -> bool:
    try:
        r = requests.get(IOBPDB_PDB_URL.format(obp_name), timeout=30)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            return True
    except:
        pass
    return False


def get_uniprot_id(obp_name: str):
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?query={obp_name}&format=json&size=1",
            timeout=10)
        results = r.json().get("results", [])
        if results:
            return results[0]["primaryAccession"]
    except:
        pass
    return None


def get_sequence(uniprot_id: str):
    try:
        r = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta", timeout=15)
        if r.status_code == 200:
            lines = r.text.strip().split("\n")
            return "".join(lines[1:])
    except:
        pass
    return None


def predict_structure(seq: str, out_pdb: str, max_len: int = 400) -> bool:
    if len(seq) > max_len:
        seq = seq[20:]
    for intent in range(3):
        r = requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/",
                          data=seq, timeout=180)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            return True
        print(f"      ESMFold intent {intent+1} error {r.status_code}, esperant...")
        time.sleep(10)
    return False


def prepare_ligand(smiles: str, out_pdbqt: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    prep = MoleculePreparation()
    prep.prepare(mol)
    prep.write_pdbqt_file(out_pdbqt)


def prepare_receptor(pdb_path: str, pdbqt_path: str):
    os.system(f"obabel {pdb_path} -O {pdbqt_path} -xr 2>nul")


def get_center(pdb_path: str):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except:
                    pass
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)


def run_vina(receptor: str, ligand: str, out: str,
             cx: float, cy: float, cz: float) -> float | None:
    cmd = [
        VINA_EXE,
        "--receptor", receptor, "--ligand", ligand, "--out", out,
        "--center_x", str(round(cx, 2)),
        "--center_y", str(round(cy, 2)),
        "--center_z", str(round(cz, 2)),
        "--size_x", "25", "--size_y", "25", "--size_z", "25",
        "--exhaustiveness", "8", "--num_modes", "5",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in result.stdout.split("\n"):
        if line.strip().startswith("1"):
            try:
                score = float(line.split()[1])
                return score if score < 0 else None
            except:
                pass
    return None


def score_to_kd_nm(score_kcal: float) -> float:
    """ΔG = RT·ln(Kd)  →  Kd = exp(ΔG/RT) · 1e9 nM"""
    RT = 0.5924  # kcal/mol a 25°C
    return math.exp(score_kcal / RT) * 1e9


def get_structure(obp_name: str, pdb_path: str) -> bool:
    """Obté estructura: primer iobpdb, després ESMFold."""
    if os.path.exists(pdb_path):
        return True
    if download_iobpdb_structure(obp_name, pdb_path):
        print(f"      ✓ Descarregat d'iobpdb")
        return True
    print(f"      Intentant ESMFold...")
    uid = get_uniprot_id(obp_name)
    if not uid:
        return False
    seq = get_sequence(uid)
    if not seq:
        return False
    return predict_structure(seq, pdb_path)


#   PRINCIPAL 

def run_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    all_results = []

    for voc in VOC_DATA:
        voc_name    = voc["name"]
        smiles      = voc["smiles"]
        iobpdb_key  = voc["iobpdb_key"]

        print(f"\n{'='*55}")
        print(f"  VOC: {voc_name.upper()}")
        print(f"  SMILES: {smiles}")
        print(f"{'='*55}")

        # TOP N OBPs experimentals
        top_obps = get_top_obps(iobpdb_key, df, TOP_N)
        if not top_obps:
            print(f"  ⚠ Sense dades experimentals, saltant.")
            continue
        print(f"  TOP {len(top_obps)} OBPs:")
        for o in top_obps:
            print(f"    {o['obp_name']:<22} Kd = {o['kd_nm_experimental']} nM")

        # Prepara lligand (una sola vegada per VOC)
        slug         = voc_name.replace("-", "_")
        ligand_pdbqt = f"{DOCKING_DIR}/ligand_{slug}.pdbqt"
        if not os.path.exists(ligand_pdbqt):
            try:
                prepare_ligand(smiles, ligand_pdbqt)
                print(f"  ✓ Lligand preparat")
            except Exception as e:
                print(f"  ⚠ Error preparant lligand: {e}")
                continue

        # Docking per cada OBP
        for obp in top_obps:
            obp_name   = obp["obp_name"]
            kd_exp     = obp["kd_nm_experimental"]
            pdb_path   = f"{DOCKING_DIR}/{obp_name}.pdb"
            pdbqt_path = f"{DOCKING_DIR}/{obp_name}.pdbqt"
            out_path   = f"{DOCKING_DIR}/{obp_name}_{slug}_out.pdbqt"

            print(f"\n  [{obp_name}]  Kd exp = {kd_exp} nM")

            if not get_structure(obp_name, pdb_path):
                print(f"      ⚠ Sense estructura, saltant.")
                continue

            prepare_receptor(pdb_path, pdbqt_path)
            cx, cy, cz = get_center(pdb_path)
            score = run_vina(pdbqt_path, ligand_pdbqt, out_path, cx, cy, cz)

            if score is None:
                print(f"      ⚠ Docking invàlid, saltant.")
                continue

            kd_pred    = score_to_kd_nm(score)
            log_kd_exp  = math.log10(kd_exp)  if kd_exp  > 0 else None
            log_kd_pred = math.log10(kd_pred) if kd_pred > 0 else None
            if log_kd_exp and log_kd_pred:
                error_pct = (log_kd_exp - log_kd_pred) / log_kd_exp * 100
            else:
                error_pct = None

            print(f"      Score:    {score:.3f} kcal/mol")
            print(f"      Kd pred:  {kd_pred:.0f} nM")
            print(f"      Kd exp:   {kd_exp} nM")
            print(f"      Error:    {error_pct:.1f}%" if error_pct is not None else "      Error:    N/A")

            all_results.append({
                "voc":                voc_name,
                "obp_name":           obp_name,
                "kd_nm_experimental": kd_exp,
                "docking_score_kcal": score,
                "kd_nm_predicted":    round(kd_pred, 1),
                "error_pct":          round(error_pct, 2) if error_pct is not None else None,
            })

    return pd.DataFrame(all_results)


# MAIN 

if __name__ == "__main__":

    print("  BENCHMARK DOCKING")
    print(f"  {len(VOC_DATA)} VOCs  x  TOP {TOP_N} OBPs per VOC")


    df_iobpdb = pd.read_csv(IOBPDB_PATH)

    print(f"\n  VOCs seleccionats:")
    for v in VOC_DATA:
        print(f"    {v['name']:<22} {v['smiles']}")

    # Executa benchmark
    df_results = run_benchmark(df_iobpdb)

    if df_results.empty:
        print("\n⚠ Sense resultats. Revisa Vina i les dependències.")
        exit()


    csv_path = f"{RESULTS_DIR}/benchmark_resultats.csv"
    df_results.to_csv(csv_path, index=False)
    print(f"\n  ✓ Resultats guardats a {csv_path}")


    print(f"\n{'='*55}")
    print("  RESUM GLOBAL")
    print(f"{'='*55}")
    print(f"  Parells VOC-OBP:    {len(df_results)}")
    print(f"  Error % medià:      {df_results['error_pct'].median():.1f}%")
    print(f"  Error % mitjà:      {df_results['error_pct'].mean():.1f}%")
    print(f"\n  Per VOC:")
    for voc, grp in df_results.groupby("voc"):
        print(f"    {voc:<22}  n={len(grp)}  "
              f"error medià={grp['error_pct'].median():.1f}%  "
              f"score mitjà={grp['docking_score_kcal'].mean():.2f} kcal/mol")

    print(f"\n  ✓ BENCHMARK COMPLETAT")    