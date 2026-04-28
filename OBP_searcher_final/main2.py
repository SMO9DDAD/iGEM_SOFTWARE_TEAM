import os
import time
import requests
import subprocess
import pandas as pd
from pathlib import Path
from pipeline.phase0_resolve import resolve_voc
from pipeline.phase4_selectivity import (
    load_interferents, compute_selectivity, ranking_selectivity
)

# ─── CONFIGURACIÓ ────────────────────────────────────────────
VOC_NAME         = "1-octen-3-ol"   # canvia aquí el VOC que vulguis
TOP_N_DOCKING    = 5
VINA_EXE         = str(Path("vina/vina_1.2.7_win.exe"))
IOBPDB_PATH      = "data/Compound_OBP_binding.csv"  # CSV complet iobpdb

os.makedirs("results", exist_ok=True)
os.makedirs("docking",  exist_ok=True)


# FASE 0 — Resolució del VOC

def fase0(voc_name):
    print("\n" + "="*55)
    print("FASE 0 — Resolució del VOC")
    print("="*55)
    try:
        voc = resolve_voc(voc_name)
        if not voc["smiles"]:
            raise ValueError("SMILES buit")
    except Exception as e:
        print(f"  Avís: PubChem no ha respost ({e}).")
        smiles = input("  Introdueix el SMILES manualment: ").strip()
        voc = {"name": voc_name, "smiles": smiles, "inchikey": ""}
    print(f"  SMILES: {voc['smiles']}")
    return voc


# FASE A — Cerca directa al CSV complet d'iobpdb

def fase_a(voc_name):
    print("\n" + "="*55)
    print("FASE A — Cerca a iobpdb")
    print("="*55)

    df = pd.read_csv(IOBPDB_PATH)

    # Busca el VOC per nom (insensible a majúscules)
    mask = df["Compound name"].str.contains(voc_name, case=False, na=False)
    row = df[mask]

    if row.empty:
        # Intent 2: busca per CAS si el nom no coincideix
        print(f"  '{voc_name}' no trobat per nom, intentant per CAS...")
        mask2 = df["CAS-number"].astype(str).str.contains(voc_name, na=False)
        row = df[mask2]

    if row.empty:
        print(f"  ⚠ VOC no trobat a iobpdb. Només es farà docking si hi ha estructures prèvies.")
        return []

    compound_name = row.iloc[0]["Compound name"]
    print(f"  Trobat: '{compound_name}'")

    # Extreu tots els valors de Kd per aquest VOC
    results = []
    for obp_name in df.columns[2:]:
        val = row.iloc[0][obp_name]
        if pd.isna(val):
            continue
        val_str = str(val).strip()
        if val_str.startswith(">"):
            continue
        try:
            kd_nm = int(round(float(val_str) * 1000))
            results.append({
                "obp_name": obp_name,
                "kd_nm":    kd_nm,
                "source":   "iobpdb"
            })
        except:
            continue

    results = sorted(results, key=lambda x: x["kd_nm"])
    print(f"  Trobades {len(results)} OBPs amb Kd mesurat")
    print(f"\n  TOP 10 per Kd:")
    for r in results[:10]:
        print(f"    {r['obp_name']:<20} Kd = {r['kd_nm']} nM")

    pd.DataFrame(results).to_csv("results/fase_a.csv", index=False)
    return results


# FASE C — Docking del VOC diana (TOP N candidats)

def get_sequence(uniprot_id):
    r = requests.get(
        f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta", timeout=15)
    if r.status_code != 200:
        return None
    lines = r.text.strip().split("\n")
    return "".join(lines[1:])

def predict_structure(seq, out_pdb, max_len=400):
    if len(seq) > max_len:
        seq = seq[20:]
    for intent in range(3):
        r = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            data=seq, timeout=180)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            return True
        print(f"    ESMFold intent {intent+1} error {r.status_code}, esperant...")
        time.sleep(10)
    return False

def prepare_receptor(pdb_path, pdbqt_path):
    os.system(f"obabel {pdb_path} -O {pdbqt_path} -xr 2>nul")

def prepare_ligand(smiles, out_pdbqt):
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

def get_center(pdb_path):
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

def run_vina(receptor, ligand, out, cx, cy, cz):
    cmd = [
        VINA_EXE,
        "--receptor", receptor,
        "--ligand",   ligand,
        "--out",      out,
        "--center_x", str(round(cx, 2)),
        "--center_y", str(round(cy, 2)),
        "--center_z", str(round(cz, 2)),
        "--size_x",   "25",
        "--size_y",   "25",
        "--size_z",   "25",
        "--exhaustiveness", "8",
        "--num_modes", "5",
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

def get_uniprot_id(obp_name):
    """Busca l'UniProt ID automàticament a UniProt REST API."""
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?query={obp_name}&format=json&size=1",
            timeout=10)
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                return results[0]["primaryAccession"]
    except:
        pass
    return None

def fase_c(obp_results, voc_smiles):
    print("\n" + "="*55)
    print(f"FASE C — Docking VOC diana (TOP {TOP_N_DOCKING})")
    print("="*55)

    # Prepara el lligand diana una sola vegada
    ligand_pdbqt = "docking/ligand_diana.pdbqt"
    if not os.path.exists(ligand_pdbqt):
        print("  Preparant lligand diana...")
        prepare_ligand(voc_smiles, ligand_pdbqt)
        print("  ✓ Lligand preparat")

    # Selecciona TOP N per Kd experimental
    top_obps = obp_results[:TOP_N_DOCKING]

    docking_results = []
    for obp in top_obps:
        name = obp["obp_name"]
        print(f"\n  [{name}]")

        pdb_path   = f"docking/{name}.pdb"
        pdbqt_path = f"docking/{name}.pdbqt"
        out_path   = f"docking/{name}_diana_out.pdbqt"

        # Obté estructura si no existeix
        if not os.path.exists(pdb_path):
            uid = get_uniprot_id(name)
            if not uid:
                print(f"    ⚠ UniProt ID no trobat, saltant...")
                continue
            print(f"    UniProt: {uid}")
            seq = get_sequence(uid)
            if not seq:
                print(f"    ⚠ Seqüència no trobada, saltant...")
                continue
            print(f"    Predint estructura ({len(seq)} aa)...")
            if not predict_structure(seq, pdb_path):
                print(f"    ⚠ ESMFold ha fallat, saltant...")
                continue
            print(f"    ✓ Estructura predita")

        prepare_receptor(pdb_path, pdbqt_path)
        cx, cy, cz = get_center(pdb_path)
        score = run_vina(pdbqt_path, ligand_pdbqt, out_path, cx, cy, cz)

        if score:
            print(f"    ✓ Score diana: {score:.3f} kcal/mol")
            docking_results.append({
                "obp_name":           name,
                "kd_nm_experimental": obp["kd_nm"],
                "docking_score_kcal": score,
            })

    return docking_results


# RANKING FINAL sense selectivitat 

def ranking_final(docking_results):
    print("\n" + "="*55)
    print("RANKING FINAL (sense selectivitat)")
    print("="*55)
    if not docking_results:
        print("  Cap resultat de docking disponible.")
        return

    df = pd.DataFrame(docking_results)
    df["score_norm"] = (df["docking_score_kcal"] - df["docking_score_kcal"].max()) / \
                       (df["docking_score_kcal"].min() - df["docking_score_kcal"].max())
    df["kd_norm"]    = (df["kd_nm_experimental"].max() - df["kd_nm_experimental"]) / \
                       (df["kd_nm_experimental"].max() - df["kd_nm_experimental"].min() + 1e-9)
    df["score_combinat"] = (0.5 * df["score_norm"] + 0.5 * df["kd_norm"]).round(3)
    df = df.sort_values("score_combinat", ascending=False)

    print(f"\n  {'#':<3} {'OBP':<20} {'Kd exp (nM)':<14} {'Docking':<14} {'Score'}")
    print("  " + "-"*65)
    for i, row in enumerate(df.itertuples(), 1):
        print(f"  {i:<3} {row.obp_name:<20} {row.kd_nm_experimental:<14} "
              f"{row.docking_score_kcal:<14.3f} {row.score_combinat:.3f}")

    df.to_csv("results/ranking_final.csv", index=False)
    print(f"\n  ✓ Ranking guardat a results/ranking_final.csv")
    print(f"\n  → MILLOR CANDIDAT: {df.iloc[0]['obp_name']}")



# MAIN

if __name__ == "__main__":
    print("\n" + "█"*55)
    print("  OBP SEARCHER — Pipeline complet")
    print(f"  VOC: {VOC_NAME}")
    print(f"  TOP {TOP_N_DOCKING} OBPs per docking")
    print("█"*55)

    # Fase 0 — Resolució del VOC
    voc = fase0(VOC_NAME)

    # Fase A — Cerca directa a iobpdb
    obp_results = fase_a(VOC_NAME)

    if not obp_results:
        print("\n⚠ Cap OBP trobat al CSV. Comprova el nom del VOC.")
        exit()

    # Fase C — Docking del VOC diana
    docking_results = fase_c(obp_results, voc["smiles"])

    if not docking_results:
        print("\n⚠ Cap resultat de docking. Revisa les estructures.")
        exit()

    # Fase D — Selectivitat o ranking simple
    interferents = load_interferents("interferents.txt")
    if interferents:
        print("\n" + "="*55)
        print("FASE D — Selectivitat")
        print("="*55)
        results_finals = compute_selectivity(
            docking_results, interferents, voc["smiles"])
        ranking_selectivity(results_finals)
    else:
        ranking_final(docking_results)