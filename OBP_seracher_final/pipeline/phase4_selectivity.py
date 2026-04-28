"""
phase4_selectivity.py
Calcula la selectivitat de cada OBP candidat respecte a una llista
de VOCs interferents. Fa docking real de cada interferent contra
cada OBP i calcula el Selectivity Index (SI).

SI = Kd_interferent_mig / Kd_diana
   (SI > 1 = OBP prefereix el diana, com més gran millor)
"""

import os
import subprocess
import time
import requests
import pandas as pd
from pathlib import Path
from pipeline.phase0_resolve import resolve_voc
from pipeline.phase3_docking import prepare_ligand, prepare_receptor, get_center

VINA_EXE = str(Path("vina/vina_1.2.7_win.exe"))


def load_interferents(path: str = "interferents.txt") -> list:
    """Carrega la llista de VOCs interferents del fitxer de configuració."""
    if not Path(path).exists():
        print(f"  ⚠ No s'ha trobat {path}. Creant fitxer d'exemple...")
        with open(path, "w") as f:
            f.write("# VOCs interferents (un per línia)\n")
            f.write("limonene\nethanol\nhexanal\nlinalool\nalpha-pinene\n")
        return []

    interferents = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                interferents.append(line)

    print(f"  Carregats {len(interferents)} VOCs interferents de {path}")
    return interferents


def run_vina_score(receptor_pdbqt, ligand_pdbqt, out_pdbqt, cx, cy, cz):
    cmd = [
        str(Path("vina/vina_1.2.7_win.exe")),
        "--receptor", receptor_pdbqt,
        "--ligand",   ligand_pdbqt,
        "--out",      out_pdbqt,
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


def score_to_kd_nm(score_kcal: float) -> float:
    """Converteix score docking (kcal/mol) a Kd aproximat (nM)."""
    import math
    RT = 0.5924  # kcal/mol a 25°C
    return math.exp(score_kcal / RT) * 1e9


def compute_selectivity(obp_results, interferents, voc_diana_smiles):
    os.makedirs("docking/interferents", exist_ok=True)

    # Resol SMILES de cada interferent
    print(f"\n  Resolent {len(interferents)} VOCs interferents...")
    interferent_data = []
    for name in interferents:
        try:
            voc = resolve_voc(name)
            if voc["smiles"]:
                interferent_data.append(voc)
                print(f"    ✓ {name}: {voc['smiles']}")
            else:
                print(f"    ⚠ {name}: SMILES no trobat, saltant")
        except Exception as e:
            print(f"    ⚠ {name}: error ({e}), saltant")
        time.sleep(0.2)

    if not interferent_data:
        print("  ⚠ Cap interferent resolt.")
        return obp_results

    results_with_selectivity = []

    for obp in obp_results:
        name = obp["obp_name"]
        print(f"\n  [Selectivitat] {name}")

        pdb_path   = f"docking/{name}.pdb"
        pdbqt_path = f"docking/{name}.pdbqt"

        if not Path(pdb_path).exists():
            print(f"    ⚠ No hi ha estructura, saltant")
            results_with_selectivity.append({
                **obp,
                "si_index":          None,
                "kd_diana_nm":       None,
                "kd_interferent_nm": None,
            })
            continue

        prepare_receptor(pdb_path, pdbqt_path)
        cx, cy, cz = get_center(pdb_path)

        score_diana = obp.get("docking_score_kcal", 0.0)
        kd_diana_nm = score_to_kd_nm(score_diana)

        kd_interferents = []
        for interf in interferent_data:
            interf_name  = interf["name"].replace(" ", "_").replace("/", "-")
            ligand_pdbqt = f"docking/interferents/{interf_name}.pdbqt"
            out_pdbqt    = f"docking/interferents/{name}_{interf_name}_out.pdbqt"

            if not Path(ligand_pdbqt).exists():
                try:
                    prepare_ligand(interf["smiles"], ligand_pdbqt)
                except Exception as e:
                    print(f"    ⚠ Error preparant {interf_name}: {e}")
                    continue

            score_interf = run_vina_score(
                pdbqt_path, ligand_pdbqt, out_pdbqt, cx, cy, cz)

            if score_interf is None:
                print(f"    {interf['name']:20} INVÀLID (saltat)")
                continue

            kd_interf_nm = score_to_kd_nm(score_interf)
            kd_interferents.append(kd_interf_nm)
            print(f"    {interf['name']:20} score: {score_interf:.3f} kcal/mol "
                  f"| Kd pred: {kd_interf_nm:.0f} nM")

        if kd_interferents:
            kd_interf_mig = sum(kd_interferents) / len(kd_interferents)
            si = kd_interf_mig / kd_diana_nm if kd_diana_nm > 0 else 0
            qualitat = ("Excel·lent" if si > 100 else
                        "Bona"       if si > 10  else
                        "Moderada"   if si > 2   else
                        "Baixa")
            print(f"    → SI = {si:.2f} ({qualitat})")
            results_with_selectivity.append({
                **obp,
                "kd_diana_nm":       round(kd_diana_nm, 1),
                "kd_interferent_nm": round(kd_interf_mig, 1),
                "si_index":          round(si, 2),
            })
        else:
            results_with_selectivity.append({
                **obp,
                "kd_diana_nm":       round(kd_diana_nm, 1),
                "kd_interferent_nm": None,
                "si_index":          None,
            })

    return results_with_selectivity


def ranking_selectivity(results):
    print("\n" + "="*65)
    print("RANKING FINAL AMB SELECTIVITAT")
    print("="*65)

    df = pd.DataFrame(results)
    df_valid = df.dropna(subset=["si_index"]).copy()

    if df_valid.empty:
        print("  ⚠ Cap OBP amb dades de selectivitat completes.")
        return

    # Normalitza els tres criteris a [0,1]
    df_valid["score_norm"] = (
        (df_valid["docking_score_kcal"] - df_valid["docking_score_kcal"].max()) /
        (df_valid["docking_score_kcal"].min() - df_valid["docking_score_kcal"].max())
    )
    df_valid["kd_norm"] = (
        (df_valid["kd_nm_experimental"].max() - df_valid["kd_nm_experimental"]) /
        (df_valid["kd_nm_experimental"].max() - df_valid["kd_nm_experimental"].min() + 1e-9)
    )
    df_valid["si_norm"] = (
        (df_valid["si_index"] - df_valid["si_index"].min()) /
        (df_valid["si_index"].max() - df_valid["si_index"].min() + 1e-9)
    )

    df_valid["score_final"] = (
        0.33 * df_valid["score_norm"] +
        0.33 * df_valid["kd_norm"] +
        0.34 * df_valid["si_norm"]
    ).round(3)

    df_valid = df_valid.sort_values("score_final", ascending=False)

    print(f"\n  {'#':<3} {'OBP':<20} {'Kd exp(nM)':<12} {'Docking':<12} {'SI':<8} {'Qualitat':<12} {'Score final'}")
    print("  " + "-"*75)
    for i, row in enumerate(df_valid.itertuples(), 1):
        si = row.si_index
        qualitat = ("Excel·lent ✓✓" if si > 100 else
                    "Bona ✓"        if si > 10  else
                    "Moderada"      if si > 2   else
                    "Baixa ✗")
        print(f"  {i:<3} {row.obp_name:<20} {row.kd_nm_experimental:<12} "
              f"{row.docking_score_kcal:<12.3f} {si:<8.1f} {qualitat:<12} {row.score_final:.3f}")

    df_valid.to_csv("results/ranking_selectivity.csv", index=False)
    print(f"\n  ✓ Ranking guardat a results/ranking_selectivity.csv")
    best = df_valid.iloc[0]
    print(f"\n  → MILLOR CANDIDAT GLOBAL: {best['obp_name']}")
    print(f"    Kd diana: {best['kd_nm_experimental']} nM | SI: {best['si_index']} | Score: {best['score_final']}")