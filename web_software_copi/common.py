# -*- coding: utf-8 -*-
"""
common.py — Funcions compartides del pipeline OBP/VOC (iGEM URV 2025)

Centralitza tot el que abans estava duplicat a main_nou.py, main3_gnina.py
i main_consulta.py:
  - Lectura de config.yaml
  - Cache persistent de SMILES (evita repetir consultes a PubChem)
  - Obtenció/predicció d'estructures d'OBP (iobpdb → ESMFold)
  - Preparació de lligand/receptor (RDKit + Meeko + OpenBabel)
  - Detecció de pocket amb P2Rank
  - Docking amb GNINA (Vina pur o CNN rescoring) i conversió CNN→Kd
  - Registre persistent de tots els dockings fets (docking/registre.csv)

Tots els valors numèrics (exhaustiveness, mides de caixa, etc.) venen de
config.yaml — no hi ha cap número "màgic" repetit pel codi.
"""

import os
import re
import csv
import json
import time
import subprocess
from datetime import datetime
from pathlib import Path

import requests
import pandas as pd
import yaml

_BASE_DIR = Path(__file__).parent
_CONFIG_PATH = _BASE_DIR / "config.yaml"

_config_cache = None


# ═══════════════════════════════════════════════════════════════════════════
# Configuració
# ═══════════════════════════════════════════════════════════════════════════

def load_config():
    """Carrega config.yaml (amb cache en memòria per no rellegir-lo sempre)."""
    global _config_cache
    if _config_cache is not None:
        return _config_cache
    if not _CONFIG_PATH.is_file():
        raise FileNotFoundError(
            f"No es troba config.yaml a {_CONFIG_PATH}. "
            "Posa'l a la mateixa carpeta que els scripts del pipeline."
        )
    with open(_CONFIG_PATH, "r", encoding="utf-8") as f:
        _config_cache = yaml.safe_load(f)
    return _config_cache


def ensure_dirs():
    cfg = load_config()
    os.makedirs(cfg["fitxers"]["carpeta_results"], exist_ok=True)
    os.makedirs(cfg["fitxers"]["carpeta_docking"], exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Cache persistent de SMILES
# ═══════════════════════════════════════════════════════════════════════════

def _cache_path():
    cfg = load_config()
    return Path(cfg["fitxers"]["cache_smiles"])


def load_smiles_cache():
    path = _cache_path()
    if path.is_file():
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def save_smiles_cache(cache):
    path = _cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def resolve_smiles(voc_name, cache=None):
    """
    Obté el SMILES d'un VOC. Ordre: cache local → PubChem → entrada manual.
    Si es passa `cache` (dict), es modifica in-place i es desa a disc.
    """
    own_cache = cache if cache is not None else load_smiles_cache()

    if voc_name in own_cache and own_cache[voc_name]:
        return own_cache[voc_name]

    cfg = load_config()
    smiles = None
    try:
        url = cfg["urls"]["pubchem_smiles"].format(voc_name)
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            smiles = r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except Exception as e:
        print(f"  Avís PubChem: {e}")

    if not smiles:
        smiles = input(f"  Introdueix el SMILES manualment per a '{voc_name}': ").strip()

    own_cache[voc_name] = smiles
    save_smiles_cache(own_cache)
    if cache is not None:
        cache.update(own_cache)
    return smiles


# ═══════════════════════════════════════════════════════════════════════════
# Estructures d'OBP (iobpdb → ESMFold)
# ═══════════════════════════════════════════════════════════════════════════

def download_iobpdb_structure(obp_name, out_pdb):
    cfg = load_config()
    try:
        url = cfg["urls"]["iobpdb_pdb"].format(obp_name)
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            print("    ✓ Estructura descarregada d'iobpdb")
            return True
        print(f"    ⚠ No trobat a iobpdb (HTTP {r.status_code})")
        return False
    except Exception as e:
        print(f"    ⚠ Error descarregant: {e}")
        return False


def get_uniprot_id(obp_name):
    cfg = load_config()
    try:
        url = cfg["urls"]["uniprot_search"].format(obp_name)
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                return results[0]["primaryAccession"]
    except Exception:
        pass
    return None


def get_sequence(uniprot_id):
    cfg = load_config()
    url = cfg["urls"]["uniprot_fasta"].format(uniprot_id)
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        return None
    lines = r.text.strip().split("\n")
    return "".join(lines[1:])


def predict_structure_esm(seq, out_pdb):
    cfg = load_config()
    max_len = cfg["docking"]["esmfold_max_len"]
    url = cfg["urls"]["esmfold"]
    if len(seq) > max_len:
        seq = seq[20:]
    for intent in range(3):
        r = requests.post(url, data=seq, timeout=180)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            return True
        print(f"    ESMFold intent {intent + 1} error {r.status_code}, esperant...")
        time.sleep(10)
    return False


def get_or_build_structure(obp_name):
    """Obté o construeix l'estructura PDB d'un OBP. Ordre: cache local → iobpdb → ESMFold."""
    cfg = load_config()
    pdb_path = os.path.join(cfg["fitxers"]["carpeta_docking"], f"{obp_name}.pdb")
    if os.path.exists(pdb_path):
        return pdb_path

    print("    Buscant estructura a iobpdb...")
    if download_iobpdb_structure(obp_name, pdb_path):
        return pdb_path

    print("    Intentant ESMFold...")
    uid = get_uniprot_id(obp_name)
    if not uid:
        print("    ⚠ UniProt ID no trobat")
        return None
    seq = get_sequence(uid)
    if not seq:
        print("    ⚠ Seqüència no trobada")
        return None
    print(f"    Predint estructura ({len(seq)} aa)...")
    if not predict_structure_esm(seq, pdb_path):
        print("    ⚠ ESMFold ha fallat")
        return None
    print("    ✓ Estructura predita per ESMFold")
    return pdb_path


# ═══════════════════════════════════════════════════════════════════════════
# Preparació de lligand / receptor
# ═══════════════════════════════════════════════════════════════════════════

def prepare_ligand_pdbqt(smiles, out_pdbqt):
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


def prepare_receptor_pdbqt(pdb_path, pdbqt_path):
    os.system(f"obabel {pdb_path} -O {pdbqt_path} -xr 2>/dev/null")


def get_center(pdb_path):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except Exception:
                    pass
    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)


def detect_pocket_p2rank(pdb_path):
    cfg = load_config()
    p2rank_exe = str(_BASE_DIR / cfg["executables"]["p2rank_exe"])
    if not os.path.exists(p2rank_exe):
        return None

    pdb_abs  = str(Path(pdb_path).resolve())
    out_dir  = str(Path(pdb_path + "_p2rank").resolve())
    p2rank_d = str(Path(p2rank_exe).parent.resolve())
    os.makedirs(out_dir, exist_ok=True)

    try:
        cmd = ["bash", p2rank_exe, "predict", "-f", pdb_abs, "-o", out_dir, "-threads", "1"]
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=cfg["docking"]["timeout_p2rank_s"], cwd=p2rank_d
        )
        if result.returncode != 0:
            return None
        candidates = list(Path(out_dir).rglob("*predictions.csv"))
        if not candidates:
            return None
        df_pockets = pd.read_csv(candidates[0])
        df_pockets.columns = df_pockets.columns.str.strip()
        if df_pockets.empty:
            return None
        best = df_pockets.iloc[0]
        cx, cy, cz = float(best["center_x"]), float(best["center_y"]), float(best["center_z"])
        p2rank_score = float(best.get("score", 0))
        print(f"    [P2Rank] Pocket: ({cx:.1f}, {cy:.1f}, {cz:.1f})  score={p2rank_score:.2f}", end="")
        if p2rank_score < cfg["docking"]["p2rank_min_score"]:
            print(f"  ⚠ score < {cfg['docking']['p2rank_min_score']}, usant centroide")
            return None
        print()
        return cx, cy, cz
    except Exception as e:
        print(f"    [P2Rank] Error: {e}")
        return None


def get_docking_box(pdb_path):
    """Retorna (cx, cy, cz, box_size) triant pocket P2Rank o centroide proteic."""
    cfg = load_config()
    pocket = detect_pocket_p2rank(pdb_path)
    if pocket:
        cx, cy, cz = pocket
        box_size = cfg["docking"]["box_with_pocket"]
        print(f"    Grid al pocket detectat ({box_size} Å)")
    else:
        cx, cy, cz = get_center(pdb_path)
        box_size = cfg["docking"]["box_fallback"]
        print(f"    Grid al centroide proteic ({box_size} Å)")
    return cx, cy, cz, box_size


# ═══════════════════════════════════════════════════════════════════════════
# Docking GNINA — Vina pur i CNN rescoring + conversió a Kd (nM)
# ═══════════════════════════════════════════════════════════════════════════

def cnn_to_kd_nm(cnn_affinity_pkd):
    """Converteix CNN affinity (pKd = -log10(Kd[M])) a Kd en nM."""
    if cnn_affinity_pkd is None:
        return None
    kd_m = 10 ** (-cnn_affinity_pkd)
    return kd_m * 1e9


def ki_um_to_nm(ki_um):
    """Normalitza un Ki/Kd experimental, donat en μM, a nM (unitat única de sortida)."""
    if ki_um is None:
        return None
    return ki_um * 1000.0


def run_gnina(receptor, ligand, out, cx, cy, cz, box_size, mode="cnn"):
    """
    Docking amb GNINA.
      mode="vina" → scoring Vina pur (--cnn_scoring none). Retorna (vina_score, None).
      mode="cnn"  → CNN rescoring (--cnn_scoring rescore). Retorna (vina_score, cnn_pkd).
    """
    cfg = load_config()
    gnina_exe = str(_BASE_DIR / cfg["executables"]["gnina_exe"])
    base_cmd = [
        gnina_exe, "--receptor", receptor, "--ligand", ligand, "--out", out,
        "--center_x", str(round(cx, 2)), "--center_y", str(round(cy, 2)),
        "--center_z", str(round(cz, 2)),
        "--size_x", str(box_size), "--size_y", str(box_size), "--size_z", str(box_size),
        "--exhaustiveness", str(cfg["docking"]["exhaustiveness"]),
        "--num_modes", str(cfg["docking"]["num_modes"]),
        "--no_gpu",
    ]
    if mode == "vina":
        cmd = base_cmd + ["--cnn_scoring", "none", "--scoring", "vina"]
    else:
        cmd = base_cmd + ["--cnn_scoring", "rescore", "--cnn", "crossdock_default2018"]

    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=cfg["docking"]["timeout_gnina_s"]
    )

    best_vina, best_cnn = None, None
    for line in result.stdout.split("\n"):
        line = line.strip()
        if re.match(r"^1\s", line):
            parts = line.split()
            try:
                best_vina = float(parts[1])
                if mode == "cnn":
                    best_cnn = float(parts[4])
            except Exception:
                pass
            break

    if best_vina is None and result.returncode != 0:
        print(f"    GNINA error: {result.stderr[-200:]}")

    return best_vina, best_cnn


def dock_pair(obp_name, voc_name, smiles_cache=None, mode="cnn"):
    """
    Fa docking complet (estructura + lligand + grid + GNINA) entre un OBP i un VOC.
    mode="cnn"  → retorna també kd_nm_predit (per a consultes puntuals).
    mode="vina" → només vina_score (ranking ràpid, com main_nou.py original).
    Retorna un dict o None si falla. Registra automàticament a docking/registre.csv.
    """
    ensure_dirs()
    cfg = load_config()
    docking_dir = cfg["fitxers"]["carpeta_docking"]

    print(f"\n  [Docking] {obp_name}  ×  {voc_name}")
    pdb_path = get_or_build_structure(obp_name)
    if not pdb_path:
        log_docking(obp_name, voc_name, mode, None, None, None, ok=False)
        return None

    pdbqt_path = os.path.join(docking_dir, f"{obp_name}.pdbqt")
    if not os.path.exists(pdbqt_path):
        prepare_receptor_pdbqt(pdb_path, pdbqt_path)

    voc_slug = re.sub(r"[^\w]+", "_", voc_name)
    ligand_pdbqt = os.path.join(docking_dir, f"ligand_{voc_slug}.pdbqt")
    if not os.path.exists(ligand_pdbqt):
        smiles = resolve_smiles(voc_name, cache=smiles_cache)
        print("    Preparant lligand...")
        prepare_ligand_pdbqt(smiles, ligand_pdbqt)
        print("    ✓ Lligand preparat")

    out_path = os.path.join(docking_dir, f"{obp_name}_{voc_slug}_{mode}_out.pdbqt")
    cx, cy, cz, box_size = get_docking_box(pdb_path)

    vina_score, cnn_pkd = run_gnina(pdbqt_path, ligand_pdbqt, out_path, cx, cy, cz, box_size, mode=mode)

    if vina_score is None:
        print("    ✗ Docking ha fallat")
        log_docking(obp_name, voc_name, mode, None, None, None, ok=False)
        return None

    kd_nm_pred = cnn_to_kd_nm(cnn_pkd) if mode == "cnn" else None

    print(f"    ✓ Vina score: {vina_score:.3f} kcal/mol")
    if kd_nm_pred is not None:
        print(f"    ✓ CNN affinity (pKd): {cnn_pkd:.3f}  →  Kd predit: {kd_nm_pred:.1f} nM")

    log_docking(obp_name, voc_name, mode, vina_score, cnn_pkd, kd_nm_pred, ok=True)

    return {
        "obp_name": obp_name,
        "voc_name": voc_name,
        "vina_score": vina_score,
        "cnn_pkd": cnn_pkd,
        "kd_nm_predit": kd_nm_pred,
        "font": "docking_gnina_cnn" if mode == "cnn" else "docking_gnina_vina",
    }


# ═══════════════════════════════════════════════════════════════════════════
# Registre persistent de docking (docking/registre.csv)
# ═══════════════════════════════════════════════════════════════════════════

_REGISTRE_HEADER = [
    "timestamp", "obp", "voc", "mode", "vina_score",
    "cnn_pkd", "kd_nm_predit", "ok",
]


def log_docking(obp, voc, mode, vina_score, cnn_pkd, kd_nm_predit, ok):
    cfg = load_config()
    path = Path(cfg["fitxers"]["registre_docking"])
    path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = path.is_file()
    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(_REGISTRE_HEADER)
        writer.writerow([
            datetime.now().isoformat(timespec="seconds"),
            obp, voc, mode,
            f"{vina_score:.3f}" if vina_score is not None else "",
            f"{cnn_pkd:.3f}" if cnn_pkd is not None else "",
            f"{kd_nm_predit:.1f}" if kd_nm_predit is not None else "",
            "1" if ok else "0",
        ])


def docking_already_done(obp, voc, mode):
    """Comprova al registre si ja s'ha fet aquest docking amb èxit (evita repetir-lo)."""
    cfg = load_config()
    path = Path(cfg["fitxers"]["registre_docking"])
    if not path.is_file():
        return None
    try:
        df = pd.read_csv(path)
    except Exception:
        return None
    match = df[(df["obp"] == obp) & (df["voc"] == voc) & (df["mode"] == mode) & (df["ok"] == 1)]
    if match.empty:
        return None
    return match.iloc[-1].to_dict()


# ═══════════════════════════════════════════════════════════════════════════
# Confirmació + estimació de temps abans de docking massiu
# ═══════════════════════════════════════════════════════════════════════════

def confirmar_docking_massiu(n_dockings):
    """
    Si n_dockings supera el llindar de config.yaml, mostra una estimació de
    temps i demana confirmació explícita. Retorna True/False.
    """
    cfg = load_config()
    llindar = cfg["seguretat"]["avis_docking_massiu"]
    minuts_per_docking = cfg["seguretat"]["minuts_estimats_per_docking"]

    if n_dockings <= llindar:
        return True

    minuts_total = n_dockings * minuts_per_docking
    hores = minuts_total / 60
    print(f"\n  ⚠ Atenció: s'han de fer {n_dockings} dockings.")
    print(f"    Estimació: ~{minuts_per_docking} min/docking → ~{minuts_total} min (~{hores:.1f} h) en total.")
    resp = input("  Vols continuar? [s/N]: ").strip().lower()
    return resp in ("s", "si", "sí", "y", "yes")


def print_progress(current, total, label=""):
    pct = (current / total) * 100 if total else 0
    print(f"\n  [Progrés: {current}/{total} — {pct:.0f}%] {label}")