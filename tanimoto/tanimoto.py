

import os
import json
import time
import requests
import pandas as pd
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   


BINDING_FILE     = "Compound_OBP_binding.csv"
CACHE_FILE       = "docking/voc_smiles_cache.json"
OUT_MATRIX_FILE  = "docking/tanimoto_matrix.npy"
OUT_INDEX_FILE   = "docking/voc_index.json"
OUT_MISSING_FILE = "docking/voc_sense_smiles.csv"

MORGAN_RADIUS    = 2
MORGAN_NBITS     = 2048

os.makedirs("docking", exist_ok=True)


PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _pubchem_resolve(identifier, max_retries=3, wait=2):

    for intent in range(max_retries):
        try:
            url_cid = f"{PUBCHEM_URL}/compound/name/{requests.utils.quote(str(identifier))}/cids/JSON"
            r = requests.get(url_cid, timeout=10)

            if r.status_code == 404:
                return None  # nom no existeix a PubChem, no té sentit reintentar
            r.raise_for_status()
            cid = r.json()["IdentifierList"]["CID"][0]

            url_props = (f"{PUBCHEM_URL}/compound/cid/{cid}/property/"
                         f"IsomericSMILES,CanonicalSMILES/JSON")
            r2 = requests.get(url_props, timeout=10)
            r2.raise_for_status()
            props = r2.json()["PropertyTable"]["Properties"][0]

            smiles = (props.get("IsomericSMILES")
                      or props.get("CanonicalSMILES")
                      or props.get("SMILES")
                      or props.get("ConnectivitySMILES"))
            return smiles

        except (requests.exceptions.RequestException, KeyError, IndexError, ValueError):
            time.sleep(wait * (intent + 1))  # backoff progressiu abans de reintentar
            continue

    return None


def resolve_smiles_for_compound(cas, compound_name, cache):

    cas_key = str(cas).strip()

    if cas_key in cache and cache[cas_key]:
        return cache[cas_key]

    smiles = None


    sinonims = [s.strip() for s in str(compound_name).split("/") if s.strip()]

    for nom in sinonims:
        smiles = _pubchem_resolve(nom)
        if smiles:
            break

    if not smiles:
        smiles = _pubchem_resolve(cas_key)

    if smiles:
        cache[cas_key] = smiles

    return smiles


#  FINGERPRINTS I TANIMOTO 

def build_fingerprint(smiles):

    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)


def build_tanimoto_matrix(fingerprints):

    n = len(fingerprints)
    matrix = np.full((n, n), np.nan, dtype=np.float32)

    valid_idx = [i for i, fp in enumerate(fingerprints) if fp is not None]
    print(f"\n  Calculant Tanimoto per a {len(valid_idx)}/{n} VOCs amb SMILES vàlid...")

    for count, i in enumerate(valid_idx, 1):
        fp_i = fingerprints[i]
        sims = DataStructs.BulkTanimotoSimilarity(
            fp_i, [fingerprints[j] for j in valid_idx]
        )
        for sim, j in zip(sims, valid_idx):
            matrix[i, j] = sim
        if count % 50 == 0 or count == len(valid_idx):
            print(f"    {count}/{len(valid_idx)} VOCs processats...")

    return matrix




def main():
    print("="*70)
    print("  GENERACIÓ MATRIU TANIMOTO VOC×VOC")
    print("="*70)

    if not os.path.isfile(BINDING_FILE):
        print(f"   ERROR: no es troba '{BINDING_FILE}' en aquesta carpeta.")
        return

    print(f"\n  Llegint {BINDING_FILE}...")
    df = pd.read_csv(BINDING_FILE)
    cas_col, name_col = df.columns[0], df.columns[1]
    print(f"  {len(df)} VOCs trobats")


    cache = {}
    if os.path.isfile(CACHE_FILE):
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            cache = json.load(f)
        print(f"   Caché de SMILES carregada: {len(cache)} entrades existents")
    else:
        print(f"   No hi ha caché prèvia, es crearà '{CACHE_FILE}'")

    print(f"\n  Resolent SMILES (caché  PubChem CID→propietats, per sinònim → per CAS)...")
    cas_list, name_list, smiles_list = [], [], []
    sense_smiles = []

    for idx, row in df.iterrows():
        cas  = row[cas_col]
        name = row[name_col]
        smiles = resolve_smiles_for_compound(cas, name, cache)

        cas_list.append(str(cas).strip())
        name_list.append(name)
        smiles_list.append(smiles)

        status = "ok" if smiles else "no"
        if (idx + 1) % 25 == 0 or smiles is None:
            print(f"    [{idx+1}/{len(df)}] {status}  {str(name)[:50]}")

        if not smiles:
            sense_smiles.append({"CAS-number": cas, "Compound name": name})

    # Desa caché actualitzada (perquè properes execucions siguin instantànies)
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)
    print(f"\n   Caché de SMILES desada a {CACHE_FILE} ({len(cache)} entrades)")

    if sense_smiles:
        pd.DataFrame(sense_smiles).to_csv(OUT_MISSING_FILE, index=False)
        print(f"   {len(sense_smiles)} VOCs sense SMILES resolt → {OUT_MISSING_FILE}")
        print(f"    (no participaran a la matriu Tanimoto; revisa'ls a mà si vols)")

    print(f"\n  Generant fingerprints Morgan (radius={MORGAN_RADIUS}, nBits={MORGAN_NBITS})...")
    fingerprints = [build_fingerprint(s) for s in smiles_list]
    n_ok = sum(1 for fp in fingerprints if fp is not None)
    print(f"  → {n_ok}/{len(fingerprints)} fingerprints generats correctament")


    matrix = build_tanimoto_matrix(fingerprints)

    np.save(OUT_MATRIX_FILE, matrix)
    print(f"\n   Matriu Tanimoto ({matrix.shape[0]}×{matrix.shape[1]}) desada a {OUT_MATRIX_FILE}")

    with open(OUT_INDEX_FILE, "w", encoding="utf-8") as f:
        json.dump(cas_list, f, indent=2, ensure_ascii=False)
    print(f"   Índex (ordre CAS) desat a {OUT_INDEX_FILE}")

    print("\n" + "="*70)
    print("  RESUM")
    print("="*70)
    print(f"  VOCs totals             : {len(df)}")
    print(f"  Amb SMILES resolt       : {n_ok}")
    print(f"  Sense SMILES (revisar)  : {len(sense_smiles)}")
    print(f"\n  Per carregar la matriu després:")
    print(f"    import numpy as np, json")
    print(f"    matrix = np.load('{OUT_MATRIX_FILE}')")
    print(f"    cas_index = json.load(open('{OUT_INDEX_FILE}'))")
    print(f"    i = cas_index.index('CAS_DEL_VOC')")
    print(f"    fila_similituds = matrix[i]   # Tanimoto contra tots els altres VOCs")


if __name__ == "__main__":
    main()