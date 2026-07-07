

import re
import pandas as pd

import common




def load_binding_table():
    cfg = common.load_config()
    df = pd.read_csv(cfg["fitxers"]["binding_csv"])
    cas_col  = df.columns[0]
    name_col = df.columns[1]
    obp_cols = list(df.columns[2:])
    return df, cas_col, name_col, obp_cols


def convert_ki_to_float(raw_value):
    if pd.isna(raw_value):
        return None
    text = str(raw_value).strip().replace('\xa0', '').replace(' ', '')
    if text.startswith('>'):
        try:
            return float(re.sub(r'[^\d.]', '', text)) * 1.1
        except ValueError:
            return None
    try:
        return float(text)
    except ValueError:
        return None


def find_voc_row(df, name_col, cas_col, query):
    q = query.strip().lower()
    mask = df[name_col].astype(str).str.lower().str.contains(q, na=False)
    row = df[mask]
    if row.empty:
        mask2 = df[cas_col].astype(str).str.lower().str.contains(q, na=False)
        row = df[mask2]
    if row.empty:
        return None
    if len(row) > 1:
        print(f"\n  Múltiples coincidències per '{query}':")
        for i, (_, r) in enumerate(row.iterrows(), start=1):
            print(f"    [{i}] {r[name_col]}")
        choice = input("  Tria un número: ").strip()
        try:
            return row.iloc[int(choice) - 1]
        except Exception:
            return row.iloc[0]
    return row.iloc[0]


def find_obp_column(obp_cols, query):
    q = query.strip().lower()
    exact = [c for c in obp_cols if c.lower() == q]
    if exact:
        return exact[0]
    partial = [c for c in obp_cols if q in c.lower()]
    if not partial:
        return None
    if len(partial) > 1:
        print(f"\n  Múltiples OBPs coincideixen amb '{query}':")
        for i, c in enumerate(partial, start=1):
            print(f"    [{i}] {c}")
        choice = input("  Tria un número: ").strip()
        try:
            return partial[int(choice) - 1]
        except Exception:
            return partial[0]
    return partial[0]


def get_experimental_value_nm(row, obp_col):
    """Retorna el Ki experimental ja normalitzat a nM (la taula original és μM), o None."""
    if row is None or obp_col not in row.index:
        return None
    ki_um = convert_ki_to_float(row[obp_col])
    return common.ki_um_to_nm(ki_um)


# Lògica de consulta: per a cada (voc, obp) → experimental, registre o docking


def resolve_pair(df, name_col, cas_col, obp_cols, voc_query, obp_query, smiles_cache):
    """
    Resol una parella (voc_query, obp_query):
      1) Mira la taula experimental (resultat en nM).
      2) Si no hi ha valor → mira el registre de docking ja fet.
      3) Si tampoc hi és → fa docking GNINA per estimar el Kd.
    """
    voc_row = find_voc_row(df, name_col, cas_col, voc_query)
    voc_name = voc_row[name_col] if voc_row is not None else voc_query
    if voc_row is None:
        print(f"   VOC '{voc_query}' no trobat a la taula. Es tractarà com a nou (caldrà SMILES).")

    obp_col = find_obp_column(obp_cols, obp_query)
    if obp_col is None:
        print(f"   OBP '{obp_query}' no trobat. Saltant.")
        return None

    kd_exp_nm = get_experimental_value_nm(voc_row, obp_col) if voc_row is not None else None
    if kd_exp_nm is not None:
        return {
            "obp_name": obp_col,
            "voc_name": voc_name,
            "kd_nm": kd_exp_nm,
            "font": "experimental",
        }

    ja_fet = common.docking_already_done(obp_col, voc_name, mode="cnn")
    if ja_fet is not None:
        print(f"\n   Ja hi ha un docking previ registrat per {obp_col} × {voc_name} (reutilitzant)")
        kd = ja_fet.get("kd_nm_predit")
        vina = ja_fet.get("vina_score")
        return {
            "obp_name": obp_col,
            "voc_name": voc_name,
            "kd_nm": float(kd) if pd.notna(kd) and kd != "" else None,
            "vina_score": float(vina) if pd.notna(vina) and vina != "" else None,
            "font": "docking_previ",
        }

    print(f"\n  Sense dada experimental per {obp_col} × {voc_name} → docking GNINA (CNN)")
    result = common.dock_pair(obp_col, voc_name, smiles_cache=smiles_cache, mode="cnn")
    if result is None:
        return None
    return {
        "obp_name": result["obp_name"],
        "voc_name": result["voc_name"],
        "kd_nm": result["kd_nm_predit"],
        "vina_score": result["vina_score"],
        "font": result["font"],
    }


def print_result_row(res):
    if res is None:
        return
    kd_str = f"{res['kd_nm']:.1f} nM" if res.get("kd_nm") is not None else "N/A"
    if res["font"] == "experimental":
        origen = "[iobpdb, experimental]"
    elif res["font"] == "docking_previ":
        origen = "[docking previ, reutilitzat]"
    else:
        origen = "[docking GNINA nou]"
    extra = f"  Vina={res['vina_score']:.2f} kcal/mol" if res.get("vina_score") is not None else ""
    print(f"  {res['obp_name']:<20} × {res['voc_name']:<30}  Kd = {kd_str}{extra}  {origen}")


def save_results(results, out_path):
    rows = []
    for r in results:
        if r is None:
            continue
        rows.append({
            "OBP": r["obp_name"],
            "VOC": r["voc_name"],
            "Kd_nM": r.get("kd_nm"),
            "Vina_score": r.get("vina_score"),
            "Font": r["font"],
        })
    if rows:
        pd.DataFrame(rows).to_csv(out_path, index=False)
        print(f"\n   Resultats desats a {out_path}  (totes les afinitats en nM)")



# MAIN


def ask_list(prompt_text):
    raw = input(prompt_text).strip()
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def comptar_pendents(df, name_col, cas_col, obp_cols, parelles):
    """Compta quantes parelles (voc, obp) NO tenen dada experimental ni docking previ."""
    pendents = 0
    for voc_q, obp_q in parelles:
        voc_row = find_voc_row(df, name_col, cas_col, voc_q)
        obp_col = find_obp_column(obp_cols, obp_q)
        if obp_col is None:
            continue
        kd_exp = get_experimental_value_nm(voc_row, obp_col) if voc_row is not None else None
        if kd_exp is not None:
            continue
        voc_name = voc_row[name_col] if voc_row is not None else voc_q
        if common.docking_already_done(obp_col, voc_name, mode="cnn") is not None:
            continue
        pendents += 1
    return pendents


def main():
    common.ensure_dirs()
    print("=" * 60)
    print("  CONSULTA RÀPIDA D'AFINITAT — VOC ↔ OBP  (iGEM URV 2025)")
    print("=" * 60)
    print("""
  Tries:
    [1] Un VOC  contra un o més OBPs
    [2] Un OBP  contra un o més VOCs
    [3] Un VOC  contra un OBP (consulta directa, 1 a 1)
""")
    mode = input("  Tria una opció [1/2/3]: ").strip()

    df, cas_col, name_col, obp_cols = load_binding_table()
    smiles_cache = common.load_smiles_cache()
    results = []

    if mode == "1":
        voc_query = input("\n  Nom o CAS del VOC: ").strip()
        obp_queries = ask_list("  Noms dels OBPs separats per comes (buit = TOTS els OBPs de la taula): ")
        if not obp_queries:
            obp_queries = obp_cols
        parelles = [(voc_query, o) for o in obp_queries]

    elif mode == "2":
        obp_query = input("\n  Nom de l'OBP: ").strip()
        voc_queries = ask_list("  Noms o CAS dels VOCs separats per comes (buit = TOTS els VOCs de la taula): ")
        if not voc_queries:
            voc_queries = df[name_col].astype(str).tolist()
        parelles = [(v, obp_query) for v in voc_queries]

    elif mode == "3":
        voc_query = input("\n  Nom o CAS del VOC: ").strip()
        obp_query = input("  Nom de l'OBP: ").strip()
        parelles = [(voc_query, obp_query)]

    else:
        print("  Opció no vàlida.")
        return

    n_pendents = comptar_pendents(df, name_col, cas_col, obp_cols, parelles)
    if n_pendents > 0 and not common.confirmar_docking_massiu(n_pendents):
        print("  Operació cancel·lada.")
        return

    print(f"\n  Consultant {len(parelles)} parella(es) (de les quals {n_pendents} necessiten docking nou)...\n")

    for i, (voc_q, obp_q) in enumerate(parelles, start=1):
        common.print_progress(i, len(parelles), label=f"{obp_q} × {voc_q}")
        res = resolve_pair(df, name_col, cas_col, obp_cols, voc_q, obp_q, smiles_cache)
        results.append(res)

    print("\n" + "=" * 60)
    print("  RESULTATS  (totes les afinitats normalitzades a nM)")
    print("=" * 60)
    for r in results:
        print_result_row(r)

    cfg = common.load_config()
    out_path = f"{cfg['fitxers']['carpeta_results']}/consulta_resultats.csv"
    save_results(results, out_path)


if __name__ == "__main__":
    main()