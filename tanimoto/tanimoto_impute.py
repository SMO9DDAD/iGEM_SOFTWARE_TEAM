"""
tanimoto_impute.py
─
Mòdul base per al protocol Tanimoto-Ki (read-across / vNN).
Implementa exactament les funcions del Tanimoto_Ki_Protocol.docx:
  - impute_ki_tanimoto()  → predicció puntual (Ki_pred + Ki_lower conservador)
  - run_loo()             → validació Leave-One-Out (RMSE)
  - fill_gaps()           → omple tots els buits de binding_table

Connecta amb la matriu Tanimoto que ja tens generada
(docking/tanimoto_matrix.npy + docking/voc_index.json), fent l'enllaç
CAS ↔ posició a la matriu ↔ fila de binding_table.

No depèn de main_nou.py: és un mòdul autònom que main_nou.py importarà.

"""

import json
import numpy as np
import pandas as pd




def load_tanimoto_data(matrix_path="docking/tanimoto_matrix.npy",
                        index_path="docking/voc_index.json"):
    """
    Carrega la matriu Tanimoto i el seu índex de CAS.
    Retorna (T_matrix, cas_to_pos) on cas_to_pos és un dict CAS -> posició
    (fila/columna) dins de T_matrix.
    """
    T_matrix = np.load(matrix_path)
    cas_index = json.load(open(index_path, encoding="utf-8"))
    cas_to_pos = {cas: i for i, cas in enumerate(cas_index)}
    return T_matrix, cas_to_pos


def build_binding_idx_to_matrix_pos(binding_table, cas_col, cas_to_pos):

    mapping = {}
    for row_idx, cas in binding_table[cas_col].astype(str).str.strip().items():
        if cas in cas_to_pos:
            mapping[row_idx] = cas_to_pos[cas]
    return mapping


#  PREDICCIÓ PUNTUAL (read-across) 

def impute_ki_tanimoto(x_idx, obp_col, binding_table, T_matrix, idx_to_pos,
                        t_min=0.35, k=5, p=2.0, z=1.0, exclude_idx=None):

    if x_idx not in idx_to_pos:
        return None  

    pos_x = idx_to_pos[x_idx]
    col = binding_table[obp_col]

    donors = []
    for cand_idx, pos_cand in idx_to_pos.items():
        if cand_idx == x_idx or cand_idx == exclude_idx:
            continue
        ki_cand = col.get(cand_idx, np.nan)
        if pd.isna(ki_cand) or ki_cand <= 0:
            continue
        t = float(T_matrix[pos_x, pos_cand])
        if np.isnan(t) or t < t_min:
            continue
        donors.append((cand_idx, ki_cand, t))

    if not donors:
        return None

    donors.sort(key=lambda d: -d[2])   # ordenar per Tanimoto descendent
    donors = donors[:k]                # top-k

    ki_vals = np.array([d[1] for d in donors], dtype=float)
    T_vals  = np.array([d[2] for d in donors], dtype=float)

    pKi = -np.log10(ki_vals * 1e-6)    # μM → M → pKi
    w   = T_vals ** p

    pKi_hat = float((w * pKi).sum() / w.sum())
    sigma   = float(np.sqrt((w * (pKi - pKi_hat) ** 2).sum() / w.sum()))

    return {
        "Ki_pred_uM":  10 ** (-pKi_hat) * 1e6,
        "Ki_lower_uM": 10 ** (-(pKi_hat + z * sigma)) * 1e6,
        "T_max":       float(T_vals.max()),
        "sigma_log":   sigma,
        "n_donors":    len(donors),
    }


#  VALIDACIÓ LEAVE-ONE-OUT (LOO) 

def run_loo(binding_table, T_matrix, idx_to_pos, obp_name_list,
            t_min=0.35, k=5, p=2.0):

    errors = []
    total_cols = len(obp_name_list)

    for col_i, obp_col in enumerate(obp_name_list, start=1):
        measured = binding_table[obp_col].dropna()
        measured = measured[measured > 0]

        for target_idx in measured.index:
            real_ki = measured[target_idx]
            result = impute_ki_tanimoto(
                target_idx, obp_col, binding_table, T_matrix, idx_to_pos,
                t_min=t_min, k=k, p=p, exclude_idx=target_idx,
            )
            if result is None:
                continue
            pred_pki = -np.log10(result["Ki_pred_uM"] * 1e-6)
            real_pki = -np.log10(real_ki * 1e-6)
            errors.append({
                "obp": obp_col,
                "voc_idx": target_idx,
                "error": pred_pki - real_pki,
                "T_max": result["T_max"],
                "n_donors": result["n_donors"],
            })

        if col_i % 50 == 0 or col_i == total_cols:
            print(f"    LOO: {col_i}/{total_cols} columnes OBP processades "
                  f"({len(errors)} prediccions acumulades)...")

    if not errors:
        print("  LOO no ha pogut validar cap cel·la (revisa t_min/k).")
        return None, pd.DataFrame()

    df_err = pd.DataFrame(errors)
    rmse = float(np.sqrt((df_err["error"] ** 2).mean()))
    return rmse, df_err


def rmse_by_tmax_bin(df_err, bins=(0.35, 0.5, 0.7, 1.01)):

    labels = [f"{bins[i]:.2f}–{bins[i+1]:.2f}" for i in range(len(bins) - 1)]
    df_err = df_err.copy()
    df_err["bin"] = pd.cut(df_err["T_max"], bins=bins, labels=labels, right=False)
    resumen = df_err.groupby("bin")["error"].apply(
        lambda e: np.sqrt((e ** 2).mean())
    ).rename("RMSE")
    counts = df_err.groupby("bin").size().rename("n_cel·les")
    return pd.concat([resumen, counts], axis=1)


#  OMPLIR TOTS ELS BUITS REALS

def fill_gaps(binding_table, T_matrix, idx_to_pos, obp_name_list,
              t_min=0.35, k=5, p=2.0, z=1.0):
    """
    Omple TOTS els buits reals de binding_table (no és LOO, és l'aplicació
    final). Treballa sobre una còpia: mai sobreescriu cel·les mesurades.

    Retorna:
      binding_imputed : còpia de binding_table amb Ki_lower_uM a les cel·les
                         imputades (valor conservador, per a selectivitat)
      binding_pred     : còpia paral·lela amb Ki_pred_uM (valor central,
                         per a mostrar/reportar)
      diag_table       : taula llarga amb una fila per cada cel·la imputada
                         (CAS, OBP, Ki_pred, Ki_lower, T_max, sigma, n_donors,
                         decision: 'confident' / 'exploratory' / 'no_impute')
    """
    binding_imputed = binding_table.copy()
    binding_pred     = binding_table.copy()
    diag_rows = []

    total_cols = len(obp_name_list)
    for col_i, obp_col in enumerate(obp_name_list, start=1):
        missing_idx = binding_table[binding_table[obp_col].isna()].index

        for x_idx in missing_idx:
            result = impute_ki_tanimoto(
                x_idx, obp_col, binding_table, T_matrix, idx_to_pos,
                t_min=t_min, k=k, p=p, z=z,
            )
            if result is None:
                continue

            # Decisió segons la taula del protocol (T_max i sigma)
            if result["T_max"] >= 0.5 and result["sigma_log"] < 0.3:
                decision = "confident"
            elif result["T_max"] >= t_min:
                decision = "exploratory"
            else:
                decision = "no_impute"

            binding_imputed.at[x_idx, obp_col] = result["Ki_lower_uM"]
            binding_pred.at[x_idx, obp_col]     = result["Ki_pred_uM"]

            diag_rows.append({
                "row_idx":     x_idx,
                "OBP":         obp_col,
                "Ki_pred_uM":  result["Ki_pred_uM"],
                "Ki_lower_uM": result["Ki_lower_uM"],
                "T_max":       result["T_max"],
                "sigma_log":   result["sigma_log"],
                "n_donors":    result["n_donors"],
                "decision":    decision,
            })

        if col_i % 50 == 0 or col_i == total_cols:
            print(f"    Omplint: {col_i}/{total_cols} columnes OBP processades "
                  f"({len(diag_rows)} cel·les imputades fins ara)...")

    diag_table = pd.DataFrame(diag_rows)
    return binding_imputed, binding_pred, diag_table