
import os
import numpy as np
import pandas as pd

import tanimoto_impute as ti

BINDING_FILE = "Compound_OBP_binding.csv"

# Graella de paràmetres a provar (pots ampliar-la si vols)
T_MIN_VALUES = [0.30, 0.35, 0.45, 0.55]
K_VALUES     = [3, 5, 8]
P_VALUES     = [1.0, 2.0]

os.makedirs("results", exist_ok=True)


def main():
    print("=" * 70)
    print("  CALIBRATGE TANIMOTO-Ki (Leave-One-Out)")
    print("=" * 70)

    print(f"\n  Llegint {BINDING_FILE}...")
    df = pd.read_csv(BINDING_FILE)
    cas_col = df.columns[0]
    obp_name_list = list(df.columns[2:])

    # Convertim Ki igual que main_nou.py (valors censurats ">X" → X*1.1)
    def convert_ki(raw):
        if pd.isna(raw):
            return np.nan
        text = str(raw).strip().replace('\xa0', '').replace(' ', '')
        if text.startswith('>'):
            try:
                import re
                return float(re.sub(r'[^\d.]', '', text)) * 1.1
            except ValueError:
                return np.nan
        try:
            return float(text)
        except ValueError:
            return np.nan

    binding_table = df.copy()
    for col in obp_name_list:
        binding_table[col] = binding_table[col].apply(convert_ki)

    n_known = int(binding_table[obp_name_list].notna().sum().sum())
    print(f"  → {n_known} cel·les amb Ki conegut (s'usaran per validar)")

    print(f"\n  Carregant matriu Tanimoto...")
    T_matrix, cas_to_pos = ti.load_tanimoto_data()
    idx_to_pos = ti.build_binding_idx_to_matrix_pos(binding_table, cas_col, cas_to_pos)
    print(f"  → {len(idx_to_pos)}/{len(binding_table)} VOCs amb posició vàlida a la matriu")


    resultats = []
    millor = {"rmse": np.inf}

    total_combos = len(T_MIN_VALUES) * len(K_VALUES) * len(P_VALUES)
    combo_i = 0

    for t_min in T_MIN_VALUES:
        for k in K_VALUES:
            for p in P_VALUES:
                combo_i += 1
                print(f"\n  [{combo_i}/{total_combos}] Provant t_min={t_min}  k={k}  p={p}")
                rmse, df_err = ti.run_loo(
                    binding_table, T_matrix, idx_to_pos, obp_name_list,
                    t_min=t_min, k=k, p=p,
                )
                if rmse is None:
                    print(f"     Sense prediccions vàlides amb aquesta combinació")
                    continue

                n_pred = len(df_err)
                cobertura = n_pred / n_known
                print(f"    → RMSE={rmse:.3f}  (n={n_pred}, cobertura={cobertura:.1%})")

                resultats.append({
                    "t_min": t_min, "k": k, "p": p,
                    "RMSE": round(rmse, 4),
                    "n_predits": n_pred,
                    "cobertura": round(cobertura, 4),
                })

                if rmse < millor["rmse"]:
                    millor = {"rmse": rmse, "t_min": t_min, "k": k, "p": p, "df_err": df_err}

    # Resum de totes les combinacions
    df_resum = pd.DataFrame(resultats).sort_values("RMSE")
    out1 = "results/loo_resum_combinacions.csv"
    df_resum.to_csv(out1, index=False)

    print("\n" + "=" * 70)
    print("  RESUM (ordenat per RMSE, més baix = millor)")
    print("=" * 70)
    print(df_resum.to_string(index=False))
    print(f"\n   Resum complet guardat a {out1}")

    if millor["rmse"] == np.inf:
        print("\n   Cap combinació ha donat resultats. Revisa la matriu Tanimoto.")
        return

    #  Detall de la combinació guanyadora 
    print("\n" + "=" * 70)
    print(f"  MILLOR COMBINACIÓ: t_min={millor['t_min']}  k={millor['k']}  "
          f"p={millor['p']}  (RMSE={millor['rmse']:.3f})")
    print("=" * 70)

    df_err = millor["df_err"]
    out2 = "results/loo_errors_detall.csv"
    df_err.to_csv(out2, index=False)
    print(f"  ✓ Errors individuals guardats a {out2}")

    resum_bins = ti.rmse_by_tmax_bin(df_err)
    out3 = "results/loo_rmse_per_tmax.csv"
    resum_bins.to_csv(out3)
    print(f"\n  RMSE desglossat per similitud del millor veí (T_max):")
    print(resum_bins.to_string())
    print(f"\n  ✓ Guardat a {out3}")

    print("\n" + "=" * 70)
    print("  INTERPRETACIÓ")
    print("=" * 70)
    print("  RMSE ≤ 0.3 (log10 Ki)  →  factor d'error ≤ ~2x   → BO")
    print("  RMSE 0.3–0.5           →  factor d'error ~2-3x   → acceptable")
    print("  RMSE > 0.5             →  factor d'error > 3x    → massa soroll")
    print(f"\n  → Paràmetres recomanats per al pas 3 (tanimoto_omplir.py):")
    print(f"      T_MIN = {millor['t_min']}")
    print(f"      K     = {millor['k']}")
    print(f"      P     = {millor['p']}")


if __name__ == "__main__":
    main()