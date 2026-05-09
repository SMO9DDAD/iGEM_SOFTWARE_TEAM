import os
import sys
import re
import pandas as pd
import numpy as np




BINDING_FILE  = "Compound_OBP_binding.csv"
INFO_FILE     = "OBP_info_new.csv"

BEST_OBP_TYPE = "Classic OBP"
BEST_CYS_NUM  = 6

BIG_KI_VALUE  = 1000.0   # μM — valor quan el paper diu ">XX"

# Tau (τ) per al criteri de selectivitat: s2 = min(1, Ki_interf / (Ki_diana × τ))
# Com més gran τ, més exigent és el criteri de selectivitat.
SELECTIVITY_TAU = 10.0



DEFAULT_WEIGHTS = {
    "w_affinity":     0.50,   # s1: Ki diana (afinitat pel VOC diana)
    "w_selectivity":  0.33,   # s2: selectivitat vs interferents
    "w_promiscuity":  0.17,   # s5: promiscuïtat (penalització per unió a molts VOCs)
    # afegir nous pesos
}


# LLEGIR I NETEJAR ELS CSV 

def convert_ki_to_float(raw_value):
    """
    Converteix un valor de Ki (afinitat) a float en μM.

    """
    if pd.isna(raw_value):
        return np.nan

    text = str(raw_value).strip().replace('\xa0', '').replace(' ', '')

    if text.startswith('>'):
        try:
            number = float(re.sub(r'[^\d.]', '', text))
            return number * 1.1
        except ValueError:
            return BIG_KI_VALUE

    try:
        return float(text)
    except ValueError:
        return np.nan


def load_csv_files(binding_file_path, info_file_path):
    """
    Llegeix els dos CSV i retorna:
      - binding_table   : matriu VOCs × OBPs amb Kis en float
      - obp_info_table  : taula de metadades de cada OBP
      - cas_col         : nom de la columna CAS
      - name_col        : nom de la columna de noms de VOC
      - obp_name_list   : llista dels noms de les columnes d'OBP
    """
    print(f"Llegint matriu de binding: {binding_file_path}")
    raw_binding = pd.read_csv(binding_file_path)
    print(f"  → {raw_binding.shape[0]} VOCs i {raw_binding.shape[1] - 2} OBPs")

    print(f"Llegint informació d'OBPs: {info_file_path}")
    obp_info_table = pd.read_csv(info_file_path)
    print(f"  → {len(obp_info_table)} OBPs amb metadades\n")

    cas_col       = raw_binding.columns[0]
    name_col      = raw_binding.columns[1]
    obp_name_list = list(raw_binding.columns[2:])

    binding_table = raw_binding.copy()
    for col_name in obp_name_list:
        binding_table[col_name] = binding_table[col_name].apply(convert_ki_to_float)

    # Rang real de Ki a tota la matriu — s'usa per normalitzar s1
    all_ki_values = binding_table[obp_name_list].values.flatten()
    all_ki_values = all_ki_values[~__import__('numpy').isnan(all_ki_values)]
    all_ki_values = all_ki_values[all_ki_values > 0]
    ki_min_matrix = float(all_ki_values.min())
    ki_max_matrix = float(all_ki_values.max())
    print(f"  → Rang Ki a la matriu: {ki_min_matrix:.3f} – {ki_max_matrix:.1f} μM")

    return binding_table, obp_info_table, cas_col, name_col, obp_name_list, ki_min_matrix, ki_max_matrix


# CERCAR UN VOC A LA MATRIU 

def find_voc_rows(binding_table, name_col, search_text):
    """
    Retorna les files on el nom del VOC contingui search_text.
    """
    is_match = binding_table[name_col].str.contains(
        search_text, case=False, na=False, regex=False
    )
    return binding_table[is_match]


#  LLEGIR EL FITXER D'INTERFERENTS 

def read_interferent_file(file_path):

    interferent_list = []
    with open(file_path, encoding='utf-8') as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line or clean_line.startswith('#'):
                continue
            interferent_list.append(clean_line)
    return interferent_list


# DEMANAR PESOS A L'USUARI ─

def ask_user_for_weights():


    print(f"  Pesos per defecte:")
    print(f"    Afinitat (s1)     : {DEFAULT_WEIGHTS['w_affinity']:.2f}")
    print(f"    Selectivitat (s2) : {DEFAULT_WEIGHTS['w_selectivity']:.2f}")
    print(f"    Promiscuïtat (s5) : {DEFAULT_WEIGHTS['w_promiscuity']:.2f}")
    print(f"  Suma: {sum(DEFAULT_WEIGHTS.values()):.2f}")


    use_default = input("Usar pesos per defecte? [S/n]: ").strip().lower()
    if use_default not in ('n', 'no'):
        return dict(DEFAULT_WEIGHTS)  # retorna còpia per no modificar l'original

    # L'usuari introdueix pesos personalitzats
    custom_weights = {}
    weight_names = {
        "w_affinity":    "Afinitat (s1)",
        "w_selectivity": "Selectivitat (s2)",
        "w_promiscuity": "Promiscuïtat (s5)",
        # Afegeix aquí la descripció de nous criteris
    }

    print("\n  Introdueix els pesos (números entre 0 i 1).")
    print("  No cal que sumin 1: el programa els normalitzarà automàticament.\n")

    for key, label in weight_names.items():
        while True:
            raw = input(f"  Pes per a {label} [per defecte {DEFAULT_WEIGHTS[key]:.2f}]: ").strip()
            if raw == "":
                custom_weights[key] = DEFAULT_WEIGHTS[key]
                break
            try:
                value = float(raw)
                if value < 0:
                    print("    El pes no pot ser negatiu.")
                    continue
                custom_weights[key] = value
                break
            except ValueError:
                print("    Introdueix un número vàlid (ex: 0.4).")

    # Normalitzem perquè sumin 1.0
    total = sum(custom_weights.values())
    if total == 0:
        print("  Tots els pesos són 0. S'usaran els pesos per defecte.")
        return dict(DEFAULT_WEIGHTS)

    normalized = {k: v / total for k, v in custom_weights.items()}

    print(f"\n  Pesos normalitzats (suma = 1.0):")
    for key, label in weight_names.items():
        print(f"    {label:<25}: {normalized[key]:.3f}")

    return normalized


#  CÀLCUL DELS SUB-SCORES 

def compute_s1_affinity(ki_diana, ki_min_matrix, ki_max_matrix):
    if pd.isna(ki_diana) or ki_diana <= 0:
        return 0.0
    if ki_min_matrix <= 0 or pd.isna(ki_min_matrix):
        return 0.0
    if ki_max_matrix <= 0 or pd.isna(ki_max_matrix):
        return 0.0

    log_diana = -np.log10(ki_diana)
    log_best  = -np.log10(ki_min_matrix)   # valor més alt (millor afinitat)
    log_worst = -np.log10(ki_max_matrix)   # valor més baix (pitjor afinitat)

    denominator = log_best - log_worst
    if denominator == 0:
        return 0.0

    s1 = (log_diana - log_worst) / denominator
    return float(np.clip(s1, 0.0, 1.0))


def compute_s2_selectivity(ki_diana, ki_min_interferent, tau=SELECTIVITY_TAU):

    if pd.isna(ki_min_interferent) or pd.isna(ki_diana) or ki_diana <= 0:
        return 1.0   # sense interferents → no penalitzem

    ratio = ki_min_interferent / (ki_diana * tau)
    return float(min(1.0, ratio))


def compute_s5_promiscuity(n_vocs_bound, n_vocs_total):

    if n_vocs_total == 0:
        return 1.0

    promiscuity_ratio = n_vocs_bound / n_vocs_total
    return float(1.0 - promiscuity_ratio)


def compute_final_score(s1, s2, s5, weights):

    score = (
        weights['w_affinity']    * s1 +
        weights['w_selectivity'] * s2 +
        weights['w_promiscuity'] * s5
        # + weights['w_nou_criteri'] * sX
    )
    return float(score)


#  CALCULAR EL RANKING

def build_obp_ranking(ki_values_diana, obp_info_table, binding_table,
                      name_col, interferent_list, obp_name_list, weights,
                      ki_min_matrix, ki_max_matrix):


    # Recollim dades bàsiques de cada OBP 
    info_by_name  = obp_info_table.set_index('Binding Protein Name')
    n_vocs_total  = len(binding_table)   # total de files = total de VOCs 

    obp_rows = []
    for obp_name in obp_name_list:

        ki_diana = ki_values_diana.get(obp_name, np.nan)

        if pd.isna(ki_diana):
            continue   # sense dada per a aquest OBP–VOC, el saltem

        # Metadades de l'OBP
        if obp_name in info_by_name.index:
            obp_row    = info_by_name.loc[obp_name]
            obp_type   = obp_row['Binding Protein Type']
            cys_count  = int(obp_row['Cystine count'])
            species    = obp_row['Species']
            uniprot_id = obp_row['UniProtID']
            alphafold  = obp_row['Alphafold']
        else:
            obp_type   = "?"
            cys_count  = np.nan
            species    = "?"
            uniprot_id = "-"
            alphafold  = "-"

        is_preferred = (obp_type == BEST_OBP_TYPE and cys_count == BEST_CYS_NUM)

        # Promiscuïtat: quants VOCs té mesurats aquest OBP a tota la matriu
        n_vocs_bound = int(binding_table[obp_name].notna().sum())

        obp_rows.append({
            'OBP':          obp_name,
            'Ki_diana_uM':  ki_diana,
            'Type':         obp_type,
            'Cystines':     cys_count,
            'Preferred':    is_preferred,
            'Species':      species,
            'UniProtID':    uniprot_id,
            'Alphafold':    alphafold,
            'N_VOCs_bound': n_vocs_bound,
        })

    if not obp_rows:
        return pd.DataFrame()

    result_table = pd.DataFrame(obp_rows)

    # Selectivitat: Ki mínima dels interferents per cada OBP 
    ki_per_interferent = {}

    if interferent_list:
        print(f"\n[·] Aplicant filtre de selectivitat ({len(interferent_list)} interferents)...")
        for interferent_name in interferent_list:
            found_rows = find_voc_rows(binding_table, name_col, interferent_name)

            if found_rows.empty:
                print(f"  Interferent '{interferent_name}' no trobat — ignorat.")
                continue
            if len(found_rows) > 1:
                print(f"  Interferent '{interferent_name}': {len(found_rows)} coincidències, s'usa la primera.")

            first_row = found_rows.iloc[0]
            ki_series = pd.Series({col: first_row[col] for col in obp_name_list})
            ki_per_interferent[interferent_name] = ki_series
            print(f"  '{interferent_name}': trobat.")

    # Per cada OBP calculem la Ki mínima dels interferents i el pitjor interferent
    min_interferent_ki_list = []
    worst_interferent_list  = []

    for _, obp_row in result_table.iterrows():
        current_obp  = obp_row['OBP']
        min_ki_found = np.nan
        worst_name   = "-"

        for interf_name, ki_series in ki_per_interferent.items():
            ki_i = ki_series.get(current_obp, np.nan)
            if pd.isna(ki_i):
                continue
            if pd.isna(min_ki_found) or ki_i < min_ki_found:
                min_ki_found = ki_i
                worst_name   = interf_name

        min_interferent_ki_list.append(min_ki_found)
        worst_interferent_list.append(worst_name)

    result_table['Min_Ki_interferent_uM'] = min_interferent_ki_list
    result_table['Worst_interferent']     = worst_interferent_list

    #  Calculem els tres sub-scores per a cada OBP 
    s1_list = []
    s2_list = []
    s5_list = []

    for _, obp_row in result_table.iterrows():
        s1 = compute_s1_affinity(
            ki_diana=obp_row['Ki_diana_uM'],
            ki_min_matrix=ki_min_matrix,
            ki_max_matrix=ki_max_matrix
        )
        s2 = compute_s2_selectivity(
            ki_diana=obp_row['Ki_diana_uM'],
            ki_min_interferent=obp_row['Min_Ki_interferent_uM'],
            tau=SELECTIVITY_TAU
        )
        s5 = compute_s5_promiscuity(
            n_vocs_bound=obp_row['N_VOCs_bound'],
            n_vocs_total=n_vocs_total
        )
        s1_list.append(s1)
        s2_list.append(s2)
        s5_list.append(s5)

    result_table['s1_affinity']    = s1_list
    result_table['s2_selectivity'] = s2_list
    result_table['s5_promiscuity'] = s5_list

    #  Score final 
    final_scores = []
    for _, obp_row in result_table.iterrows():
        score = compute_final_score(
            s1=obp_row['s1_affinity'],
            s2=obp_row['s2_selectivity'],
            s5=obp_row['s5_promiscuity'],
            weights=weights
        )
        final_scores.append(score)

    result_table['Score'] = final_scores

    #  Ordenació: Score descendent (més alt = millor) 
    result_table = result_table.sort_values(
        by='Score', ascending=False
    ).reset_index(drop=True)

    return result_table


#  MOSTRAR ELS RESULTATS 

def show_results(result_table, voc_name, how_many, weights):
    """Imprimeix el ranking per pantalla."""

    separator = "═" * 82
    print(f"\n{separator}")
    print(f"  RANKING OBP per a: {voc_name}")
    print(f"  Pesos: afinitat={weights['w_affinity']:.2f}  "
          f"selectivitat={weights['w_selectivity']:.2f}  "
          f"promiscuïtat={weights['w_promiscuity']:.2f}")
    print(separator)

    best_obp = result_table.iloc[0]
    cys_display = int(best_obp['Cystines']) if not pd.isna(best_obp['Cystines']) else "?"

    print(f"\n  MILLOR CANDIDAT: {best_obp['OBP']}")
    print(f"    Score final  : {best_obp['Score']:.4f}  "
          f"(s1={best_obp['s1_affinity']:.3f}  "
          f"s2={best_obp['s2_selectivity']:.3f}  "
          f"s5={best_obp['s5_promiscuity']:.3f})")
    print(f"    Ki diana     : {best_obp['Ki_diana_uM']:.2f} μM")
    print(f"    Tipus OBP    : {best_obp['Type']}")
    print(f"    Cisteïnes    : {cys_display}")
    print(f"    Espècie      : {best_obp['Species']}")
    print(f"    Preferit     : {'Sí (Classic OBP + 6 Cys)' if best_obp['Preferred'] else 'No'}")
    print(f"    VOCs units   : {best_obp['N_VOCs_bound']} de {len(result_table)} mesurats")
    if not pd.isna(best_obp['Min_Ki_interferent_uM']):
        print(f"    Pitjor interf: {best_obp['Worst_interferent']} "
              f"(Ki={best_obp['Min_Ki_interferent_uM']:.2f} μM)")
    if best_obp['UniProtID'] not in (None, '-', 'nan', ''):
        print(f"    UniProt ID   : {best_obp['UniProtID']}")
    if best_obp['Alphafold'] not in (None, '-', 'nan', ''):
        print(f"    AlphaFold    : {best_obp['Alphafold']}")

    print(f"\n{'─'*82}")
    print(f"  TOP {how_many} OBPs — Score alt = millor candidat")
    print(f"{'─'*82}")
    print(f"  {'#':>3}  {'OBP':<20} {'Score':>7}  {'s1':>5}  {'s2':>5}  {'s5':>5}  "
          f"{'Ki(μM)':>8}  {'Cys':>4}  {'Pref':>4}  {'Ki_interf':>9}")
    print(f"{'─'*82}")

    for position, (_, row) in enumerate(result_table.head(how_many).iterrows()):
        pref_mark = "*" if row['Preferred'] else " "
        cys_val   = int(row['Cystines']) if not pd.isna(row['Cystines']) else "?"
        ki_interf = (f"{row['Min_Ki_interferent_uM']:.2f}"
                     if not pd.isna(row['Min_Ki_interferent_uM']) else "   —  ")
        print(f"  {position+1:>3}. {row['OBP']:<20} {row['Score']:>7.4f}  "
              f"{row['s1_affinity']:>5.3f}  {row['s2_selectivity']:>5.3f}  "
              f"{row['s5_promiscuity']:>5.3f}  {row['Ki_diana_uM']:>8.2f}  "
              f"{str(cys_val):>4}  {pref_mark:>4}  {ki_interf:>9}")

    print(f"{'─'*82}")

    total       = len(result_table)
    n_preferred = result_table['Preferred'].sum()
    print(f"\n  Total OBPs amb dades : {total}")
    print(f"  Classic OBP + 6 Cys  : {n_preferred}")
    print(f"\n  Columnes sub-score:")
    print(f"    s1 = afinitat [-log10(Ki) normalitzat a [0,1]]")
    print(f"    s2 = selectivitat [min(1, Ki_interf / (Ki_diana × τ={SELECTIVITY_TAU:.0f}))]")
    print(f"    s5 = 1 − promiscuïtat [1 − VOCs_units/VOCs_totals]")
    print(f"    *  = Classic OBP amb 6 cisteïnes")
    print(separator)



# MAIN 

def main():

    print("OBP FINDER — iGEM URV 2025")
    print("Selecció del millor candidat OBP per a un VOC diana\n")

    # Comprovem que els fitxers CSV existeixen
    for csv_path in (BINDING_FILE, INFO_FILE):
        if not os.path.isfile(csv_path):
            print(f"ERROR: No es troba el fitxer '{csv_path}'.")
            print("Posa els CSV a la mateixa carpeta que aquest script.")
            sys.exit(1)

    # Llegim els CSV
    binding_table, obp_info_table, cas_col, name_col, obp_name_list, ki_min_matrix, ki_max_matrix = load_csv_files(
        BINDING_FILE, INFO_FILE
    )

    # Demanem el VOC diana a l'usuari
    while True:
        user_query = input("Nom del VOC diana (o part del nom): ").strip()
        if not user_query:
            print("  Cal escriure alguna cosa.")
            continue

        matches = find_voc_rows(binding_table, name_col, user_query)

        if matches.empty:
            print(f"  No s'ha trobat cap VOC amb '{user_query}'. Torna-ho a intentar.")
            continue

        if len(matches) == 1:
            chosen_voc = matches.iloc[0]
            print(f"\n  VOC seleccionat: {chosen_voc[name_col]}")
            print(f"  CAS: {chosen_voc[cas_col]}")
            break

        print(f"\n  {len(matches)} coincidències trobades:")
        for i, (_, row) in enumerate(matches.iterrows()):
            print(f"    [{i+1}]  {row[name_col][:70]}  (CAS: {row[cas_col]})")
        while True:
            user_choice = input(f"  Tria un número [1-{len(matches)}]: ").strip()
            if user_choice.isdigit() and 1 <= int(user_choice) <= len(matches):
                chosen_voc = matches.iloc[int(user_choice) - 1]
                print(f"\n  VOC seleccionat: {chosen_voc[name_col]}")
                break
            print("  Número invàlid.")
        break

    voc_display   = chosen_voc[name_col]
    voc_name_safe = re.sub(r'[^\w]+', '_', voc_display)[:40]

    ki_values_diana = pd.Series({col: chosen_voc[col] for col in obp_name_list})

    # Demanem el fitxer d'interferents
    interf_path = input(
        "\nFitxer d'interferents (deixa buit per saltar): "
    ).strip().strip('"').strip("'")

    interferent_list = []
    if interf_path:
        if not os.path.isfile(interf_path):
            print(f"  Fitxer '{interf_path}' no trobat. Continuem sense interferents.")
        else:
            interferent_list = read_interferent_file(interf_path)
            print(f"  {len(interferent_list)} interferents carregats: {interferent_list}")

    # Demanem els pesos del score
    weights = ask_user_for_weights()

    # Quants resultats vol veure l'usuari
    top_input = input("\nQuants candidats vols veure? [per defecte: 10]: ").strip()
    how_many  = int(top_input) if top_input.isdigit() and int(top_input) > 0 else 10

    # Calculem el ranking
    print("\nCalculant ranking...")
    result_table = build_obp_ranking(
        ki_values_diana=ki_values_diana,
        obp_info_table=obp_info_table,
        binding_table=binding_table,
        name_col=name_col,
        interferent_list=interferent_list,
        obp_name_list=obp_name_list,
        weights=weights,
        ki_min_matrix=ki_min_matrix,
        ki_max_matrix=ki_max_matrix,
    )

    if result_table.empty:
        print("  Cap OBP té dades per a aquest VOC.")
        sys.exit(0)

    # Mostrem els resultats
    show_results(result_table, voc_display, how_many, weights)


if __name__ == "__main__":
    main()