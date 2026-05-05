import os
import sys
import re
import pandas as pd
import numpy as np


# CONSTANTS 

BINDING_FILE  = "Compound_OBP_binding.csv"
INFO_FILE     = "OBP_info_new.csv"

BEST_OBP_TYPE = "Classic OBP"
BEST_CYS_NUM  = 6

BIG_KI_VALUE  = 1000.0   # μM — valor quan el paper diu ">XX", significa "no s'uneix gaire"

# Si un interferent té Ki menys que (Ki_diana * aquest factor), descartem l'OBP
# 1.0 = estrictament millor que el diana; 2.0 = tolera fins a 2 vegades millor
SELECTIVITY_LIMIT = 1.0


#  LLEGIR I NETEJAR ELS CSV 

def convert_ki_to_float(raw_value):
    """
    Converteix un valor de Ki (afinitat) a float.
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

    """
    print(f"Llegint matriu de binding: {binding_file_path}")
    raw_binding = pd.read_csv(binding_file_path)
    print(f"  → {raw_binding.shape[0]} VOCs i {raw_binding.shape[1] - 2} OBPs")

    print(f"Llegint informació d'OBPs: {info_file_path}")
    obp_info_table = pd.read_csv(info_file_path)
    print(f"  → {len(obp_info_table)} OBPs amb metadades\n")

    cas_col      = raw_binding.columns[0]   # primera columna: número CAS
    name_col     = raw_binding.columns[1]   # segona columna: nom del compost
    obp_name_list = list(raw_binding.columns[2:])  # la resta: una columna per OBP

    # Convertim totes les cel·les de les columnes d'OBP a float
    binding_table = raw_binding.copy()
    for col_name in obp_name_list:
        binding_table[col_name] = binding_table[col_name].apply(convert_ki_to_float)

    return binding_table, obp_info_table, cas_col, name_col, obp_name_list


#  CERCAR UN VOC A LA MATRIU 

def find_voc_rows(binding_table, name_col, search_text):
    """
    Busca files a la matriu on el nom del VOC contingui search_text.
    La cerca no distingeix majúscules/minúscules.
    Retorna un DataFrame amb totes les coincidències.
    """
    is_match = binding_table[name_col].str.contains(
        search_text, case=False, na=False, regex=False
    )
    return binding_table[is_match]


#  LLEGIR EL FITXER D'INTERFERENTS 

def read_interferent_file(file_path):
    """
    Llegeix un fitxer de text amb un interferent per línia.
    Les línies que comencen per # o estan buides s'ignoren.
    Retorna una llista de strings.
    """
    interferent_list = []
    with open(file_path, encoding='utf-8') as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line or clean_line.startswith('#'):
                continue
            interferent_list.append(clean_line)
    return interferent_list


#  CALCULAR EL RANKING 

def build_obp_ranking(ki_values_diana, obp_info_table, binding_table,
                      name_col, interferent_list, obp_name_list):
    """
    Construeix el ranking d'OBPs per al VOC diana.

    Arguments:
      ki_values_diana  : pd.Series amb {nom_OBP → Ki en μM} del VOC diana
      obp_info_table   : DataFrame amb metadades de cada OBP
      binding_table    : matriu completa (per buscar interferents)
      name_col         : nom de la columna de noms de VOC
      interferent_list : llista de noms d'interferents
      obp_name_list    : llista de totes les columnes d'OBP

    Retorna un DataFrame ordenat amb les columnes:
      OBP, Ki_diana_uM, Type, Cystines, Preferred,
      Selectivity_OK, Min_Ki_interferent_uM, Worst_interferent, Score
    """

    # Recollim dades bàsiques de cada OBP 
    # Posem el nom de l'OBP com a índex per buscar ràpid
    info_by_name = obp_info_table.set_index('Binding Protein Name')

    obp_rows = []
    for obp_name in obp_name_list:

        ki_diana = ki_values_diana.get(obp_name, np.nan)

        # Si no hi ha Ki mesurada per a aquest VOC, saltem l'OBP
        if pd.isna(ki_diana):
            continue

        # Busquem les metadades a la taula d'info
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

        # L'OBP és "preferida" si és Classic OBP amb exactament 6 cisteïnes
        is_preferred = (obp_type == BEST_OBP_TYPE and cys_count == BEST_CYS_NUM)

        obp_rows.append({
            'OBP':         obp_name,
            'Ki_diana_uM': ki_diana,
            'Type':        obp_type,
            'Cystines':    cys_count,
            'Preferred':   is_preferred,
            'Species':     species,
            'UniProtID':   uniprot_id,
            'Alphafold':   alphafold,
        })

    if not obp_rows:
        return pd.DataFrame()  # cap OBP té dades per a aquest VOC

    result_table = pd.DataFrame(obp_rows)

    # Filtre de selectivitat 
    if interferent_list:

       
        ki_per_interferent = {}   # {nom_interferent → pd.Series de Kis per OBP}

        for interferent_name in interferent_list:
            found_rows = find_voc_rows(binding_table, name_col, interferent_name)

            if found_rows.empty:
                print(f"  Interferent '{interferent_name}' no trobat — ignorat.")
                continue

            if len(found_rows) > 1:
                print(f"  Interferent '{interferent_name}' té {len(found_rows)} coincidències, s'usa la primera.")

            first_row = found_rows.iloc[0]
            ki_series = pd.Series({col: first_row[col] for col in obp_name_list})
            ki_per_interferent[interferent_name] = ki_series
            print(f"  Interferent '{interferent_name}': trobat.")

        # Per cada OBP, comprovem si algun interferent s'hi uneix millor que el diana
        selectivity_ok_list     = []
        min_interferent_ki_list = []
        worst_interferent_list  = []

        for _, obp_row in result_table.iterrows():
            current_obp   = obp_row['OBP']
            ki_diana_val  = obp_row['Ki_diana_uM']

            passes_filter = True
            min_ki_found  = np.nan
            worst_name    = "-"

            for interferent_name, ki_series in ki_per_interferent.items():
                ki_interferent = ki_series.get(current_obp, np.nan)

                # Si no hi ha dada per a aquest parell, no penalitzem
                if pd.isna(ki_interferent):
                    continue

                # Guardem el pitjor interferent (Ki més baixa = s'uneix més fort)
                if pd.isna(min_ki_found) or ki_interferent < min_ki_found:
                    min_ki_found = ki_interferent
                    worst_name   = interferent_name

                # Si l'interferent té Ki menor que el diana → l'OBP no és selectiva
                if ki_interferent < ki_diana_val * SELECTIVITY_LIMIT:
                    passes_filter = False

            selectivity_ok_list.append(passes_filter)
            min_interferent_ki_list.append(min_ki_found)
            worst_interferent_list.append(worst_name)

        result_table['Selectivity_OK']        = selectivity_ok_list
        result_table['Min_Ki_interferent_uM'] = min_interferent_ki_list
        result_table['Worst_interferent']     = worst_interferent_list

    else:
        
        result_table['Selectivity_OK']        = True
        result_table['Min_Ki_interferent_uM'] = np.nan
        result_table['Worst_interferent']     = "-"

    result_table['Score'] = result_table['Ki_diana_uM'].copy()
    preferred_mask = result_table['Preferred']
    result_table.loc[preferred_mask, 'Score'] *= 0.70

    #  Ordenació
    # Primer les que passen selectivitat (True > False)
    # Dins de cada grup, les de Score més baix primer
    result_table = result_table.sort_values(
        by=['Selectivity_OK', 'Score'],
        ascending=[False, True]
    ).reset_index(drop=True)

    return result_table


# MOSTRAR ELS RESULTATS 

def show_results(result_table, voc_name, how_many):
    """Imprimeix el ranking per pantalla."""

    separator = "═" * 78
    print(f"\n{separator}")
    print(f"  RANKING OBP per a: {voc_name}")
    print(separator)

    best_obp = result_table.iloc[0]
    print(f"\n  MILLOR CANDIDAT: {best_obp['OBP']}")
    print(f"    Ki diana     : {best_obp['Ki_diana_uM']:.2f} μM")
    print(f"    Tipus OBP    : {best_obp['Type']}")
    cys_display = int(best_obp['Cystines']) if not pd.isna(best_obp['Cystines']) else "?"
    print(f"    Cisteïnes    : {cys_display}")
    print(f"    Espècie      : {best_obp['Species']}")
    print(f"    Preferit     : {'Sí (Classic OBP + 6 Cys)' if best_obp['Preferred'] else 'No'}")
    print(f"    Selectivitat : {'OK' if best_obp['Selectivity_OK'] else 'Falsos positius possibles'}")
    if not pd.isna(best_obp['Min_Ki_interferent_uM']):
        print(f"    Pitjor interf: {best_obp['Worst_interferent']} "
              f"(Ki={best_obp['Min_Ki_interferent_uM']:.2f} μM)")
    if best_obp['UniProtID'] not in (None, '-', 'nan', ''):
        print(f"    UniProt ID   : {best_obp['UniProtID']}")
    if best_obp['Alphafold'] not in (None, '-', 'nan', ''):
        print(f"    AlphaFold    : {best_obp['Alphafold']}")

    print(f"\n{'─'*78}")
    print(f"  TOP {how_many} OBPs (score baix = millor candidat):")
    print(f"{'─'*78}")
    print(f"  {'#':>3}  {'OBP':<20} {'Ki(μM)':>8}  {'Tipus':<15} "
          f"{'Cys':>4}  {'Pref':>5}  {'Sel':>4}  {'Ki_interf':>10}")
    print(f"{'─'*78}")

    for position, (_, row) in enumerate(result_table.head(how_many).iterrows()):
        pref_mark = "*" if row['Preferred'] else " "
        sel_mark  = "OK" if row['Selectivity_OK'] else "X"
        cys_val   = int(row['Cystines']) if not pd.isna(row['Cystines']) else "?"
        ki_interf = (f"{row['Min_Ki_interferent_uM']:.2f}"
                     if not pd.isna(row['Min_Ki_interferent_uM']) else "  —   ")
        print(f"  {position+1:>3}. {row['OBP']:<20} {row['Ki_diana_uM']:>8.2f}  "
              f"{str(row['Type']):<15} {str(cys_val):>4}  {pref_mark:>5}  "
              f"{sel_mark:>4}  {ki_interf:>10}")

    print(f"{'─'*78}")
    total         = len(result_table)
    n_preferred   = result_table['Preferred'].sum()
    n_selective   = result_table['Selectivity_OK'].sum()
    print(f"\n  Total OBPs amb dades : {total}")
    print(f"  Classic OBP + 6 Cys  : {n_preferred}")
    print(f"  Passen selectivitat  : {n_selective}")
    print(f"\n  Nota: * = Classic OBP amb 6 Cys. Score = Ki × 0.70 si preferida.")
    print(separator)


#  GUARDAR EL CSV 

def save_ranking_to_csv(result_table, voc_name_safe):
    """Guarda el ranking complet en un CSV."""
    output_filename = f"ranking_{voc_name_safe}.csv"
    result_table.to_csv(output_filename, index=False)
    print(f"\nRanking guardat a: {output_filename}")


# PROGRAMA PRINCIPAL 

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
    binding_table, obp_info_table, cas_col, name_col, obp_name_list = load_csv_files(
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

        # Diverses coincidències: l'usuari tria
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

    voc_display  = chosen_voc[name_col]
    voc_name_safe = re.sub(r'[^\w]+', '_', voc_display)[:40]  # per al nom del fitxer CSV

    # Extraiem la fila de Kis del VOC diana com a pd.Series
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
    )

    if result_table.empty:
        print("  Cap OBP té dades per a aquest VOC.")
        sys.exit(0)

    # Mostrem els resultats
    show_results(result_table, voc_display, how_many)

    # Guardem el CSV si l'usuari vol
    save_answer = input("\nVols guardar el ranking en CSV? [S/N]: ").strip().lower()
    if save_answer not in ('n', 'no'):
        save_ranking_to_csv(result_table, voc_name_safe)

    print("\nFet!\n")


if __name__ == "__main__":
    main()