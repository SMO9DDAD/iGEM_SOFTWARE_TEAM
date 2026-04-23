import argparse
import sys
import io
import csv
import math
from pathlib import Path

# Windows UTF-8 fix
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")



CAS_A_NOM = {
    # Alcohols
    "3391-86-4": "1-octen-3-ol",
    "111-27-3":  "1-hexanol",
    "111-70-6":  "1-heptanol",
    "111-87-5":  "1-octanol",
    "143-08-8":  "1-nonanol",
    "112-30-1":  "1-decanol",
    "112-42-5":  "1-undecanol",
    "112-53-8":  "1-dodecanol",
    "112-72-1":  "1-tetradecanol",
    "36653-82-4":"1-hexadecanol",
    "543-49-7":  "2-heptanol",
    "584-02-1":  "3-pentanol",
    "623-37-0":  "3-hexanol",
    "104-76-7":  "2-ethyl-1-hexanol",
    "616-25-1":  "1-penten-3-ol",
    "928-94-9":  "trans-2-hexen-1-ol",
    "928-95-0":  "trans-2-hexen-1-ol",
    "928-96-1":  "cis-3-hexen-1-ol",
    "6789-88-4": "1-octen-3-ol",
    "6789-80-6": "cis-3-hexen-1-ol",
    "89-78-1":   "menthol",
    "60-12-8":   "2-phenylethanol",
    "100-51-6":  "benzyl alcohol",
    "106-22-9":  "citronellol",
    "106-24-1":  "geraniol",
    "106-25-2":  "nerol",
    "78-70-6":   "linalool",
    "7785-70-8": "alpha-terpineol",
    "98-55-5":   "alpha-terpineol",
    "562-74-3":  "4-terpineol",
    "77-53-2":   "cedrol",
    "7212-44-4": "nerolidol",
    "4602-84-0": "farnesol",
    "81-13-0":   "panthenol",
    "68-26-8":   "retinol",
    # Aldehids
    "66-25-1":   "hexanal",
    "111-71-7":  "heptanal",
    "124-13-0":  "octanal",
    "124-19-6":  "nonanal",
    "112-31-2":  "decanal",
    "112-54-9":  "dodecanal",
    "110-62-3":  "valeraldehyde",
    "96-17-3":   "2-methylbutanal",
    "122-78-1":  "phenylacetaldehyde",
    "5392-40-5": "citral",
    "116-26-7":  "safranal",
    "6728-26-3": "trans-2-hexenal",
    "18829-56-6":"trans-2-nonenal",
    "4340-84-5": "2,6-nonadienal",
    "104-55-2":  "cinnamaldehyde",
    "100-52-7":  "benzaldehyde",
    "90-02-8":   "salicylaldehyde",
    "106-23-0":  "citronellal",
    "20126-76-5":"neral",
    # Cetones
    "589-38-8":  "3-hexanone",
    "591-78-6":  "2-hexanone",
    "110-43-0":  "2-heptanone",
    "111-13-7":  "2-octanone",
    "502-56-7":  "5-nonanone",
    "2345-28-0": "2-decanone",
    "3016-19-1": "2-undecanone",
    "2173-56-0": "2-dodecanone",
    "2445-76-3": "2-tridecanone",
    "110-93-0":  "6-methyl-5-hepten-2-one",
    "98-86-2":   "acetophenone",
    "127-41-3":  "alpha-ionone",
    "79-77-6":   "beta-ionone",
    "76-22-2":   "camphor",
    "2244-16-8": "carvone",
    "490-10-8":  "pulegone",
    "34010-20-3":"isophorone",
    "34010-21-4":"4-oxoisophorone",
    "40642-40-8":"nepetalactone",
    # Terpens
    "138-86-3":  "limonene",
    "5989-27-5": "limonene",
    "5989-54-8": "limonene",
    "80-56-8":   "alpha-pinene",
    "7785-26-4": "alpha-pinene",
    "127-91-3":  "beta-pinene",
    "13877-91-3":"beta-pinene",
    "13466-78-9":"delta-3-carene",
    "87-44-5":   "beta-caryophyllene",
    "502-61-4":  "alpha-farnesene",
    "30507-70-1":"alpha-farnesene",
    "123-35-3":  "myrcene",
    "6186-98-7": "beta-ocimene",
    "99-83-2":   "alpha-phellandrene",
    "99-85-4":   "gamma-terpinene",
    "99-86-5":   "alpha-terpinene",
    "99-87-6":   "para-cymene",
    "586-62-9":  "terpinolene",
    "470-82-6":  "eucalyptol",
    "79-92-5":   "camphene",
    "469-61-4":  "alpha-cedrene",
    "6753-98-6": "humulene",
    "18794-84-8":"delta-cadinene",
    "765-17-3":  "bicyclogermacrene",
    "3387-41-5": "sabinene",
    "5989-08-2": "alpha-terpineol",
    "10482-56-1":"alpha-terpineol",
    "35153-15-2":"beta-bisabolol",
    "35153-20-9":"alpha-bisabolol",
    "68279-24-3":"beta-farnesene",
    "28079-04-1":"alpha-humulene",
    "22679-54-5":"alpha-copaene",
    "3790-78-1": "alpha-copaene",
    # Esters
    "141-78-6":  "ethyl acetate",
    "123-86-4":  "butyl acetate",
    "142-92-7":  "hexyl acetate",
    "105-54-4":  "ethyl butyrate",
    "123-92-2":  "isoamyl acetate",
    "109-21-7":  "butyl butyrate",
    "4748-78-1": "ethyl hexanoate",
    "106-30-9":  "ethyl heptanoate",
    "5973-71-7": "octyl acetate",
    "105-87-3":  "geranyl acetate",
    "53939-28-9":"linalyl acetate",
    "53939-27-8":"cis-3-hexenyl acetate",
    "3681-71-8": "cis-3-hexenyl acetate",
    "93-89-0":   "ethyl benzoate",
    "93-58-3":   "methyl benzoate",
    "119-36-8":  "methyl salicylate",
    "94-50-8":   "benzyl acetate",
    "101-41-7":  "methyl phenylacetate",
    "120-51-4":  "benzyl benzoate",
    "136-60-7":  "butyl benzoate",
    "19902-08-0":"citronellyl acetate",
    "16725-53-4":"neryl acetate",
    "23192-82-7":"citronellyl formate",
    "20711-10-8":"geranyl formate",
    "14959-86-5":"geranyl propionate",
    "33189-72-9":"geranyl butyrate",
    "61301-56-2":"neryl butyrate",
    "38363-29-0":"citronellyl butyrate",
    "31501-11-8":"hexyl hexanoate",
    "2639-63-6": "hexyl butyrate",
    "58594-45-9":"alpha-terpinyl acetate",
    "72698-30-7":"methyl jasmonate",
    "14816-18-3":"methyl jasmonate",
    "876-02-8":  "methyl dihydrojasmonate",
    "1009-61-6": "methyl decanoate",
    "143-13-5":  "methyl nonanoate",
    "27465-51-6":"methyl laurate",
    "7061-54-3": "ethyl myristate",
    "2027-47-6": "dodecyl acetate",
    "502-99-8":  "dodecanyl acetate",
    "112-06-1":  "heptyl acetate",
    "94444-18-5":"hexyl salicylate",
    "4695-62-9": "methyl cinnamate",
    "101-97-3":  "ethyl phenylacetate",
    # Aromatiques i fenols
    "97-53-0":   "eugenol",
    "93-15-2":   "methyleugenol",
    "63408-45-7":"methyl eugenol",
    "937-30-4":  "methylchavicol",
    "121-33-5":  "vanillin",
    "120-72-9":  "indole",
    "91-20-3":   "naphthalene",
    "83-34-1":   "skatole",
    "100-47-0":  "benzonitrile",
    "7786-61-0": "2,6-dimethylphenol",
    "18172-67-3":"2-phenoxyethanol",
    "64-04-0":   "phenylethylamine",
    "5146-66-7": "coumarin",
    "56219-04-6":"myristicin",
    # Alcans
    "109-66-0":  "pentane",
    "110-54-3":  "hexane",
    "111-65-9":  "octane",
    "111-84-2":  "nonane",
    "124-18-5":  "decane",
    "1120-21-4": "undecane",
    "112-40-3":  "dodecane",
    "629-50-5":  "tridecane",
    "629-59-4":  "tetradecane",
    "629-62-9":  "pentadecane",
    "629-78-7":  "heptadecane",
    "629-92-5":  "nonadecane",
    "544-76-3":  "hexadecane",
    "112-95-8":  "eicosane",
    "629-80-1":  "hexadecanol",
    "10378-01-5":"tridecane",
    # Acids grassos
    "57-10-3":   "palmitic acid",
    "60-33-3":   "linoleic acid",
    "544-63-8":  "myristic acid",
    "65-85-0":   "benzoic acid",
    "143-07-7":  "dodecanoic acid",
    # Compostos de sofre
    "624-92-0":  "dimethyl disulfide",
    "3658-80-8": "dimethyl trisulfide",
    # Miscel·lania
    "116-31-4":  "phytol",
    "2497-18-9": "squalene",
    "84-74-2":   "dibutyl phthalate",
    "131-11-3":  "dimethyl phthalate",
    "128-37-0":  "BHT",
    "2921-88-2": "chlorpyrifos",
    "53398-83-7":"cis-jasmone",
    "6909-30-4": "hexadecatrienoic acid",
    "56683-54-6":"methyl linoleate",
    "503-74-2":  "isovaleric acid",
    "1576-95-0": "trans-2-penten-1-ol",
    "35237-64-0":"methyl heptadecanoate",
    "40716-66-3":"decanal dimethyl acetal",
    "16974-10-0":"alpha-ionol",
    "16491-36-4":"alpha-terpinyl formate",
    "67446-07-5":"alpha-isomethylionone",
    "104-67-6":  "gamma-undecalactone",
    "6378-65-0": "2-nonanol",
    "3025-30-7": "3-nonanol",
    "3913-71-1": "1-octen-3-ol",
    "3650-28-0": "nonanal",
    "1139-30-6": "alpha-cedrene",
    "5445-77-2": "dihydrojasmone",
    "10486-19-8":"citronellyl propionate",
    "1211-29-6": "cyclohexyl butyrate",
    "1135-66-6": "myrcenol",
    "12400-00-7":"linalool oxide",
    "38963-29-0":"citronellyl butyrate",
}

NOM_A_CAS = {}
for cas, nom in CAS_A_NOM.items():
    NOM_A_CAS[nom.lower()] = cas
    NOM_A_CAS[nom.lower().replace("-", " ")] = cas


def resoldre_a_cas(consulta, llista_vocs):

    text = consulta.strip()
    text_minus = text.lower()

    # 1. Coincidencia exacta al dataset (funciona si l'usuari escriu el CAS)
    for voc in llista_vocs:
        if voc.lower() == text_minus:
            return voc

    # 2. Nom exacte al diccionari
    if text_minus in NOM_A_CAS:
        cas = NOM_A_CAS[text_minus]
        if cas in llista_vocs:
            return cas

    # 3. Coincidencia parcial del nom als valors del diccionari
    candidats = []
    for clau_nom, cas in NOM_A_CAS.items():
        if text_minus in clau_nom and cas in llista_vocs:
            candidats.append((clau_nom, cas))
    if candidats:
        # preferim el nom mes curt (coincidencia mes especifica)
        candidats.sort(key=lambda x: len(x[0]))
        return candidats[0][1]

    # 4. Subcadena directament a llista_vocs (ultim recurs)
    coincidencies = [voc for voc in llista_vocs if text_minus in voc.lower()]
    if coincidencies:
        return sorted(coincidencies, key=len)[0]

    return None


def mostrar_cas(cas, mostrar_nom=True):
    """Retorna 'nom comu (CAS)' o nomes el CAS si el nom no es conegut."""
    if mostrar_nom and cas in CAS_A_NOM:
        return f"{CAS_A_NOM[cas]} ({cas})"
    return cas


def puntuacio_afinitat(ki):

    return round(50.0 * math.exp(-ki / 8.0), 2)


def puntuacio_especificitat(ki_objectiu, kis_competidors):

    num_competidors = len(kis_competidors)
    punts = 35.0 - num_competidors * 3.0
    if num_competidors > 0:
        mitjana_comp = sum(kis_competidors) / num_competidors
        ratio = mitjana_comp / ki_objectiu
        if ratio >= 5.0:
            punts += 5.0
        elif ratio < 1.0:
            punts -= 8.0
    return round(max(0.0, min(35.0, punts)), 2)


def puntuacio_fiabilitat(num_estudis):

    return round(15.0 * (1.0 - math.exp(-num_estudis / 4.0)), 2)


def puntuacio_cisteines(num_cys):

    return round(max(0.0, 5.0 - abs(num_cys - 6)), 2)


def calcular_puntuacio(ki_objectiu, kis_competidors, num_estudis, num_cys=6,
                       pes_aff=1.0, pes_spec=1.0, pes_rel=1.0):
    aff  = min(puntuacio_afinitat(ki_objectiu)                        * pes_aff,  50.0)
    spec = min(puntuacio_especificitat(ki_objectiu, kis_competidors)  * pes_spec, 35.0)
    rel  = min(puntuacio_fiabilitat(num_estudis)                      * pes_rel,  15.0)
    cys  = puntuacio_cisteines(num_cys)
    return {
        "afinitat":      round(aff,  1),
        "especificitat": round(spec, 1),
        "fiabilitat":    round(rel,  1),
        "bonus_cys":     round(cys,  1),
        "total":         round(aff + spec + rel + cys, 1),
    }


def etiqueta_rang(puntuacio):
    if puntuacio >= 80: return "EXCELLENT"
    if puntuacio >= 60: return "GOOD"
    if puntuacio >= 40: return "MODERATE"
    return "POOR"


def parsejar_ki(valor):

    if valor is None:
        return None
    text = str(valor).strip()
    if text in ("-", "", "nan", "NaN", "N/A", "NA"):
        return None
    try:
        num = float(text)
        return num if num > 0 else None
    except ValueError:
        return None


def carregar_matriu_binding(ruta):

    llista_vocs = []
    llista_obps = []
    matriu      = {}

    with open(ruta, newline="", encoding="utf-8-sig") as fitxer:
        lector      = csv.reader(fitxer)
        capçalera   = next(lector)
        llista_obps = [col.strip() for col in capçalera[1:]]

        for fila in lector:
            if not fila or not fila[0].strip():
                continue
            voc = fila[0].strip()
            llista_vocs.append(voc)
            matriu[voc] = {}
            for i, obp in enumerate(llista_obps):
                valor = fila[i + 1] if i + 1 < len(fila) else None
                matriu[voc][obp] = parsejar_ki(valor)

    return llista_vocs, llista_obps, matriu


def analitzar(voc_objectiu, llista_vocs, llista_obps, matriu,
              pes_aff, pes_spec, pes_rel, mapa_cys):

    voc_trobat = resoldre_a_cas(voc_objectiu, llista_vocs)
    if voc_trobat is None:
        return None, None

    resultats = []
    for obp in llista_obps:
        ki_objectiu = matriu[voc_trobat].get(obp)
        if ki_objectiu is None:
            continue

        kis_competidors = [
            matriu[voc][obp]
            for voc in llista_vocs
            if voc != voc_trobat and matriu[voc].get(obp) is not None
        ]

        num_estudis = sum(1 for voc in llista_vocs if matriu[voc].get(obp) is not None)
        num_cys     = mapa_cys.get(obp, 6)

        punts        = calcular_puntuacio(ki_objectiu, kis_competidors, num_estudis, num_cys,
                                          pes_aff, pes_spec, pes_rel)
        mitjana_comp = sum(kis_competidors) / len(kis_competidors) if kis_competidors else None
        ratio_sel    = round(mitjana_comp / ki_objectiu, 1) if mitjana_comp is not None else None

        resultats.append({
            "obp":               obp,
            "ki_objectiu":       ki_objectiu,
            "num_competidors":   len(kis_competidors),
            "mitjana_comp_ki":   round(mitjana_comp, 2) if mitjana_comp else None,
            "ratio_selectivitat":ratio_sel,
            "num_estudis":       num_estudis,
            "num_cys":           num_cys,
            "punt_afinitat":     punts["afinitat"],
            "punt_especificitat":punts["especificitat"],
            "punt_fiabilitat":   punts["fiabilitat"],
            "punt_cys":          punts["bonus_cys"],
            "punt_total":        punts["total"],
            "rang":              etiqueta_rang(punts["total"]),
        })

    resultats.sort(key=lambda r: r["punt_total"], reverse=True)
    return resultats, voc_trobat


def mostrar_resultats(resultats, voc_trobat, top_n, detallat):
    num_mostrar = min(top_n, len(resultats))
    text_display = mostrar_cas(voc_trobat)
    print()
    print(f"  Target VOC : {text_display}")
    print(f"  OBPs found : {len(resultats)} with binding data")
    print(f"  Showing top: {num_mostrar}")
    print()

    amplades_col = [4, 16, 10, 7, 10, 7, 7, 10, 10, 10, 7, 8, 12]
    capçaleres   = ["#", "OBP", "Ki(uM)", "N_comp", "AvgComp", "Sel.x", "N_Cys",
                    "Aff/50", "Spec/35", "Rel/15", "Cys/5", "SCORE", "RANK"]

    print("  " + "  ".join(h.ljust(amplades_col[i]) for i, h in enumerate(capçaleres)))
    print("  " + "-" * (sum(amplades_col) + len(amplades_col) * 2))

    for i, r in enumerate(resultats[:num_mostrar], start=1):
        sel       = f"{r['ratio_selectivitat']}x" if r["ratio_selectivitat"] else "--"
        mitjana_c = f"{r['mitjana_comp_ki']:.2f}" if r["mitjana_comp_ki"] else "--"
        valors_fila = [
            str(i), r["obp"], str(r["ki_objectiu"]),
            str(r["num_competidors"]), mitjana_c, sel, str(r["num_cys"]),
            str(r["punt_afinitat"]), str(r["punt_especificitat"]),
            str(r["punt_fiabilitat"]), str(r["punt_cys"]),
            str(r["punt_total"]), r["rang"],
        ]
        linia = "  " + "  ".join(v.ljust(amplades_col[j]) for j, v in enumerate(valors_fila))
        print(linia)

    print()
    num_excellents = sum(1 for r in resultats if r["rang"] == "EXCELLENT")
    num_bons       = sum(1 for r in resultats if r["rang"] == "GOOD")


def parsejar_arguments():
    parser = argparse.ArgumentParser(
        prog="scentinel.py",
        description="Scentinel-code -- Rank OBPs as biosensor candidates for a target VOC.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--csv",    required=True,
                        help="Path to Compound_OBP_binding.csv (iOBPdb)")
    parser.add_argument("--voc",
                        help="Target VOC: common name, partial name, or CAS number")
    parser.add_argument("--mode",   default="automatic",
                        choices=["automatic", "custom"],
                        help="Scoring mode (default: automatic)")
    parser.add_argument("--w-aff",  type=float, default=1.0, help="Affinity weight [custom]")
    parser.add_argument("--w-spec", type=float, default=1.0, help="Specificity weight [custom]")
    parser.add_argument("--w-rel",  type=float, default=1.0, help="Reliability weight [custom]")

    return parser.parse_args()


def main():
    arguments = parsejar_arguments()

    print()
    print("  Scentinel-code . iGEM URV 2025")
    print("  OBP Biosensor Design Tool")

    # Carregar la matriu
    ruta_csv = Path(arguments.csv)
    if not ruta_csv.exists():
        print(f"\n  ERROR: file not found: {arguments.csv}\n")
        sys.exit(1)

    print(f"  Loading : {ruta_csv.name}")
    try:
        llista_vocs, llista_obps, matriu = carregar_matriu_binding(arguments.csv)
    except Exception as error:
        print(f"\n  ERROR reading CSV: {error}\n")
        sys.exit(1)
    print(f"  Dataset : {len(llista_vocs)} VOCs x {len(llista_obps)} OBPs")

    if not arguments.voc:
        print("  ERROR: --voc is required.\n")
        sys.exit(1)

    # Pesos
    if arguments.mode == "automatic":
        pes_aff = pes_spec = pes_rel = 1.0
    else:
        pes_aff, pes_spec, pes_rel = arguments.w_aff, arguments.w_spec, arguments.w_rel
        print(f"  Mode    : CUSTOM  (w_aff={pes_aff}  w_spec={pes_spec}  w_rel={pes_rel})")

    # Executar
    resultats, voc_trobat = analitzar(
        voc_objectiu=arguments.voc,
        llista_vocs=llista_vocs,
        llista_obps=llista_obps,
        matriu=matriu,
        pes_aff=pes_aff, pes_spec=pes_spec, pes_rel=pes_rel,
        mapa_cys={},
    )

    if resultats is None:
        print(f"  ERROR: VOC '{arguments.voc}' not found in dataset.")
        print(f"  Tip: check the compound name or try the CAS number directly.")
        print(f"  You can type common names (e.g. 'hexanal') or CAS numbers (e.g. '66-25-1').\n")
        sys.exit(1)

    if not resultats:
        print(f"  No OBPs have binding data for '{mostrar_cas(voc_trobat)}'.\n")
        sys.exit(0)

    mostrar_resultats(resultats, voc_trobat, top_n=15, detallat=False)


if __name__ == "__main__":
    main()
