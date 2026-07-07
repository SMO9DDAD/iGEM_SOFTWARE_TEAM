
import os, re, sys, time, subprocess, shutil
import requests
import pandas as pd
from pathlib import Path

# ── CONFIG ────────────────────────────────────────────────────────────────────

TOP_MUTANTS_PER_OBP = 5
BINDING_CUTOFF_A    = 4.0
DDG_THRESHOLD       = -0.5
EXHAUSTIVENESS      = 16
NUM_MODES           = 9
BOX_SIZE            = 34     



FOLDX_EXE   = str(Path(__file__).parent / "foldx" / "foldx")
FOLDX_WORK  = "/tmp/foldx_work"          # directori sense espais per FoldX
GNINA_EXE   = str(Path(__file__).parent / "gnina" / "gnina")
DOCKING_DIR = str(Path(__file__).parent / "docking")
RESULTS_DIR = str(Path(__file__).parent / "results")
MUTANT_DIR  = str(Path(__file__).parent / "mutants")

os.makedirs(DOCKING_DIR,  exist_ok=True)
os.makedirs(RESULTS_DIR,  exist_ok=True)
os.makedirs(MUTANT_DIR,   exist_ok=True)
os.makedirs(FOLDX_WORK,   exist_ok=True)

# Copiem rotabase.txt i altres fitxers de dades de FoldX al directori de treball
FOLDX_SRC = str(Path(__file__).parent / "foldx")
for _f in Path(FOLDX_SRC).glob("*.txt"):
    _dst = Path(FOLDX_WORK) / _f.name
    if not _dst.exists():
        shutil.copy(str(_f), str(_dst))

# ── UTILS ─────────────────────────────────────────────────────────────────────

def log(msg, indent=0):
    print("  " * indent + msg)

def sep(title="", w=62):
    if title:
        p = max(2, (w - len(title) - 2) // 2)
        print(f"\n{'='*p} {title} {'='*p}")
    else:
        print("\n" + "="*w)

def run_foldx(args):
    """Executa FoldX sempre des de /tmp/foldx_work (sense espais a la ruta)."""
    return subprocess.run(
        [FOLDX_EXE] + args,
        capture_output=True, text=True,
        cwd=FOLDX_WORK, timeout=1800,
    )

# FASE A: Binding site 

THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
    "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
    "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
    "TYR":"Y","VAL":"V",
}

def parse_pdb_atoms(pdb_path):
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                atoms.append({
                    "record":  rec,
                    "resname": line[17:20].strip(),
                    "chain":   line[21].strip(),
                    "resseq":  int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                })
            except:
                pass
    return atoms

def pdbqt_ligand_coords(pdbqt_path):
    coords, in1 = [], True
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL") and line.strip() != "MODEL 1":
                in1 = False
            if not in1:
                continue
            if line[:6].strip() in ("ATOM", "HETATM"):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except:
                    pass
    return coords

def dist3d(a, b):
    return sum((a[i]-b[i])**2 for i in range(3)) ** 0.5

def find_binding_residues(receptor_pdb, ligand_pdbqt):
    ratoms  = [a for a in parse_pdb_atoms(receptor_pdb) if a["record"] == "ATOM"]
    lcoords = pdbqt_ligand_coords(ligand_pdbqt)
    if not lcoords:
        log("AVIS: no s'han llegit coordenades del lligand", 2)
        return []
    close = set()
    for ra in ratoms:
        p = (ra["x"], ra["y"], ra["z"])
        for lc in lcoords:
            if dist3d(p, lc) <= BINDING_CUTOFF_A:
                one = THREE_TO_ONE.get(ra["resname"], "X")
                close.add((ra["chain"], ra["resseq"], ra["resname"], one))
                break
    return sorted(close, key=lambda x: (x[0], x[1]))

#  FASE B: FoldX 
def foldx_repair(pdb_path, work_dir):
    pdb_abs  = Path(pdb_path).resolve()
    stem     = pdb_abs.stem
    out_abs  = Path(work_dir).resolve()
    repaired = out_abs / f"{stem}_Repair.pdb"

    if repaired.exists():
        log(f"RepairPDB ja existeix: {repaired.name}", 1)
        return str(repaired)

    # Copiem el PDB a /tmp/foldx_work (sense espais)
    tmp_pdb = Path(FOLDX_WORK) / pdb_abs.name
    shutil.copy(str(pdb_abs), str(tmp_pdb))

    res = run_foldx([
        "--command=RepairPDB",
        f"--pdb={pdb_abs.name}",
    ])

    tmp_pdb.unlink(missing_ok=True)

    # FoldX escriu <stem>_Repair.pdb al seu CWD (/tmp/foldx_work)
    repaired_tmp = Path(FOLDX_WORK) / f"{stem}_Repair.pdb"
    if repaired_tmp.exists():
        shutil.move(str(repaired_tmp), str(repaired))
        return str(repaired)

    log(f"RepairPDB error:\n{res.stdout[-500:]}", 1)
    return None


def foldx_position_scan(repaired_pdb, binding_residues, work_dir, obp_name):
    sep(f"PositionScan — {obp_name}")

    rep_abs = Path(repaired_pdb).resolve()
    stem    = rep_abs.stem
    out_abs = Path(work_dir).resolve()

    # Auto-detecció: comprova si ja existeixen fitxers PS_ per tots els residus
    existing = list(out_abs.glob(f"PS_{stem}*.fxout")) + \
               list(out_abs.glob(f"PS_{stem}*.txt"))

    residus_a_escanar = [
        (chain, resseq, resname, wt_one)
        for chain, resseq, resname, wt_one in binding_residues
        if wt_one not in ("C", "X")
    ]

    residus_amb_fitxer = set()
    for f in existing:
        tag = f.stem[len(f"PS_{stem}_"):]
        residus_amb_fitxer.add(tag)

    tots_coberts = all(
        f"{wt_one}{chain}{resseq}" in residus_amb_fitxer
        for chain, resseq, resname, wt_one in residus_a_escanar
    )

    if existing and tots_coberts:
        log(f"Scan previ complet trobat ({len(existing)} fitxers) — reutilitzant", 1)
        mutations_found = _parse_ps_files(existing, binding_residues, obp_name)
        return _filter_mutations(mutations_found)
    elif existing and not tots_coberts:
        log(f"Scan previ incomplet — executant FoldX per als residus que falten", 1)
    else:
        log(f"Cap scan previ — executant FoldX (15-30 min)...", 1)

    # Execució normal: copiem a /tmp/foldx_work i executem
    tmp_pdb = Path(FOLDX_WORK) / rep_abs.name
    if not tmp_pdb.exists():
        shutil.copy(str(rep_abs), str(tmp_pdb))

    mutations_found = []

    for chain, resseq, resname, wt_one in binding_residues:
        if wt_one in ("C", "X"):
            log(f"Saltant {chain}{resseq}{resname} ({'Cys' if wt_one=='C' else 'no standard'})", 1)
            continue

        scan_code = f"{wt_one}{chain}{resseq}a"
        log(f"Scan {chain}{resseq}{resname}({wt_one}) → {scan_code} ...", 1)

        run_foldx([
            "--command=PositionScan",
            f"--pdb={rep_abs.name}",
            f"--positions={scan_code}",
        ])

        # Movem fitxers PS_ generats a out_abs reanomenant-los amb el residu
        # per evitar que cada scan sobreescrigui el fitxer de l'anterior
        res_tag   = f"{wt_one}{chain}{resseq}"
        res_files = []
        for ps_orig in list(Path(FOLDX_WORK).glob(f"PS_{stem}*.fxout")) + \
                       list(Path(FOLDX_WORK).glob(f"PS_{stem}*.txt")):
            suffix = ps_orig.suffix
            dest   = out_abs / f"PS_{stem}_{res_tag}{suffix}"
            shutil.move(str(ps_orig), str(dest))
            res_files.append(dest)

        # Parsegem només els fitxers d'aquest residu
        new = _parse_ps_files(res_files, [(chain, resseq, resname, wt_one)], obp_name)
        mutations_found.extend(new)

    tmp_pdb.unlink(missing_ok=True)

    return _filter_mutations(mutations_found)


def _parse_ps_files(ps_files, binding_residues, obp_name):
    """Parseja fitxers PS_ de FoldX i retorna llista de mutacions."""
    mutations = []
    res_lookup = {n: (c, o) for c, n, r, o in binding_residues}

    for ps_file in ps_files:
        with open(ps_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(("Pdb", "#", "total")):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                try:
                    ddg = float(parts[1])
                    m   = re.search(r'[A-Z][A-Z](\d+)([A-Z])', parts[0])
                    if not m:
                        continue
                    resseq = int(m.group(1))
                    mut_aa = m.group(2)
                    if resseq not in res_lookup:
                        continue
                    chain, wt_aa = res_lookup[resseq]
                    if mut_aa in (wt_aa, "C"):
                        continue
                    mutations.append({
                        "obp_name":  obp_name,
                        "chain":     chain,
                        "resseq":    resseq,
                        "wt_aa":     wt_aa,
                        "mut_aa":    mut_aa,
                        "mutation":  f"{wt_aa}{resseq}{mut_aa}",
                        "ddg_foldx": ddg,
                    })
                except (ValueError, IndexError):
                    continue
    return mutations


def _filter_mutations(mutations_found):
    """Filtra per DDG_THRESHOLD, o agafa les 5 millors si no n'hi ha cap."""
    good = [m for m in mutations_found if m["ddg_foldx"] < DDG_THRESHOLD]
    good = sorted(good, key=lambda x: x["ddg_foldx"])
    log(f"{len(mutations_found)} escanejades, {len(good)} passen ΔΔG<{DDG_THRESHOLD}", 1)
    if not good and mutations_found:
        log(f"Relaxant filtre: agafant les 5 millors per ΔΔG", 1)
        good = sorted(mutations_found, key=lambda x: x["ddg_foldx"])[:5]
    return good


def foldx_build_mutant(repaired_pdb, mut, work_dir):
    rep_abs   = Path(repaired_pdb).resolve()
    stem      = rep_abs.stem
    out_abs   = Path(work_dir).resolve()
    mut_label = mut["mutation"]
    mut_code  = f"{mut['wt_aa']}{mut['chain']}{mut['resseq']}{mut['mut_aa']}"

    # individual_list.txt ha d'estar a /tmp/foldx_work
    ilist = Path(FOLDX_WORK) / "individual_list.txt"
    ilist.write_text(mut_code + ";\n")

    tmp_pdb = Path(FOLDX_WORK) / rep_abs.name
    if not tmp_pdb.exists():
        shutil.copy(str(rep_abs), str(tmp_pdb))

    res = run_foldx([
        "--command=BuildModel",
        f"--pdb={rep_abs.name}",
        "--mutant-file=individual_list.txt",
        "--numberOfRuns=1",
    ])

    ilist.unlink(missing_ok=True)
    tmp_pdb.unlink(missing_ok=True)

    # FoldX escriu <stem>_1.pdb a /tmp/foldx_work
    final = out_abs / f"{stem}_{mut_label}.pdb"
    candidate = Path(FOLDX_WORK) / f"{stem}_1.pdb"
    if candidate.exists():
        shutil.move(str(candidate), str(final))
        for wt_file in Path(FOLDX_WORK).glob(f"{stem}_WT_*.pdb"):
            wt_file.unlink()
        return str(final)

    log(f"BuildModel no ha generat PDB. stdout: {res.stdout[-300:]}", 2)
    return None

# FASE D: Receptor PDBQT 
def fix_pdb_elements(pdb_in, pdb_out):
    """Afegeix la columna d'element si falta.
    FoldX sovint deixa buida la columna 77-78 i RDKit/Meeko peta amb
    'Element not found'."""
    with open(pdb_in) as f:
        lines = f.readlines()
    fixed = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            if len(line) < 78 or line[76:78].strip() == "":
                atom_name = line[12:16].strip()
                elem = re.sub(r'[0-9]', '', atom_name).strip()
                elem = elem[0].upper() if elem else "C"
                line = line.rstrip('\n').ljust(78)
                line = line[:76] + f"{elem:>2}" + "\n"
        fixed.append(line)
    with open(pdb_out, "w") as f:
        f.writelines(fixed)

def prepare_receptor_pdbqt(pdb, pdbqt):
    if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 0:
        return True

    tmp_pdb   = "/tmp/receptor_tmp.pdb"
    tmp_fixed = "/tmp/receptor_fixed.pdb"
    tmp_pdbqt = "/tmp/receptor_tmp.pdbqt"

    shutil.copy(pdb, tmp_pdb)
    fix_pdb_elements(tmp_pdb, tmp_fixed)

    r = subprocess.run(
        ["mk_prepare_receptor.py", "--read_pdb", tmp_fixed, "-p", tmp_pdbqt],
        capture_output=True, text=True
    )
    if os.path.exists(tmp_pdbqt) and os.path.getsize(tmp_pdbqt) > 0:
        shutil.copy(tmp_pdbqt, pdbqt)
        return True

    print(f"    meeko error: {r.stderr[-300:]}")
    return False

#  FASE D: Docking 

def pose_center(pose_pdbqt):
    xs, ys, zs = [], [], []
    with open(pose_pdbqt) as f:
        for line in f:
            if line.startswith("MODEL") and line.strip() != "MODEL 1":
                break
            if line[:6].strip() in ("ATOM", "HETATM"):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except:
                    pass
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)) if xs else (0, 0, 0)

def run_gnina(receptor, ligand, out, cx, cy, cz):
    """Docking pur Vina — sense CNN scoring."""
    cmd = [
        GNINA_EXE, "--receptor", receptor, "--ligand", ligand, "--out", out,
        "--center_x", str(round(cx, 2)), "--center_y", str(round(cy, 2)), "--center_z", str(round(cz, 2)),
        "--size_x", str(BOX_SIZE), "--size_y", str(BOX_SIZE), "--size_z", str(BOX_SIZE),
        "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
        "--no_gpu", "--cnn_scoring", "none", "--scoring", "vina",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    vina = None
    for line in res.stdout.split("\n"):
        line = line.strip()
        if re.match(r"^1\s", line):
            parts = line.split()
            try:
                vina = float(parts[1])
            except:
                pass
            break
    if vina is None and res.returncode != 0:
        print(f"    GNINA stderr: {res.stderr[-200:]}")
    return vina



#  SMILES 

def resolve_smiles(voc_name):
    """Obté SMILES via PubChem. Si falla, atura el programa amb error."""
    try:
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
               f"{voc_name}/property/IsomericSMILES/JSON")
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            smiles = r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
            if smiles:
                return smiles
        raise ValueError(f"PubChem HTTP {r.status_code}")
    except Exception as e:
        log(f"ERROR: no s'ha pogut obtenir el SMILES de '{voc_name}' ({e})")
        log(f"Comprova el nom del VOC o la connexio a internet.")
        sys.exit(1)

#  ESTRUCTURA WT 

IOBPDB_PDB_URL = "https://raw.githubusercontent.com/sshuklz/iobpdb_app/master/alpha_pdbs/{}.pdb"
P2RANK_EXE     = str(Path(__file__).parent / "p2rank" / "prank")
P2RANK_MIN_SCORE = 10
BOX_WITH_POCKET  = 26

def _download_iobpdb(obp_name, out_pdb):
    try:
        r = requests.get(IOBPDB_PDB_URL.format(obp_name), timeout=30)
        if r.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(r.text)
            log("Estructura descarregada d'iobpdb", 2)
            return True
        log(f"No trobat a iobpdb (HTTP {r.status_code})", 2)
        return False
    except Exception as e:
        log(f"Error descarregant: {e}", 2)
        return False

def _esmfold(obp_name, out_pdb):
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?query={obp_name}&format=json&size=1",
            timeout=10)
        if r.status_code != 200:
            return False
        results = r.json().get("results", [])
        if not results:
            return False
        uid = results[0]["primaryAccession"]
        r2  = requests.get(f"https://www.uniprot.org/uniprot/{uid}.fasta", timeout=15)
        if r2.status_code != 200:
            return False
        seq = "".join(r2.text.strip().split("\n")[1:])
        if len(seq) > 400:
            seq = seq[20:]
        log(f"Predint estructura ESMFold ({len(seq)} aa)...", 2)
        for i in range(3):
            r3 = requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/",
                               data=seq, timeout=180)
            if r3.status_code == 200:
                with open(out_pdb, "w") as f:
                    f.write(r3.text)
                log("Estructura predita per ESMFold", 2)
                return True
            log(f"ESMFold intent {i+1} error {r3.status_code}", 2)
            time.sleep(10)
    except Exception as e:
        log(f"ESMFold error: {e}", 2)
    return False

def _get_structure(obp_name):
    pdb = str(Path(DOCKING_DIR) / f"{obp_name}.pdb")
    if os.path.exists(pdb):
        return pdb
    log("Buscant estructura a iobpdb...", 2)
    if _download_iobpdb(obp_name, pdb):
        return pdb
    log("Intentant ESMFold...", 2)
    if _esmfold(obp_name, pdb):
        return pdb
    log("No s'ha pogut obtenir l'estructura", 2)
    return None

def _get_center(pdb):
    xs, ys, zs = [], [], []
    with open(pdb) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except:
                    pass
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)

def _detect_pocket(pdb):
    if not os.path.exists(P2RANK_EXE):
        return None
    pdb_abs = str(Path(pdb).resolve())
    out_dir = str(Path(pdb + "_p2rank").resolve())
    p2rank_d = str(Path(P2RANK_EXE).parent.resolve())
    os.makedirs(out_dir, exist_ok=True)
    try:
        res = subprocess.run(
            ["bash", P2RANK_EXE, "predict", "-f", pdb_abs, "-o", out_dir, "-threads", "1"],
            capture_output=True, text=True, timeout=120, cwd=p2rank_d)
        if res.returncode != 0:
            return None
        csvs = list(Path(out_dir).rglob("*predictions.csv"))
        if not csvs:
            return None
        df_p = pd.read_csv(csvs[0])
        df_p.columns = df_p.columns.str.strip()
        if df_p.empty:
            return None
        best = df_p.iloc[0]
        score = float(best.get("score", 0))
        if score < P2RANK_MIN_SCORE:
            return None
        return float(best["center_x"]), float(best["center_y"]), float(best["center_z"])
    except Exception as e:
        log(f"P2Rank error: {e}", 2)
        return None

def _prepare_ligand(smiles, out_pdbqt):
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

def _prepare_receptor_obabel(pdb, pdbqt):
    os.system(f'obabel "{pdb}" -O "{pdbqt}" -xr 2>/dev/null')

def _dock_wt(obp_name, voc_name, voc_smiles):
    """Fa docking del wildtype i retorna vina_score. None si falla."""
    slug      = voc_name.replace(" ","_").replace("-","_")
    pdb       = str(Path(DOCKING_DIR) / f"{obp_name}.pdb")
    pdbqt_rec = str(Path(DOCKING_DIR) / f"{obp_name}.pdbqt")
    lig_pdbqt = str(Path(DOCKING_DIR) / f"ligand_{slug}.pdbqt")
    pose_out  = str(Path(DOCKING_DIR) / f"{obp_name}_{slug}_p1_out.pdbqt")

    pdb = _get_structure(obp_name)
    if not pdb:
        return None

    if not os.path.exists(lig_pdbqt):
        log("Preparant lligand...", 2)
        _prepare_ligand(voc_smiles, lig_pdbqt)

    _prepare_receptor_obabel(pdb, pdbqt_rec)

    pocket = _detect_pocket(pdb)
    if pocket:
        cx, cy, cz = pocket
        box = BOX_WITH_POCKET
        log(f"Grid al pocket ({box} Å)", 2)
    else:
        cx, cy, cz = _get_center(pdb)
        box = BOX_SIZE
        log(f"Grid al centroide ({box} Å)", 2)

    vina = run_gnina(pdbqt_rec, lig_pdbqt, pose_out, cx, cy, cz)
    if vina is not None:
        log(f"Vina WT: {vina:.3f} kcal/mol", 2)
    return vina

# ── CARREGA WT (CSV o docking automàtic) ───────────────────────────────────────

def load_wt(voc_name, obp_names_needed=None):

    slug = voc_name.replace(" ", "_").replace("-", "_")
    csv  = Path(RESULTS_DIR) / f"docking_complet_{slug}.csv"

    wt_all = {}
    if csv.exists():
        df = pd.read_csv(csv)
        score_col = next((c for c in ("vina_score","docking_score_vina") if c in df.columns), None)
        kd_col    = next((c for c in ("kd_nm_exp","kd_exp") if c in df.columns), None)
        for _, r in df.iterrows():
            wt_all[r["obp_name"]] = {
                "vina_wt": r.get(score_col) if score_col else None,
                "kd_exp":  r.get(kd_col)    if kd_col    else None,
            }

    if obp_names_needed:
        missing = [n for n in obp_names_needed if n not in wt_all or wt_all[n]["vina_wt"] is None]
        if missing:
            log(f"OBPs sense score WT al CSV: {missing} — fent docking automàtic...", 1)
            voc_smiles = resolve_smiles(voc_name)
            for obp_name in missing:
                log(f"Docking WT de {obp_name}...", 1)
                vina = _dock_wt(obp_name, voc_name, voc_smiles)
                wt_all[obp_name] = {"vina_wt": vina, "kd_exp": wt_all.get(obp_name, {}).get("kd_exp")}

    return wt_all

def run_mutant_pipeline(obp_name, voc_name, wt_data):
    slug       = voc_name.replace(" ", "_").replace("-", "_")
    pdb_wt     = str(Path(DOCKING_DIR) / f"{obp_name}.pdb")
    lig_pdbqt  = str(Path(DOCKING_DIR) / f"ligand_{slug}.pdbqt")
    pose_pdbqt = str(Path(DOCKING_DIR) / f"{obp_name}_{slug}_p1_out.pdbqt")
    obp_mutdir = str(Path(MUTANT_DIR)  / obp_name)
    os.makedirs(obp_mutdir, exist_ok=True)

    vina_wt = wt_data.get("vina_wt")

    sep(f"OBP: {obp_name}")
    vina_str = f"{vina_wt:.3f}" if vina_wt is not None else "N/A"
    log(f"Wildtype: Vina={vina_str} kcal/mol", 1)

    for path, desc in [(pdb_wt, "PDB WT"), (lig_pdbqt, "Lligand"), (pose_pdbqt, "Pose GNINA")]:
        if not os.path.exists(path):
            log(f"No trobo {desc}: {path}", 1)
            return []

    # A — Binding site
    log("Fase A — Binding site...", 1)
    bres = find_binding_residues(pdb_wt, pose_pdbqt)
    if not bres:
        log("Sense residus al binding site.", 1)
        return []
    log(f"  {len(bres)} residus: {[f'{c}{n}{r}({o})' for c,n,r,o in bres]}", 1)

    # B — Repair
    log("Fase B — RepairPDB (2-5 min)...", 1)
    repaired = foldx_repair(pdb_wt, obp_mutdir)
    if not repaired:
        log("RepairPDB ha fallat.", 1)
        return []
    log(f"  OK: {Path(repaired).name}", 1)

    # B — PositionScan
    log("Fase B — PositionScan (15-30 min)...", 1)
    mutations = foldx_position_scan(repaired, bres, obp_mutdir, obp_name)
    if not mutations:
        log("Cap mutació. Prova DDG_THRESHOLD=0.0", 1)
        return []

    top        = mutations[:TOP_MUTANTS_PER_OBP]
    cx, cy, cz = pose_center(pose_pdbqt)
    results    = []

    for mut in top:
        ml = mut["mutation"]
        log(f"\nFase C — BuildModel {ml} (ΔΔG={mut['ddg_foldx']:.3f})...", 1)
        mpdb = foldx_build_mutant(repaired, mut, obp_mutdir)
        if not mpdb:
            continue

        mpdbqt = mpdb.replace(".pdb", ".pdbqt")
        if not prepare_receptor_pdbqt(mpdb, mpdbqt):
            log(f"  prepare_receptor ha fallat", 2)
            continue

        log(f"Fase D — Docking {ml}...", 1)
        out_dock = str(Path(obp_mutdir) / f"{obp_name}_{ml}_{slug}_out.pdbqt")
        vina = run_gnina(mpdbqt, lig_pdbqt, out_dock, cx, cy, cz)
        if vina is None:
            log(f"  Docking sense resultat", 2)
            continue

        dvina = (vina - vina_wt) if (vina is not None and vina_wt is not None) else None
        dvina_str = f"{dvina:+.3f}" if dvina is not None else "N/A"
        log(f"  Vina={vina:.3f} kcal/mol (Δ={dvina_str})", 2)

        results.append({
            "obp_name":    obp_name,    "mutation":     ml,
            "chain":       mut["chain"],"resseq":       mut["resseq"],
            "wt_aa":       mut["wt_aa"],"mut_aa":       mut["mut_aa"],
            "ddg_foldx":   mut["ddg_foldx"],
            "vina_wt":     vina_wt,    "vina_mutant":  vina,  "delta_vina": dvina,
            "kd_exp":      wt_data.get("kd_exp"),
        })

    return results

def mostrar_ranking(all_results, voc_name):
    if not all_results:
        log("\nSense resultats de mutagènesi.")
        return
    df   = pd.DataFrame(all_results)
    slug = voc_name.replace(" ", "_").replace("-", "_")
    sep("RANKING FINAL MUTANTS")
    log("Ordenat per ΔVina (negatiu = millora d'afinitat)\n")
    df = df.sort_values("delta_vina", ascending=True).reset_index(drop=True)
    log(f"  {'#':<4} {'OBP':<18} {'Mut':<8} {'ΔΔG FoldX':>10} {'Vina_wt':>8} {'Vina_mut':>8} {'ΔVina':>8}")
    log("  " + "─"*65)
    for i, r in df.iterrows():
        s = lambda v, fmt: (fmt % v) if v is not None else "  N/A"
        log(f"  {i+1:<4} {r['obp_name']:<18} {r['mutation']:<8} "
            f"{r['ddg_foldx']:>10.3f} {s(r['vina_wt'],'%8.3f')} {s(r['vina_mutant'],'%8.3f')} "
            f"{s(r['delta_vina'],'%+8.3f')}")

    out = f"{RESULTS_DIR}/mutants_{slug}.csv"
    df.to_csv(out, index=False)
    log(f"\n  Guardat: {out}")

    best = df.iloc[0]
    sep("MILLOR MUTANT")
    log(f"OBP:       {best['obp_name']}")
    log(f"Mutació:   {best['mutation']}  ({best['wt_aa']}{best['resseq']} → {best['mut_aa']})")
    log(f"ΔΔG FoldX: {best['ddg_foldx']:.3f} kcal/mol")
    if best["delta_vina"] is not None:
        log(f"ΔVina:     {best['delta_vina']:+.3f} kcal/mol  ({best['vina_wt']:.3f} → {best['vina_mutant']:.3f})")
    sep()

    interpretacio(df, voc_name)


def interpretacio(df, voc_name):
    """Interpreta els resultats de mutagènesi (Vina) i dona recomanacions."""
    sep("INTERPRETACIÓ")

    VINA_MILLORA    = -0.3   # ΔVina per considerar millora real (kcal/mol)
    VINA_ACCEPTABLE =  0.3   # ΔVina màxim per considerar "quasi igual"
    DDG_ESTABLE     = -1.0   # ΔΔG prou negatiu per considerar estabilització

    milloren    = df[df["delta_vina"] <= VINA_MILLORA]
    acceptables = df[(df["delta_vina"] > VINA_MILLORA) & (df["delta_vina"] <= VINA_ACCEPTABLE)]
    estables    = df[df["ddg_foldx"] <= DDG_ESTABLE]

    # ── Cas 1: hi ha mutacions que milloren l'afinitat Vina ─────────────────
    if not milloren.empty:
        log(f"TROBADES {len(milloren)} MUTACIÓ(NS) QUE MILLOREN L'AFINITAT VINA:\n")
        for _, r in milloren.iterrows():
            log(f"  {r['mutation']:<8}  ΔVina={r['delta_vina']:+.3f}  ΔΔG={r['ddg_foldx']:+.3f}", 1)
            if r["ddg_foldx"] <= DDG_ESTABLE:
                log(f"  → Proteïna més estable (ΔΔG={r['ddg_foldx']:.3f}) i millor afinitat: candidat sòlid ✓✓", 1)
            elif r["ddg_foldx"] > 0:
                log(f"  → ΔΔG positiu: desestabilitza l'estructura. Millora d'afinitat però valida estabilitat", 1)
            else:
                log(f"  → Candidat a validar experimentalment ✓", 1)

    # ── Cas 2: cap millora Vina ──────────────────────────────────────────────
    else:
        log("Cap mutació millora l'afinitat Vina respecte al wildtype.\n")
        log(f"  Wildtype Vina = {df.iloc[0]['vina_wt']:.3f} kcal/mol  —  ja força optimitzat.\n")

        if not acceptables.empty:
            log(f"  Mutacions amb pèrdua mínima (ΔVina ≤ {VINA_ACCEPTABLE}):")
            for _, r in acceptables.iterrows():
                estab = f"  ΔΔG={r['ddg_foldx']:+.3f} (més estable)" if r["ddg_foldx"] <= DDG_ESTABLE else f"  ΔΔG={r['ddg_foldx']:+.3f}"
                log(f"    {r['mutation']:<8}  ΔVina={r['delta_vina']:+.3f}{estab}", 1)
            log("")

        if not estables.empty:
            log(f"  Mutacions que estabilitzen molt la proteïna (ΔΔG ≤ {DDG_ESTABLE}) tot i perdre afinitat:")
            for _, r in estables.sort_values("ddg_foldx").iterrows():
                log(f"    {r['mutation']:<8}  ΔΔG={r['ddg_foldx']:+.3f}  ΔVina={r['delta_vina']:+.3f}", 1)
            log(f"  → Interessants per aplicacions que requereixin proteïna estable (biosensors, etc.)\n")

    # ── Recomanacions ────────────────────────────────────────────────────────
    sep("RECOMANACIONS")

    if milloren.empty:
        log("1. El wildtype sembla evolutivament optimitzat per aquest VOC.")
        log("   Considera provar altres OBPs (augmenta TOP_OBPS_FROM_MAIN3).\n")
        log("2. El pipeline usa la pose del wildtype com a referència per al docking")
        log("   dels mutants. Si la mutació reorienta el lligand, Vina pot")
        log("   subestimar la millora real. Validació experimental recomanada.\n")
        if not estables.empty:
            best_e = estables.sort_values("ddg_foldx").iloc[0]
            log(f"3. Candidat prioritari per validació: {best_e['mutation']}")
            log(f"   ΔΔG={best_e['ddg_foldx']:.3f} (molt estable)  ΔVina={best_e['delta_vina']:+.3f}")
            log(f"   Proteïna estable que uneix quasi igual → bon candidat pràctic.")
    else:
        best_m = milloren.sort_values("delta_vina", ascending=True).iloc[0]
        log(f"1. Candidat prioritari: {best_m['mutation']}  (ΔVina={best_m['delta_vina']:+.3f} kcal/mol)")
        log(f"   Valida experimentalment amb fluorescència de desplaçament de lligand.")
        if best_m["ddg_foldx"] <= DDG_ESTABLE:
            log(f"   Proteïna estable i millor afinitat → alta confiança.\n")
        else:
            log(f"   Comprova l'estabilitat estructural (ΔΔG={best_m['ddg_foldx']:+.3f}).\n")
        log("2. Si vols explorar més opcions, augmenta TOP_MUTANTS_PER_OBP.")
    sep()

#  MAIN 
if __name__ == "__main__":
    sep("OBP MUTANT DESIGNER")


    VOC_NAME = input("VOC diana (ex: 1-octen-3-ol): ").strip()
    if not VOC_NAME:
        log("ERROR: cal introduir un VOC.")
        exit(1)


    print("\nMode de selecció d'OBPs:")
    print("  [1] Automàtic — agafa el top N del CSV de main3 (per defecte)")
    print("  [2] Manual    — especifica un OBP concret")
    mode = input("Tria [1/2] (per defecte 1): ").strip()

    slug     = VOC_NAME.replace(" ", "_").replace("-", "_")
    csv_path = Path(RESULTS_DIR) / f"docking_complet_{slug}.csv"

    if mode == "2":

        if csv_path.exists():
            df_disp = pd.read_csv(csv_path)
            score_col = next((c for c in ("vina_score","docking_score_vina") if c in df_disp.columns), None)
            if score_col:
                df_disp = df_disp.dropna(subset=[score_col]).sort_values(score_col)
            log(f"\nOBPs disponibles al CSV ({len(df_disp)} total):")
            for _, r in df_disp.head(15).iterrows():
                sc = f"{r[score_col]:.3f}" if score_col else "—"
                log(f"  {r['obp_name']:<25}  Vina={sc}  Kd_exp={r.get('kd_nm_exp', r.get('kd_exp','—'))} nM", 1)
            if len(df_disp) > 15:
                log(f"  ... i {len(df_disp)-15} més al CSV.", 1)

        obp_input = input("\nNom de l'OBP a millorar (ex: CpinOBP2): ").strip()
        if not obp_input:
            log("ERROR: cal introduir un nom d'OBP.")
            exit(1)
        obps_a_processar = [obp_input]
        log(f"\nMode manual → OBP: {obp_input}")

    else:

        top_n_str = input("Quants OBPs del top vols usar? [3]: ").strip()
        TOP_OBPS_FROM_MAIN3 = int(top_n_str) if top_n_str.isdigit() and int(top_n_str) > 0 else 3

        if not csv_path.exists():
            log(f"No trobo CSV de docking a {csv_path}.")
            log(f"S'hauran de dockar tots els OBPs seleccionats des de zero.")
            obp_input = input("Introdueix el nom de l'OBP manualment: ").strip()
            if not obp_input:
                log("ERROR: cal introduir un OBP.")
                exit(1)
            obps_a_processar = [obp_input]
        else:
            df_wt = pd.read_csv(csv_path)
            score_col = next((c for c in ("vina_score","docking_score_vina") if c in df_wt.columns), None)
            if score_col:
                df_wt = df_wt.dropna(subset=[score_col]).sort_values(score_col)
            df_wt = df_wt.head(TOP_OBPS_FROM_MAIN3)
            obps_a_processar = list(df_wt["obp_name"])
            log(f"\nMode automàtic → Top {TOP_OBPS_FROM_MAIN3} OBPs:")
            for _, r in df_wt.iterrows():
                sc = f"{r[score_col]:.3f}" if score_col else "—"
                log(f"  {r['obp_name']:<25}  Vina={sc}  Kd_exp={r.get('kd_nm_exp', r.get('kd_exp','—'))} nM", 1)


    wt_all = load_wt(VOC_NAME, obp_names_needed=obps_a_processar)

    sep("OBP MUTANT DESIGNER")
    log(f"VOC: {VOC_NAME}  |  OBPs: {obps_a_processar}  |  Top {TOP_MUTANTS_PER_OBP} mutants/OBP")
    sep()

    all_results = []
    for obp_name in obps_a_processar:
        all_results.extend(run_mutant_pipeline(obp_name, VOC_NAME, wt_all.get(obp_name, {})))

    mostrar_ranking(all_results, VOC_NAME)