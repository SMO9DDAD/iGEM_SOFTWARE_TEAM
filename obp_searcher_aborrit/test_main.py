"""
Tests unitaris per a main.py — OBP Finder (iGEM URV 2025)
Executa amb: python -m pytest test_main.py -v
"""

import math
import numpy as np
import pandas as pd
import pytest

# Importem les funcions que volem testar
from main import (
    BIG_KI_VALUE,
    DEFAULT_WEIGHTS,
    SELECTIVITY_TAU,
    convert_ki_to_float,
    find_voc_rows,
    read_interferent_file,
    compute_s1_affinity,
    compute_s2_selectivity,
    compute_s5_promiscuity,
    compute_final_score,
    build_obp_ranking,
)


# ──────────────────────────────────────────────
# 1. convert_ki_to_float
# ──────────────────────────────────────────────

class TestConvertKiToFloat:

    def test_valor_numeric_simple(self):
        assert convert_ki_to_float("3.5") == pytest.approx(3.5)

    def test_valor_numeric_enter(self):
        assert convert_ki_to_float("100") == pytest.approx(100.0)

    def test_valor_amb_signe_major(self):
        # ">10" → 10 * 1.1 = 11.0
        result = convert_ki_to_float(">10")
        assert result == pytest.approx(11.0)

    def test_valor_amb_signe_major_gran(self):
        # ">500" → 550.0
        result = convert_ki_to_float(">500")
        assert result == pytest.approx(550.0)

    def test_nan_retorna_nan(self):
        assert math.isnan(convert_ki_to_float(np.nan))

    def test_text_invalid_retorna_nan(self):
        assert math.isnan(convert_ki_to_float("no_es_un_numero"))

    def test_signe_major_invalit_retorna_big_ki(self):
        # ">abc" → no es pot extreure nombre → BIG_KI_VALUE
        result = convert_ki_to_float(">abc")
        assert result == pytest.approx(BIG_KI_VALUE)

    def test_espais_blancs(self):
        assert convert_ki_to_float("  7.2  ") == pytest.approx(7.2)

    def test_espai_no_separable(self):
        # '\xa0' és un espai no separable
        assert convert_ki_to_float("\xa0 5.0 \xa0") == pytest.approx(5.0)

    def test_float_ja_es_float(self):
        assert convert_ki_to_float(42.0) == pytest.approx(42.0)


# ──────────────────────────────────────────────
# 2. find_voc_rows
# ──────────────────────────────────────────────

class TestFindVocRows:

    @pytest.fixture
    def taula(self):
        return pd.DataFrame({
            'Nom': ['Ethanol', 'Ethyl acetate', 'Benzaldehyde', 'Limonene'],
            'CAS': ['64-17-5', '141-78-6', '100-52-7', '138-86-3'],
            'OBP1': [1.0, 2.0, 3.0, 4.0],
        })

    def test_troba_coincidencia_parcial(self, taula):
        result = find_voc_rows(taula, 'Nom', 'Ethyl')
        assert len(result) == 1
        assert result.iloc[0]['Nom'] == 'Ethyl acetate'

    def test_cerca_insensible_majuscules(self, taula):
        result = find_voc_rows(taula, 'Nom', 'ethanol')
        assert len(result) == 1

    def test_no_troba_res(self, taula):
        result = find_voc_rows(taula, 'Nom', 'xyzzy')
        assert result.empty

    def test_troba_multiples(self, taula):
        result = find_voc_rows(taula, 'Nom', 'eth')
        assert len(result) == 2   # Ethanol i Ethyl acetate

    def test_cerca_completa(self, taula):
        result = find_voc_rows(taula, 'Nom', 'Limonene')
        assert len(result) == 1


# ──────────────────────────────────────────────
# 3. read_interferent_file
# ──────────────────────────────────────────────

class TestReadInterferentFile:

    def test_llegeix_linies_valides(self, tmp_path):
        f = tmp_path / "interf.txt"
        f.write_text("Ethanol\nAcetone\nBenzaldehyde\n", encoding='utf-8')
        result = read_interferent_file(str(f))
        assert result == ['Ethanol', 'Acetone', 'Benzaldehyde']

    def test_ignora_comentaris(self, tmp_path):
        f = tmp_path / "interf.txt"
        f.write_text("# comentari\nEthanol\n# altre\nAcetone\n", encoding='utf-8')
        result = read_interferent_file(str(f))
        assert result == ['Ethanol', 'Acetone']

    def test_ignora_linies_buides(self, tmp_path):
        f = tmp_path / "interf.txt"
        f.write_text("\nEthanol\n\nAcetone\n\n", encoding='utf-8')
        result = read_interferent_file(str(f))
        assert result == ['Ethanol', 'Acetone']

    def test_fitxer_buit(self, tmp_path):
        f = tmp_path / "buit.txt"
        f.write_text("", encoding='utf-8')
        result = read_interferent_file(str(f))
        assert result == []

    def test_nomes_comentaris(self, tmp_path):
        f = tmp_path / "c.txt"
        f.write_text("# línia 1\n# línia 2\n", encoding='utf-8')
        result = read_interferent_file(str(f))
        assert result == []


# ──────────────────────────────────────────────
# 4. compute_s1_affinity
# ──────────────────────────────────────────────

class TestComputeS1Affinity:

    def test_ki_igual_al_minim_retorna_1(self):
        # ki_diana == ki_min → log_diana == log_best → s1 = 1
        s1 = compute_s1_affinity(ki_diana=1.0, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert s1 == pytest.approx(1.0)

    def test_ki_igual_al_maxim_retorna_0(self):
        s1 = compute_s1_affinity(ki_diana=100.0, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert s1 == pytest.approx(0.0)

    def test_ki_nan_retorna_0(self):
        s1 = compute_s1_affinity(ki_diana=np.nan, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert s1 == 0.0

    def test_ki_zero_retorna_0(self):
        s1 = compute_s1_affinity(ki_diana=0.0, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert s1 == 0.0

    def test_ki_min_zero_retorna_0(self):
        s1 = compute_s1_affinity(ki_diana=5.0, ki_min_matrix=0.0, ki_max_matrix=100.0)
        assert s1 == 0.0

    def test_ki_min_igual_max_retorna_0(self):
        # denominator = 0 → retorna 0
        s1 = compute_s1_affinity(ki_diana=5.0, ki_min_matrix=10.0, ki_max_matrix=10.0)
        assert s1 == 0.0

    def test_valor_intermedi_entre_0_i_1(self):
        s1 = compute_s1_affinity(ki_diana=10.0, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert 0.0 < s1 < 1.0

    def test_resultat_sempre_clipat_entre_0_i_1(self):
        # Ki molt gran (pitjor que el màxim de la matriu)
        s1 = compute_s1_affinity(ki_diana=999999.0, ki_min_matrix=1.0, ki_max_matrix=100.0)
        assert 0.0 <= s1 <= 1.0


# ──────────────────────────────────────────────
# 5. compute_s2_selectivity
# ──────────────────────────────────────────────

class TestComputeS2Selectivity:

    def test_sense_interferents_retorna_1(self):
        s2 = compute_s2_selectivity(ki_diana=5.0, ki_min_interferent=np.nan)
        assert s2 == 1.0

    def test_ki_diana_nan_retorna_1(self):
        s2 = compute_s2_selectivity(ki_diana=np.nan, ki_min_interferent=50.0)
        assert s2 == 1.0

    def test_interferent_molt_mes_gran_retorna_1(self):
        # ratio = 1000 / (1 * 10) = 100 → min(1, 100) = 1
        s2 = compute_s2_selectivity(ki_diana=1.0, ki_min_interferent=1000.0, tau=10.0)
        assert s2 == pytest.approx(1.0)

    def test_interferent_igual_diana_tau10_retorna_01(self):
        # ratio = 5 / (5 * 10) = 0.1
        s2 = compute_s2_selectivity(ki_diana=5.0, ki_min_interferent=5.0, tau=10.0)
        assert s2 == pytest.approx(0.1)

    def test_valor_clipat_a_1(self):
        s2 = compute_s2_selectivity(ki_diana=1.0, ki_min_interferent=999.0, tau=1.0)
        assert s2 <= 1.0

    def test_tau_personalitzada(self):
        # tau=1 → ratio = 10 / (5 * 1) = 2 → min(1, 2) = 1
        s2 = compute_s2_selectivity(ki_diana=5.0, ki_min_interferent=10.0, tau=1.0)
        assert s2 == pytest.approx(1.0)


# ──────────────────────────────────────────────
# 6. compute_s5_promiscuity
# ──────────────────────────────────────────────

class TestComputeS5Promiscuity:

    def test_cap_voc_unit_retorna_1(self):
        s5 = compute_s5_promiscuity(n_vocs_bound=0, n_vocs_total=50)
        assert s5 == pytest.approx(1.0)

    def test_tots_vocs_units_retorna_0(self):
        s5 = compute_s5_promiscuity(n_vocs_bound=50, n_vocs_total=50)
        assert s5 == pytest.approx(0.0)

    def test_meitat_units(self):
        s5 = compute_s5_promiscuity(n_vocs_bound=25, n_vocs_total=50)
        assert s5 == pytest.approx(0.5)

    def test_n_vocs_total_zero_retorna_1(self):
        s5 = compute_s5_promiscuity(n_vocs_bound=0, n_vocs_total=0)
        assert s5 == pytest.approx(1.0)

    def test_resultat_entre_0_i_1(self):
        s5 = compute_s5_promiscuity(n_vocs_bound=10, n_vocs_total=100)
        assert 0.0 <= s5 <= 1.0


# ──────────────────────────────────────────────
# 7. compute_final_score
# ──────────────────────────────────────────────

class TestComputeFinalScore:

    def test_score_amb_pesos_defecte(self):
        score = compute_final_score(
            s1=1.0, s2=1.0, s5=1.0,
            weights=DEFAULT_WEIGHTS
        )
        assert score == pytest.approx(1.0, abs=1e-6)

    def test_score_zeros(self):
        score = compute_final_score(
            s1=0.0, s2=0.0, s5=0.0,
            weights=DEFAULT_WEIGHTS
        )
        assert score == pytest.approx(0.0)

    def test_score_parcial(self):
        weights = {'w_affinity': 0.5, 'w_selectivity': 0.3, 'w_promiscuity': 0.2}
        score = compute_final_score(s1=1.0, s2=0.0, s5=0.0, weights=weights)
        assert score == pytest.approx(0.5)

    def test_score_retorna_float(self):
        score = compute_final_score(s1=0.8, s2=0.6, s5=0.4, weights=DEFAULT_WEIGHTS)
        assert isinstance(score, float)

    def test_score_entre_0_i_1_amb_subscores_valids(self):
        score = compute_final_score(s1=0.7, s2=0.5, s5=0.9, weights=DEFAULT_WEIGHTS)
        assert 0.0 <= score <= 1.0


# ──────────────────────────────────────────────
# 8. build_obp_ranking  (test d'integració)
# ──────────────────────────────────────────────

@pytest.fixture
def dades_exemple():
    """DataFrames mínims per provar build_obp_ranking."""
    binding = pd.DataFrame({
        'CAS':  ['64-17-5', '141-78-6', '100-52-7'],
        'Name': ['Ethanol', 'Ethyl acetate', 'Benzaldehyde'],
        'OBP_A': [5.0,  np.nan, 20.0],
        'OBP_B': [10.0, 8.0,   np.nan],
        'OBP_C': [2.0,  15.0,  3.0],
    })

    obp_info = pd.DataFrame({
        'Binding Protein Name': ['OBP_A', 'OBP_B', 'OBP_C'],
        'Binding Protein Type': ['Classic OBP', 'Other OBP', 'Classic OBP'],
        'Cystine count':        [6, 4, 6],
        'Species':              ['Apis mellifera', 'Bombyx mori', 'Apis mellifera'],
        'UniProtID':            ['P12345', '-', 'P67890'],
        'Alphafold':            ['-', '-', '-'],
    })

    ki_diana = pd.Series({'OBP_A': 5.0, 'OBP_B': 10.0, 'OBP_C': 2.0})
    weights  = {'w_affinity': 0.50, 'w_selectivity': 0.33, 'w_promiscuity': 0.17}

    return binding, obp_info, ki_diana, weights


class TestBuildObpRanking:

    def test_retorna_dataframe(self, dades_exemple):
        binding, obp_info, ki_diana, weights = dades_exemple
        result = build_obp_ranking(
            ki_values_diana=ki_diana,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        assert isinstance(result, pd.DataFrame)

    def test_columnes_esperades(self, dades_exemple):
        binding, obp_info, ki_diana, weights = dades_exemple
        result = build_obp_ranking(
            ki_values_diana=ki_diana,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        for col in ('OBP', 'Ki_diana_uM', 's1_affinity', 's2_selectivity',
                    's5_promiscuity', 'Score', 'Preferred'):
            assert col in result.columns, f"Falta la columna: {col}"

    def test_ordenat_per_score_descendent(self, dades_exemple):
        binding, obp_info, ki_diana, weights = dades_exemple
        result = build_obp_ranking(
            ki_values_diana=ki_diana,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        scores = result['Score'].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_preferred_classic_obp_6cys(self, dades_exemple):
        binding, obp_info, ki_diana, weights = dades_exemple
        result = build_obp_ranking(
            ki_values_diana=ki_diana,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        # OBP_B té Cystine count=4 → preferred=False
        obp_b_row = result[result['OBP'] == 'OBP_B']
        assert not obp_b_row.empty
        assert obp_b_row.iloc[0]['Preferred'] == False

        # OBP_A i OBP_C: Classic OBP + 6 Cys → preferred=True
        for obp in ('OBP_A', 'OBP_C'):
            row = result[result['OBP'] == obp]
            assert row.iloc[0]['Preferred'] == True

    def test_sense_dades_retorna_buit(self, dades_exemple):
        binding, obp_info, _, weights = dades_exemple
        # Ki diana tot NaN → cap OBP ha de sortir
        ki_diana_nan = pd.Series({'OBP_A': np.nan, 'OBP_B': np.nan, 'OBP_C': np.nan})
        result = build_obp_ranking(
            ki_values_diana=ki_diana_nan,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        assert result.empty

    def test_scores_entre_0_i_1(self, dades_exemple):
        binding, obp_info, ki_diana, weights = dades_exemple
        result = build_obp_ranking(
            ki_values_diana=ki_diana,
            obp_info_table=obp_info,
            binding_table=binding,
            name_col='Name',
            interferent_list=[],
            obp_name_list=['OBP_A', 'OBP_B', 'OBP_C'],
            weights=weights,
            ki_min_matrix=2.0,
            ki_max_matrix=20.0,
        )
        assert (result['Score'] >= 0.0).all()
        assert (result['Score'] <= 1.0).all()
