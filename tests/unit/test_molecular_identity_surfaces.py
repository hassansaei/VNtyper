"""Molecular-identity cells on positive caller and diagnostic surfaces."""

from __future__ import annotations

import json
from pathlib import Path
from typing import cast

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.modules.advntr import advntr_genotyping
from vntyper.scripts import kestrel_genotyping
from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    IDENTITY_SELECTION_COLUMNS,
    PersistedIdentityCandidate,
)
from vntyper.scripts.identity_candidates import LEGACY_GATE_COLUMNS, IdentityTranslationComponent
from vntyper.scripts.molecular_identity import (
    IdentityTranslation,
    make_coding_edit,
    make_molecular_identity,
    serialize_molecular_identity,
)
from vntyper.scripts.molecular_identity_presentation import (
    IDENTITY_COLUMNS,
    IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS,
    identity_result_cells,
    identity_translation_diagnostic_cells,
    persisted_identity_result_cells,
)
from vntyper.scripts.nomenclature import load_nomenclature_config
from vntyper.scripts.nomenclature_annotate import annotate_kestrel_frame

pytestmark = pytest.mark.unit

_DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_OTHER = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))
_RESOLVED = IdentityTranslation(_DUPC, "resolved", None, False)
_RESOLVED_OTHER = IdentityTranslation(_OTHER, "resolved", None, False)
_UNRESOLVED = IdentityTranslation(None, "unresolved", "motif-anchor-mismatch", False)

_NEGATIVE_KESTREL_COLUMNS = (
    "Motif",
    "Variant",
    "POS",
    "REF",
    "ALT",
    "Motif_sequence",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "Depth_Score",
    "Confidence",
)

_NEGATIVE_ADVNTR_COLUMNS = (
    "VID",
    "Variant",
    "NumberOfSupportingReads",
    "MeanCoverage",
    "Pvalue",
    "RU",
    "POS",
    "REF",
    "ALT",
    "Flag",
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
)


def test_identity_columns_have_the_exact_public_order_and_spelling() -> None:
    assert IDENTITY_COLUMNS == (
        "Molecular_Identity",
        "Molecular_Identity_Status",
        "Equivalent_Representation_Count",
        "Identity_Hypothesis_Count",
    )


def test_unique_positive_row_has_one_identity_and_representation() -> None:
    assert identity_result_cells(_RESOLVED) == {
        "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
        "Molecular_Identity_Status": "unique",
        "Equivalent_Representation_Count": 1,
        "Identity_Hypothesis_Count": 1,
    }


def test_selected_identity_among_multiple_hypotheses_is_explicit() -> None:
    cells = identity_result_cells(_RESOLVED, _RESOLVED_OTHER)
    assert cells["Molecular_Identity_Status"] == "legacy-selected-among-multiple"
    assert cells["Equivalent_Representation_Count"] == 1
    assert cells["Identity_Hypothesis_Count"] == 2


def test_unresolved_positive_row_has_consistent_cells() -> None:
    cells = identity_result_cells(_UNRESOLVED, _RESOLVED_OTHER)
    assert cells == {
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 1,
    }


def test_unresolved_positive_row_without_other_hypotheses_reports_zero() -> None:
    cells = identity_result_cells(_UNRESOLVED)
    assert cells["Molecular_Identity"] == ""
    assert cells["Equivalent_Representation_Count"] == 0
    assert cells["Identity_Hypothesis_Count"] == 0


def test_equivalent_positive_representations_share_one_hypothesis() -> None:
    cells = identity_result_cells(_RESOLVED, _RESOLVED)
    assert cells["Molecular_Identity_Status"] == "unique"
    assert cells["Equivalent_Representation_Count"] == 2
    assert cells["Identity_Hypothesis_Count"] == 1


def test_presentation_rejects_non_translation_inputs() -> None:
    with pytest.raises(ValueError, match="Identity presentation requires IdentityTranslation"):
        identity_result_cells(_RESOLVED, "not-a-translation")  # type: ignore[arg-type]


def test_unresolved_translation_diagnostics_are_empty_and_explicit() -> None:
    assert identity_translation_diagnostic_cells(_UNRESOLVED) == {
        "Molecular_Identity": "",
        "Molecular_Identity_Translation_Status": "unresolved",
        "Molecular_Identity_Translation_Failure": "motif-anchor-mismatch",
        "Molecular_Identity_Context_Diverges": False,
    }


def test_translation_diagnostics_reject_non_translation_input() -> None:
    with pytest.raises(ValueError, match="diagnostics require an IdentityTranslation"):
        identity_translation_diagnostic_cells("not-a-translation")  # type: ignore[arg-type]


def test_persisted_presentation_rejects_non_candidate_input() -> None:
    with pytest.raises(ValueError, match="requires PersistedIdentityCandidate"):
        persisted_identity_result_cells("not-a-candidate")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("equivalent_count", "hypothesis_count", "message"),
    [(0, 1, "positive equivalent"), (1, 0, "positive hypothesis")],
)
def test_persisted_resolved_presentation_rejects_impossible_counts(
    equivalent_count: int,
    hypothesis_count: int,
    message: str,
) -> None:
    persisted = PersistedIdentityCandidate(
        translation=_RESOLVED,
        selected_row_key=_persisted_key(),
        selected_observation_ordinal=0,
        equivalent_representation_count=equivalent_count,
        identity_hypothesis_count=hypothesis_count,
        blocking_gates=frozenset(),
        flags=frozenset(),
        group_context_diverges=False,
    )
    with pytest.raises(ValueError, match=message):
        persisted_identity_result_cells(persisted)


def test_persisted_unresolved_presentation_rejects_nonzero_equivalent_count() -> None:
    persisted = PersistedIdentityCandidate(
        translation=_UNRESOLVED,
        selected_row_key=_persisted_key(),
        selected_observation_ordinal=0,
        equivalent_representation_count=1,
        identity_hypothesis_count=0,
        blocking_gates=frozenset(),
        flags=frozenset(),
        group_context_diverges=False,
    )
    with pytest.raises(ValueError, match="unresolved identity requires zero equivalent"):
        persisted_identity_result_cells(persisted)


def test_persisted_unresolved_presentation_retains_other_hypothesis_count() -> None:
    persisted = PersistedIdentityCandidate(
        translation=_UNRESOLVED,
        selected_row_key=_persisted_key(),
        selected_observation_ordinal=0,
        equivalent_representation_count=0,
        identity_hypothesis_count=1,
        blocking_gates=frozenset(),
        flags=frozenset(),
        group_context_diverges=False,
    )
    assert persisted_identity_result_cells(persisted) == {
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 1,
    }


@pytest.mark.parametrize("invalid_count", [True, -1])
def test_persisted_presentation_rejects_noncanonical_counts(invalid_count: object) -> None:
    persisted = PersistedIdentityCandidate(
        translation=_RESOLVED,
        selected_row_key=_persisted_key(),
        selected_observation_ordinal=0,
        equivalent_representation_count=invalid_count,  # type: ignore[arg-type]
        identity_hypothesis_count=1,
        blocking_gates=frozenset(),
        flags=frozenset(),
        group_context_diverges=False,
    )
    with pytest.raises(ValueError, match="non-negative integer"):
        persisted_identity_result_cells(persisted)


def test_positive_kestrel_annotation_appends_the_public_quartet() -> None:
    frame = kestrel_stage_frame("final")
    for column, value in _persisted_kestrel_cells().items():
        frame[column] = value

    annotated = annotate_kestrel_frame(frame)

    assert tuple(annotated.columns[-4:]) == IDENTITY_COLUMNS
    assert annotated.loc[0, list(IDENTITY_COLUMNS)].to_dict() == identity_result_cells(_RESOLVED)


def test_kestrel_pre_result_has_only_the_four_approved_identity_diagnostics(tmp_path: Path) -> None:
    config = load_nomenclature_config()
    motifs = config["motifs"]
    rows = pd.DataFrame(
        [
            {
                **kestrel_stage_frame(
                    "raw", rows=1, motifs="S-C", pos=67, ref="G", alt="GG", depth_alt=7, depth_region=500
                )
                .iloc[0]
                .to_dict(),
                "Motif_sequence": motifs["C"] + motifs["S"],
            },
            {
                **kestrel_stage_frame(
                    "raw", rows=1, motifs="A-J", pos=67, ref="C", alt="CG", depth_alt=1, depth_region=500
                )
                .iloc[0]
                .to_dict(),
                "Motif_sequence": motifs["J"] + motifs["A"],
            },
        ]
    )
    merged = pd.DataFrame({"Motif": ["S", "A"], "Motif_sequence": [motifs["S"], motifs["A"]]})

    result = kestrel_genotyping.process_kmer_results(
        rows,
        merged,
        str(tmp_path),
        kestrel_config(),
        identity_component=cast(IdentityTranslationComponent, _DistinctKestrelTranslations()),
    )
    annotated = annotate_kestrel_frame(result)
    pre_result = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t", keep_default_na=False)

    assert len(pre_result) == 2
    assert pre_result["Motifs"].tolist() == ["S-C", "A-J"]
    assert set(LEGACY_GATE_COLUMNS).issubset(pre_result.columns)
    assert set(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS).issubset(pre_result.columns)
    assert pre_result["Molecular_Identity_Translation_Status"].tolist() == ["resolved", "resolved"]
    assert not {
        "Molecular_Identity_Status",
        "Equivalent_Representation_Count",
        "Identity_Hypothesis_Count",
    }.intersection(pre_result.columns)
    assert tuple(annotated.columns[-4:]) == IDENTITY_COLUMNS
    assert annotated.loc[0, "Molecular_Identity_Status"] == "unique"
    assert annotated.loc[0, "Identity_Hypothesis_Count"] == 1


def test_positive_advntr_publication_appends_all_identity_fields(tmp_path: Path) -> None:
    source = tmp_path / "output_adVNTR.vcf"
    source.write_text(
        "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n25561\tI22_2_G_LEN1\t11\t153.98\t0.0001\n",
        encoding="utf-8",
    )

    advntr_genotyping.process_advntr_output(str(source), str(tmp_path), "output")

    result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", keep_default_na=False)
    assert tuple(result.columns[-4:]) == IDENTITY_COLUMNS
    assert result.loc[0, "Molecular_Identity"] == "MUC1-X-60-coding-v1|60|59|-|C"
    assert result.loc[0, "Molecular_Identity_Status"] == "unique"
    assert result.loc[0, "Equivalent_Representation_Count"] == 1
    assert result.loc[0, "Identity_Hypothesis_Count"] == 1


def test_negative_caller_schemas_have_no_identity_fields_and_keep_one_kestrel_note(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    note = "One candidate met all gates except the reporting floor."
    kestrel_genotyping.output_empty_result(kestrel_dir, ["## VNtyper Kestrel result"], note=note)
    kestrel_text = (kestrel_dir / "kestrel_result.tsv").read_text(encoding="utf-8")
    kestrel = pd.read_csv(kestrel_dir / "kestrel_result.tsv", sep="\t", comment="#", dtype=str)

    advntr_source = tmp_path / "output_adVNTR.vcf"
    advntr_source.write_text(
        "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n",
        encoding="utf-8",
    )
    advntr_genotyping.process_advntr_output(str(advntr_source), str(tmp_path), "output")
    advntr = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)

    assert tuple(kestrel.columns) == _NEGATIVE_KESTREL_COLUMNS
    assert kestrel_text.count(f"## {note}\n") == 1
    assert "Molecular_Identity" not in kestrel_text
    assert tuple(advntr.columns) == _NEGATIVE_ADVNTR_COLUMNS
    assert not set(IDENTITY_COLUMNS).intersection(advntr.columns)


def _persisted_key():
    from vntyper.scripts.identity_candidates import RawRepresentationKey

    return RawRepresentationKey("kestrel", ("X-5", 67, "G", "GG"))


def _persisted_kestrel_cells() -> dict[str, str]:
    raw_key = json.dumps({"source": "kestrel", "values": ["X-5", 67, "G", "GG"]}, separators=(",", ":"), sort_keys=True)
    return {
        IDENTITY_CAPTURE_COLUMNS[0]: raw_key,
        IDENTITY_CAPTURE_COLUMNS[1]: serialize_molecular_identity(_DUPC),
        IDENTITY_CAPTURE_COLUMNS[2]: "resolved",
        IDENTITY_CAPTURE_COLUMNS[3]: "absent",
        IDENTITY_CAPTURE_COLUMNS[4]: "false",
        IDENTITY_CAPTURE_COLUMNS[5]: "0",
        IDENTITY_SELECTION_COLUMNS[0]: raw_key,
        IDENTITY_SELECTION_COLUMNS[1]: "1",
        IDENTITY_SELECTION_COLUMNS[2]: "1",
        IDENTITY_SELECTION_COLUMNS[3]: "[]",
        IDENTITY_SELECTION_COLUMNS[4]: "[]",
        IDENTITY_SELECTION_COLUMNS[5]: "0",
        IDENTITY_SELECTION_COLUMNS[6]: "false",
    }


class _DistinctKestrelTranslations:
    """Resolve the two test motif contexts to distinct identities."""

    def translate_kestrel(self, representation) -> IdentityTranslation:
        return _RESOLVED if representation.motifs == "S-C" else _RESOLVED_OTHER
