"""Independent truth and public-artifact oracle for molecular identity.

This module deliberately uses only the Python standard library. It derives truth from
the simulator's mutation classes and a frozen canonical X-repeat sequence; it does not
import VNtyper normalization, naming, grouping, reconciliation, or decision code.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

REFERENCE = "MUC1-X-60-coding-v1"
CANONICAL_X = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"
EXPERIMENTS = ("experiment1_dupC", "experiment2_atypical")
IDENTITY_COLUMNS = (
    "Molecular_Identity",
    "Molecular_Identity_Status",
    "Equivalent_Representation_Count",
    "Identity_Hypothesis_Count",
)
SELECTED_PROJECTION_COLUMNS = (
    "Motifs",
    "POS",
    "REF",
    "ALT",
    "Sample",
    "Motif_sequence",
    "Variant",
    "Del",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "ref_len",
    "alt_len",
    "Frame_Score",
    "is_frameshift",
    "direction",
    "frameshift_amount",
    "is_valid_frameshift",
    "Depth_Score",
    "Confidence",
    "depth_confidence_pass",
    "haplo_count",
    "alt_filter_pass",
    "motif_filter_pass",
    "Motif_fasta",
    "POS_fasta",
    "Motif",
    "Flag",
)


@dataclass(frozen=True)
class ExpectedIdentity:
    """Truth-derived stable identity and its independent display projection."""

    identity: str
    name: str


@dataclass(frozen=True)
class SampleExpectation:
    """One mutated sample's truth class and expected molecular identity."""

    mutation: str
    expected: ExpectedIdentity


@dataclass(frozen=True)
class DisplayCounts:
    """Displayed, truth-exact, and wrong public names."""

    displayed: int
    exact: int
    wrong: int


@dataclass(frozen=True)
class IdentityCounts:
    """Resolved identity accuracy on rows whose caller-local name matches truth."""

    resolved: int
    exact: int
    wrong: int


@dataclass(frozen=True)
class GoldenCorpus:
    """Independently loaded truth, historical projection, and public-row contract."""

    sim_root: Path
    advntr_root: Path
    expected_by_class: dict[str, ExpectedIdentity]
    expected_by_sample: dict[str, SampleExpectation]
    mutated_keys: frozenset[str]
    control_keys: frozenset[str]
    mutated_samples: int
    control_samples: int
    total: DisplayCounts
    by_tier: dict[str, DisplayCounts]
    tier_keys: dict[str, frozenset[str]]
    control_findings: int
    dupc_vcf_names: tuple[str, ...]
    identity_on_truth_exact_names: IdentityCounts
    identity_on_truth_exact_names_by_tier: dict[str, IdentityCounts]
    identity_outcome_keys: dict[str, frozenset[str]]
    identity_contract_violations: tuple[str, ...]
    identity_projection_by_sample: dict[str, tuple[tuple[str, ...], ...]]
    selected_projection_by_sample: dict[str, tuple[tuple[str, ...], ...]]
    policy_stable_projection_by_sample: dict[str, tuple[tuple[tuple[str, str], ...], ...]]
    unaffected_public_projection_by_sample: dict[str, tuple[tuple[tuple[str, str], ...], ...]]
    policy_explanation_by_sample: dict[str, tuple[tuple[str, str, str, str], ...]]
    complete_bam_evidence_keys: frozenset[str]
    bam_truth_match_keys: frozenset[str]
    public_truth_identity_keys: frozenset[str]
    recurrent_state_collisions: tuple[tuple[str, int, str], ...]


@dataclass(frozen=True)
class _Mutation:
    kind: str
    start: int
    end: int
    inserted: str = ""


# These are transcribed from the simulator's frozen mutation definitions, not from
# VNtyper's nomenclature configuration. Coordinates are one-based; an insertion's
# ``start`` is the base immediately to its right.
_MUTATIONS = {
    "dupC": _Mutation("insert", 60, 61, "C"),
    "insCCCC": _Mutation("insert", 60, 61, "CCCC"),
    "insG": _Mutation("insert", 59, 60, "G"),
    "dupA": _Mutation("insert", 60, 61, "A"),
    "delinsAT": _Mutation("delete_insert", 54, 56, "AT"),
    "delGCCCA": _Mutation("delete", 1, 5),
    "ins25bp": _Mutation("insert", 31, 32, "CAGGCCGGCCCCGGGCTCCGGACAC"),
    "insC_pos23": _Mutation("insert", 23, 24, "C"),
    "insG_pos58": _Mutation("insert", 58, 59, "G"),
    "insG_pos54": _Mutation("insert", 54, 55, "G"),
    "insA_pos54": _Mutation("insert", 54, 55, "A"),
}

# Independent literal transcription of the governed #267 recurrent-State evidence.
# Never import the production artifact here: this oracle must detect its drift.
RECURRENT_STATE_EVIDENCE = frozenset(
    {
        "I10_2_A_LEN1",
        "D8_2&D9_2&I9_2_A_LEN9",
        "D2_2&I2_2_C_LEN5",
        "I39_2_A_LEN4",
        "I52_2_A_LEN7",
        "I45_2_A_LEN4",
        "D45_2&I45_2_A_LEN2",
        "D14_2&I14_2_G_LEN14",
        "D58_2&D59_2",
        "I60_2_A_LEN10",
        "I14_2_G_LEN16",
        "I18_2_T_LEN1",
        "I21_2_G_LEN4",
        "D29_2&I29_2_A_LEN2",
        "D8_2&I8_2_A_LEN20",
        "D20_2&D21_2",
        "D21_2&D22_2",
        "I14_2_A_LEN1",
        "I11_2_G_LEN1",
        "I26_7_A_LEN25",
        "D17_2&D18_2&D19_2&D20_2&D21_2",
        "I14_2_C_LEN4",
        "I23_6_G_LEN1",
        "I21_2_T_LEN1",
    }
)
GOVERNED_ASSERTION = "A carried-forward recurrent adVNTR State is insufficient for molecular identity."
_POLICY_EXPLANATION_COLUMNS = frozenset({"Nomenclature_Flags", "Nomenclature_Note", "Evidence_Disposition"})
_NONPUBLIC_PERSISTENCE_COLUMNS = frozenset({"__Reconciled_Molecular_Identity"})


def assert_independent_import_closure(entrypoint: Path, repository_root: Path) -> tuple[Path, ...]:
    """Recursively AST-scan local imports and reject every production import.

    Args:
        entrypoint: Python file at the root of the closure.
        repository_root: Directory within which imports count as local.

    Returns:
        tuple[Path, ...]: Every recursively scanned local source file.

    Raises:
        AssertionError: If a file is absent or any closure member imports ``vntyper``.
    """
    root = repository_root.resolve()
    pending = [entrypoint.resolve()]
    scanned: set[Path] = set()
    forbidden: list[str] = []

    while pending:
        source = pending.pop()
        if source in scanned:
            continue
        if not source.is_file():
            raise AssertionError(f"oracle import-closure source does not exist: {source}")
        try:
            source.relative_to(root)
        except ValueError as error:
            raise AssertionError(f"oracle import-closure source escapes repository root: {source}") from error
        scanned.add(source)
        tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
        importlib_aliases = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
            if alias.name == "importlib"
        }
        import_module_aliases = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module == "importlib"
            for alias in node.names
            if alias.name == "import_module"
        }

        for node in ast.walk(tree):
            modules: list[tuple[str, int]] = []
            if isinstance(node, ast.Import):
                modules.extend((alias.name, 0) for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                modules.append((node.module or "", node.level))
                modules.extend(
                    (f"{node.module}.{alias.name}" if node.module else alias.name, node.level)
                    for alias in node.names
                    if alias.name != "*"
                )
            elif isinstance(node, ast.Call):
                dynamic_module = _literal_dynamic_import(node, importlib_aliases, import_module_aliases)
                if dynamic_module is not None:
                    modules.append((dynamic_module, 0))
            else:
                continue

            for module, level in modules:
                if module == "vntyper" or module.startswith("vntyper."):
                    forbidden.append(f"{source}: imports {module}")
                    continue
                pending.extend(
                    path for path in _resolve_local_imports(source, root, module, level) if path not in scanned
                )

    if forbidden:
        raise AssertionError("forbidden production import in oracle closure: " + "; ".join(sorted(forbidden)))
    return tuple(sorted(scanned))


def _literal_dynamic_import(
    node: ast.Call,
    importlib_aliases: set[str],
    import_module_aliases: set[str],
) -> str | None:
    """Return the literal target of a supported dynamic-import spelling."""
    if not node.args or not isinstance(node.args[0], ast.Constant) or not isinstance(node.args[0].value, str):
        return None
    if isinstance(node.func, ast.Name) and node.func.id == "__import__":
        return node.args[0].value
    if isinstance(node.func, ast.Name) and node.func.id in import_module_aliases:
        return node.args[0].value
    if (
        isinstance(node.func, ast.Attribute)
        and node.func.attr == "import_module"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id in importlib_aliases
    ):
        return node.args[0].value
    return None


def load_golden_corpus(sim_root: Path, advntr_root: Path) -> GoldenCorpus:
    """Load both explicit roots and derive truth and historical public outcomes.

    Args:
        sim_root: Root containing simulator ground truth and Kestrel artifacts.
        advntr_root: Root containing the paired final caller artifacts.

    Returns:
        GoldenCorpus: Independently derived truth and historical observations.

    Raises:
        AssertionError: If either root or any required corpus artifact is missing.
    """
    sim = _require_directory(sim_root, "simulation")
    advntr = _require_directory(advntr_root, "adVNTR")
    expected_by_class = {mutation: _expected_identity(spec) for mutation, spec in _MUTATIONS.items()}
    expected_by_sample: dict[str, SampleExpectation] = {}
    normal_keys: set[str] = set()

    for experiment in EXPERIMENTS:
        for row in _read_csv(_require_file(sim / experiment / "ground_truth.csv")):
            key = f"{experiment}/{row['pair_id']}"
            if row["condition"] == "mutated":
                mutation = row["mutation"]
                if mutation not in expected_by_class:
                    raise AssertionError(f"unrecognized simulator mutation for {key}: {mutation}")
                expected_by_sample[key] = SampleExpectation(mutation, expected_by_class[mutation])
            elif row["condition"] == "normal":
                normal_keys.add(key)

    if len(expected_by_sample) != 200 or len(normal_keys) != 200:
        raise AssertionError(
            f"incomplete truth corpus: {len(expected_by_sample)} mutated and {len(normal_keys)} controls"
        )

    tier_keys: dict[str, set[str]] = {"A": set(), "B": set(), "C": set()}
    tier_counts = {tier: [0, 0, 0] for tier in tier_keys}
    total = [0, 0, 0]
    control_findings = 0
    dupc_vcf_names: list[str] = []
    identity_counts = [0, 0, 0]
    identity_counts_by_tier = {tier: [0, 0, 0] for tier in tier_keys}
    identity_outcome_keys: dict[str, set[str]] = {
        outcome: set() for outcome in ("agreement", "disagreement", "unresolved")
    }
    violations: list[str] = []
    identity_projection: dict[str, tuple[tuple[str, ...], ...]] = {}
    selected_projection: dict[str, tuple[tuple[str, ...], ...]] = {}
    policy_stable_projection: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
    unaffected_public_projection: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
    policy_explanation: dict[str, tuple[tuple[str, str, str, str], ...]] = {}
    complete_bam_evidence_keys: set[str] = set()
    bam_truth_match_keys: set[str] = set()
    public_truth_identity_keys: set[str] = set()
    recurrent_state_collisions: list[tuple[str, int, str]] = []

    for key, expectation in sorted(expected_by_sample.items()):
        experiment, pair_id = key.split("/", maxsplit=1)
        public_rows = _public_rows(advntr, experiment, pair_id, "mutated")
        recurrent_state_collisions.extend(_recurrent_state_collisions(key, "mutated", public_rows))
        _record_identity_projection(identity_projection, key, "mutated", public_rows)
        _record_policy_projections(
            policy_stable_projection,
            unaffected_public_projection,
            policy_explanation,
            key,
            "mutated",
            public_rows,
        )
        row_violations, row_identity_counts = _identity_observations(
            key,
            "mutated",
            public_rows,
            expectation.expected,
        )
        violations.extend(row_violations)
        identity_counts[0] += row_identity_counts.resolved
        identity_counts[1] += row_identity_counts.exact
        identity_counts[2] += row_identity_counts.wrong
        for tier, counts in _identity_counts_by_tier(public_rows, expectation.expected).items():
            identity_counts_by_tier[tier][0] += counts.resolved
            identity_counts_by_tier[tier][1] += counts.exact
            identity_counts_by_tier[tier][2] += counts.wrong
        identity_outcome_keys[_identity_outcome(public_rows)].add(key)
        if _public_has_identity(public_rows, expectation.expected.identity):
            public_truth_identity_keys.add(key)
        complete_bam, bam_truth_match = _bam_replay_observation(
            advntr / experiment / pair_id / "mutated" / "kestrel",
            expectation.expected.identity,
        )
        if complete_bam:
            complete_bam_evidence_keys.add(key)
        if bam_truth_match:
            bam_truth_match_keys.add(key)
        verdicts = _displayed_verdicts(public_rows)
        if len(verdicts) > 1:
            raise AssertionError(f"conflicting public verdicts for {key}: {sorted(verdicts)}")
        if verdicts:
            name, tier = next(iter(verdicts))
            if tier not in tier_keys:
                raise AssertionError(f"unknown historical tier for {key}: {tier}")
            exact = name == expectation.expected.name
            total[0] += 1
            total[1] += int(exact)
            total[2] += int(not exact)
            tier_counts[tier][0] += 1
            tier_counts[tier][1] += int(exact)
            tier_counts[tier][2] += int(not exact)
            tier_keys[tier].add(key)

        _simulation_kestrel_rows(sim, experiment, pair_id, "mutated")
        selected_projection[f"{key}/mutated"] = _selected_projection(public_rows["kestrel"])
        if expectation.mutation == "dupC":
            for row in public_rows["kestrel"]:
                name = (row.get("Nomenclature_Kestrel") or "").strip()
                if name and name[0].isdigit():
                    dupc_vcf_names.append(name)

    for key in sorted(normal_keys):
        experiment, pair_id = key.split("/", maxsplit=1)
        public_rows = _public_rows(advntr, experiment, pair_id, "normal")
        recurrent_state_collisions.extend(_recurrent_state_collisions(key, "normal", public_rows))
        _record_identity_projection(identity_projection, key, "normal", public_rows)
        _record_policy_projections(
            policy_stable_projection,
            unaffected_public_projection,
            policy_explanation,
            key,
            "normal",
            public_rows,
        )
        row_violations, _ = _identity_observations(key, "normal", public_rows, None)
        violations.extend(row_violations)
        control_findings += int(bool(_displayed_verdicts(public_rows)))
        _simulation_kestrel_rows(sim, experiment, pair_id, "normal")
        selected_projection[f"{key}/normal"] = _selected_projection(public_rows["kestrel"])

    return GoldenCorpus(
        sim_root=sim,
        advntr_root=advntr,
        expected_by_class=expected_by_class,
        expected_by_sample=expected_by_sample,
        mutated_keys=frozenset(expected_by_sample),
        control_keys=frozenset(normal_keys),
        mutated_samples=len(expected_by_sample),
        control_samples=len(normal_keys),
        total=DisplayCounts(*total),
        by_tier={tier: DisplayCounts(*counts) for tier, counts in tier_counts.items()},
        tier_keys={tier: frozenset(keys) for tier, keys in tier_keys.items()},
        control_findings=control_findings,
        dupc_vcf_names=tuple(dupc_vcf_names),
        identity_on_truth_exact_names=IdentityCounts(*identity_counts),
        identity_on_truth_exact_names_by_tier={
            tier: IdentityCounts(*counts) for tier, counts in identity_counts_by_tier.items()
        },
        identity_outcome_keys={outcome: frozenset(keys) for outcome, keys in identity_outcome_keys.items()},
        identity_contract_violations=tuple(violations),
        identity_projection_by_sample=identity_projection,
        selected_projection_by_sample=selected_projection,
        policy_stable_projection_by_sample=policy_stable_projection,
        unaffected_public_projection_by_sample=unaffected_public_projection,
        policy_explanation_by_sample=policy_explanation,
        complete_bam_evidence_keys=frozenset(complete_bam_evidence_keys),
        bam_truth_match_keys=frozenset(bam_truth_match_keys),
        public_truth_identity_keys=frozenset(public_truth_identity_keys),
        recurrent_state_collisions=tuple(recurrent_state_collisions),
    )


def _recurrent_state_collisions(
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
) -> tuple[tuple[str, int, str], ...]:
    """Return every row matching the independent recurrent-State evidence tuple."""
    decision_key = f"{key}/{condition}"
    return tuple(
        (decision_key, ordinal, state)
        for ordinal, row in enumerate(public_rows["advntr"])
        if (state := (row.get("Variant") or "").strip()) in RECURRENT_STATE_EVIDENCE
    )


def selected_projection_fingerprint(projection: Mapping[str, tuple[tuple[str, ...], ...]]) -> str:
    """Return a deterministic fingerprint of every literal selected-row cell."""
    payload = json.dumps(projection, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def public_projection_fingerprint(projection: Mapping[str, object]) -> str:
    """Return a deterministic fingerprint for a nested public-row projection."""
    payload = json.dumps(projection, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def sample_sets_fingerprint(sample_sets: dict[str, frozenset[str]]) -> str:
    """Return a deterministic fingerprint of named literal sample-key sets."""
    payload = json.dumps(
        {name: sorted(keys) for name, keys in sample_sets.items()},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _record_identity_projection(
    projection: dict[str, tuple[tuple[str, ...], ...]],
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
) -> None:
    for caller, rows in public_rows.items():
        projection[f"{key}/{condition}/{caller}"] = tuple(
            tuple(row.get(column, "<missing>") for column in IDENTITY_COLUMNS) for row in rows
        )


def _record_policy_projections(
    stable: dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    unaffected: dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    explanations: dict[str, tuple[tuple[str, str, str, str], ...]],
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
) -> None:
    """Record B1-stable cells and independently scoped B2 explanation cells."""
    decision_key = f"{key}/{condition}"
    collision = any((row.get("Variant") or "").strip() in RECURRENT_STATE_EVIDENCE for row in public_rows["advntr"])
    for caller, rows in public_rows.items():
        projection_key = f"{decision_key}/{caller}"
        stable[projection_key] = tuple(
            tuple(
                sorted(
                    (column, value)
                    for column, value in row.items()
                    if column not in _POLICY_EXPLANATION_COLUMNS | _NONPUBLIC_PERSISTENCE_COLUMNS
                )
            )
            for row in rows
        )
        unaffected[projection_key] = tuple(
            tuple(
                sorted(
                    (column, value)
                    for column, value in row.items()
                    if column != "Evidence_Disposition" and column not in _NONPUBLIC_PERSISTENCE_COLUMNS
                )
            )
            for row in rows
            if not (
                (caller == "advntr" and (row.get("Variant") or "").strip() in RECURRENT_STATE_EVIDENCE)
                or (caller == "kestrel" and collision)
            )
        )
        explanations[projection_key] = tuple(
            (
                (row.get("Variant") or "").strip(),
                (row.get("Nomenclature_Flags") or "").strip(),
                (row.get("Nomenclature_Note") or "").strip(),
                row.get("Evidence_Disposition", "<missing>").strip(),
            )
            for row in rows
        )


def _resolve_local_imports(source: Path, root: Path, module: str, level: int) -> tuple[Path, ...]:
    if level:
        base = source.parent
        for _ in range(level - 1):
            base = base.parent
        parts = module.split(".") if module else []
        candidates = [(base, parts)]
    else:
        parts = module.split(".") if module else []
        candidates = [(root, parts), (source.parent, parts)]

    for base, candidate_parts in candidates:
        candidate = base.joinpath(*candidate_parts)
        file_candidate = candidate.with_suffix(".py")
        package_candidate = candidate / "__init__.py"
        if file_candidate.is_file():
            return _package_initializers(base, candidate_parts[:-1]) + (file_candidate.resolve(),)
        if package_candidate.is_file():
            return _package_initializers(base, candidate_parts)
    return ()


def _package_initializers(base: Path, parts: list[str]) -> tuple[Path, ...]:
    initializers: list[Path] = []
    for index in range(1, len(parts) + 1):
        initializer = base.joinpath(*parts[:index], "__init__.py")
        if initializer.is_file():
            initializers.append(initializer.resolve())
    return tuple(initializers)


def _expected_identity(spec: _Mutation) -> ExpectedIdentity:
    if spec.kind == "insert":
        start, inserted = _normalize_insertion(spec.start, spec.inserted)
        end = start - 1
        deleted = ""
    elif spec.kind in {"delete", "delete_insert"}:
        start, end = spec.start, spec.end
        deleted = CANONICAL_X[start - 1 : end]
        inserted = spec.inserted
    else:
        raise AssertionError(f"unsupported simulator mutation kind: {spec.kind}")

    if not deleted and not inserted:
        raise AssertionError("truth mutation cannot normalize to an empty edit")
    identity = f"{REFERENCE}|{start}|{end}|{deleted or '-'}|{inserted or '-'}"
    return ExpectedIdentity(identity, _render_name(start, end, deleted, inserted))


def _normalize_insertion(start: int, inserted: str) -> tuple[int, str]:
    if not inserted or set(inserted) - set("ACGT"):
        raise AssertionError(f"invalid truth insertion allele: {inserted!r}")
    boundary = start
    rotated = inserted
    # The closed identity reference has insertion starts 1..60. Do not roll an
    # end-of-unit duplication across the repeat junction to the invalid start 61.
    while boundary < len(CANONICAL_X) and CANONICAL_X[boundary - 1] == rotated[0]:
        rotated = rotated[1:] + rotated[0]
        boundary += 1
    return boundary, rotated


def _render_name(start: int, end: int, deleted: str, inserted: str) -> str:
    if deleted and inserted:
        return f"{_span(start, end)}delins{inserted}"
    if deleted:
        return f"{_span(start, end)}del{deleted}"

    preceding_end = start - 1
    preceding_start = preceding_end - len(inserted) + 1
    preceding = CANONICAL_X[preceding_start - 1 : preceding_end]
    if preceding_start >= 1 and preceding == inserted:
        return f"{_span(preceding_start, preceding_end)}dup{inserted}"

    following_end = start + len(inserted) - 1
    following = CANONICAL_X[start - 1 : following_end]
    if following_end <= len(CANONICAL_X) and following == inserted:
        return f"{_span(start, following_end)}dup{inserted}"
    return f"{end}_{start}ins{inserted}"


def _span(start: int, end: int) -> str:
    return str(start) if start == end else f"{start}_{end}"


def _parse_canonical_identity(value: object) -> tuple[tuple[int, int, str, str], ...]:
    if not isinstance(value, str):
        raise AssertionError("molecular identity must be a string")
    reference, separator, encoded_edits = value.partition("|")
    if reference != REFERENCE or not separator or not encoded_edits:
        raise AssertionError("molecular identity requires the exact versioned reference")
    edits: list[tuple[int, int, str, str]] = []
    for encoded_edit in encoded_edits.split(";"):
        fields = encoded_edit.split("|")
        if len(fields) != 4 or any(not field for field in fields):
            raise AssertionError("molecular identity edit requires four nonempty fields")
        start_text, end_text, deleted_text, inserted_text = fields
        start = _canonical_coordinate(start_text)
        end = _canonical_coordinate(end_text, allow_zero=True)
        deleted = _canonical_allele(deleted_text)
        inserted = _canonical_allele(inserted_text)
        if not 1 <= start <= len(CANONICAL_X):
            raise AssertionError("molecular identity start is outside the canonical repeat")
        if not deleted:
            if end != start - 1 or not inserted:
                raise AssertionError("molecular identity insertion coordinates are inconsistent")
        else:
            if not start <= end <= len(CANONICAL_X) or len(deleted) != end - start + 1:
                raise AssertionError("molecular identity deletion interval is inconsistent")
            if CANONICAL_X[start - 1 : end] != deleted:
                raise AssertionError("molecular identity deletion does not match the canonical repeat")
        if not deleted and not inserted:
            raise AssertionError("molecular identity edit cannot be empty")
        canonical = _canonical_edit(start, end, deleted, inserted)
        if canonical != (start, end, deleted, inserted):
            raise AssertionError("molecular identity edit is not minimal and 3-prime-most")
        edits.append(canonical)
    if edits != sorted(edits, key=lambda edit: (edit[0], edit[1])):
        raise AssertionError("molecular identity edits are not in canonical order")
    for previous, current in zip(edits, edits[1:], strict=False):
        if _identity_edits_collide(previous, current):
            raise AssertionError("molecular identity edits overlap")
    serialized = ";".join(
        f"{start}|{end}|{deleted or '-'}|{inserted or '-'}" for start, end, deleted, inserted in edits
    )
    if value != f"{REFERENCE}|{serialized}":
        raise AssertionError("molecular identity is not canonically serialized")
    return tuple(edits)


def _canonical_coordinate(value: str, *, allow_zero: bool = False) -> int:
    if not value.isascii() or not value.isdecimal() or (len(value) > 1 and value.startswith("0")):
        raise AssertionError("molecular identity coordinate is not a canonical decimal")
    coordinate = int(value)
    if coordinate < int(not allow_zero):
        raise AssertionError("molecular identity coordinate is outside the canonical repeat")
    return coordinate


def _canonical_allele(value: str) -> str:
    if value == "-":
        return ""
    if not value or any(base not in "ACGT" for base in value):
        raise AssertionError("molecular identity allele is not canonical uppercase DNA")
    return value


def _canonical_edit(start: int, end: int, deleted: str, inserted: str) -> tuple[int, int, str, str]:
    if not deleted:
        normalized_start, normalized_inserted = _normalize_insertion(start, inserted)
        return normalized_start, normalized_start - 1, "", normalized_inserted
    alternate = CANONICAL_X[: start - 1] + inserted + CANONICAL_X[end:]
    prefix = 0
    while prefix < len(CANONICAL_X) and prefix < len(alternate) and CANONICAL_X[prefix] == alternate[prefix]:
        prefix += 1
    reference_end = len(CANONICAL_X)
    alternate_end = len(alternate)
    while (
        reference_end > prefix
        and alternate_end > prefix
        and CANONICAL_X[reference_end - 1] == alternate[alternate_end - 1]
    ):
        reference_end -= 1
        alternate_end -= 1
    canonical_deleted = CANONICAL_X[prefix:reference_end]
    canonical_inserted = alternate[prefix:alternate_end]
    if not canonical_deleted and not canonical_inserted:
        raise AssertionError("molecular identity edit does not change the canonical repeat")
    canonical_start = prefix + 1
    canonical_end = reference_end if canonical_deleted else canonical_start - 1
    return canonical_start, canonical_end, canonical_deleted, canonical_inserted


def _identity_edits_collide(
    previous: tuple[int, int, str, str],
    current: tuple[int, int, str, str],
) -> bool:
    previous_start, previous_end, previous_deleted, _ = previous
    current_start, _, current_deleted, _ = current
    if not previous_deleted and not current_deleted:
        return previous_start == current_start
    if not previous_deleted:
        return current_start <= previous_start
    if not current_deleted:
        return current_start <= previous_end + 1
    return current_start <= previous_end


def _require_directory(path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_dir():
        raise AssertionError(f"{label} corpus root not found: {resolved}")
    return resolved


def _require_file(path: Path) -> Path:
    if not path.is_file():
        raise AssertionError(f"required golden artifact not found: {path}")
    return path


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_tsv(path: Path, *, comments: bool) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with _require_file(path).open(encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not comments or not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    return tuple(reader.fieldnames or ()), list(reader)


def _positive_rows(path: Path, *, comments: bool) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    fields, rows = _read_tsv(path, comments=comments)
    return fields, [row for row in rows if not _is_negative(row)]


def _is_negative(row: dict[str, str]) -> bool:
    if row.get("Confidence", "") == "Negative" or row.get("VID", "") == "Negative":
        return True
    return row.get("Motif", "") == "None" and "Motifs" not in row


def _public_rows(
    root: Path,
    experiment: str,
    pair_id: str,
    condition: str,
) -> dict[str, tuple[dict[str, str], ...]]:
    sample = root / experiment / pair_id / condition
    _, kestrel = _positive_rows(sample / "kestrel" / "kestrel_result.tsv", comments=True)
    _, advntr = _positive_rows(sample / "advntr" / "output_adVNTR_result.tsv", comments=False)
    return {"kestrel": tuple(kestrel), "advntr": tuple(advntr)}


def _simulation_kestrel_rows(root: Path, experiment: str, pair_id: str, condition: str) -> list[dict[str, str]]:
    path = root / experiment / "vntyper" / pair_id / condition / "kestrel" / "kestrel_result.tsv"
    _, rows = _positive_rows(path, comments=True)
    return rows


def _displayed_verdicts(public_rows: dict[str, tuple[dict[str, str], ...]]) -> set[tuple[str, str]]:
    verdicts: set[tuple[str, str]] = set()
    for rows in public_rows.values():
        for row in rows:
            name = (row.get("Nomenclature") or "").strip()
            tier = (row.get("Nomenclature_Tier") or "").strip()
            if name and name[0].isdigit():
                verdicts.add((name, tier))
    return verdicts


def _identity_observations(
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
    expected: ExpectedIdentity | None,
) -> tuple[list[str], IdentityCounts]:
    violations: list[str] = []
    resolved = 0
    exact_count = 0
    wrong = 0
    for caller, rows in public_rows.items():
        parsed_rows: list[tuple[str, str, str, int, int]] = []
        for index, row in enumerate(rows):
            row_key = f"{key}/{condition}/{caller}[{index}]"
            missing = tuple(column for column in IDENTITY_COLUMNS if column not in row)
            if missing:
                violations.append(f"{row_key}: missing {missing}")
                continue

            identity = (row["Molecular_Identity"] or "").strip()
            status = (row["Molecular_Identity_Status"] or "").strip()
            representations = _canonical_nonnegative_decimal(row["Equivalent_Representation_Count"])
            hypotheses = _canonical_nonnegative_decimal(row["Identity_Hypothesis_Count"])
            if representations is None or hypotheses is None:
                violations.append(f"{row_key}: identity counts are not integers")
                continue

            if status == "unresolved":
                if identity or representations != 0:
                    violations.append(f"{row_key}: inconsistent unresolved identity quartet")
            elif status in {"unique", "legacy-selected-among-multiple"}:
                if not identity or representations < 1 or hypotheses < 1:
                    violations.append(f"{row_key}: inconsistent resolved identity quartet")
                elif not _is_canonical_identity(identity):
                    violations.append(f"{row_key}: noncanonical molecular identity {identity!r}")
                elif status == "unique" and hypotheses != 1:
                    violations.append(f"{row_key}: unique status requires identity hypothesis count 1")
                elif status == "legacy-selected-among-multiple" and hypotheses <= 1:
                    violations.append(
                        f"{row_key}: legacy-selected-among-multiple status requires identity hypothesis count > 1"
                    )
                elif expected is not None:
                    caller_name_column = "Nomenclature_Kestrel" if caller == "kestrel" else "Nomenclature_adVNTR"
                    caller_name = (row.get(caller_name_column) or "").strip()
                    if caller_name == expected.name:
                        exact = identity == expected.identity
                        resolved += 1
                        exact_count += int(exact)
                        wrong += int(not exact)
                        if not exact:
                            violations.append(
                                f"{row_key}: truth-exact caller name {caller_name!r} has identity {identity!r}, "
                                f"expected {expected.identity!r}"
                            )
            else:
                violations.append(f"{row_key}: unknown identity status {status!r}")
            parsed_rows.append((row_key, identity, status, representations, hypotheses))

        hypothesis_counts = {hypotheses for _, _, _, _, hypotheses in parsed_rows}
        if len(hypothesis_counts) > 1:
            violations.append(f"{key}/{condition}/{caller}: inconsistent caller-wide hypothesis count")
        resolved_identities = {
            identity
            for _, identity, status, _, _ in parsed_rows
            if identity and status in {"unique", "legacy-selected-among-multiple"} and _is_canonical_identity(identity)
        }
        for row_key, identity, _status, representations, hypotheses in parsed_rows:
            if caller == "advntr":
                if hypotheses != len(resolved_identities):
                    violations.append(f"{row_key}: identity hypothesis count does not match emitted adVNTR identities")
                expected_representations = sum(other_identity == identity for _, other_identity, _, _, _ in parsed_rows)
                if identity and representations != expected_representations:
                    violations.append(f"{row_key}: equivalent representation count does not match emitted adVNTR rows")
            else:
                if hypotheses < len(resolved_identities):
                    violations.append(f"{row_key}: identity hypothesis count is below emitted Kestrel identities")
                observable_representations = sum(
                    other_identity == identity for _, other_identity, _, _, _ in parsed_rows
                )
                if identity and representations < observable_representations:
                    violations.append(
                        f"{row_key}: equivalent representation count is below emitted Kestrel representations"
                    )
    return violations, IdentityCounts(resolved, exact_count, wrong)


def _canonical_nonnegative_decimal(value: object) -> int | None:
    if not isinstance(value, str) or not value.isascii() or not value.isdecimal():
        return None
    if len(value) > 1 and value.startswith("0"):
        return None
    return int(value)


def _is_canonical_identity(value: str) -> bool:
    try:
        _parse_canonical_identity(value)
    except AssertionError:
        return False
    return True


def _identity_counts_by_tier(
    public_rows: dict[str, tuple[dict[str, str], ...]],
    expected: ExpectedIdentity,
) -> dict[str, IdentityCounts]:
    counts = {tier: [0, 0, 0] for tier in ("A", "B", "C")}
    for caller, rows in public_rows.items():
        caller_column = "Nomenclature_Kestrel" if caller == "kestrel" else "Nomenclature_adVNTR"
        for row in rows:
            if any(column not in row for column in IDENTITY_COLUMNS):
                continue
            identity = (row.get("Molecular_Identity") or "").strip()
            status = (row.get("Molecular_Identity_Status") or "").strip()
            caller_name = (row.get(caller_column) or "").strip()
            tier = (row.get("Nomenclature_Tier") or "").strip()
            if tier not in counts or status not in {"unique", "legacy-selected-among-multiple"}:
                continue
            if caller_name != expected.name or not identity:
                continue
            exact = identity == expected.identity
            counts[tier][0] += 1
            counts[tier][1] += int(exact)
            counts[tier][2] += int(not exact)
    return {tier: IdentityCounts(*values) for tier, values in counts.items()}


def _identity_outcome(public_rows: dict[str, tuple[dict[str, str], ...]]) -> str:
    identities: dict[str, set[str]] = {}
    for caller, rows in public_rows.items():
        identities[caller] = {
            (row.get("Molecular_Identity") or "").strip()
            for row in rows
            if (row.get("Molecular_Identity_Status") or "").strip() in {"unique", "legacy-selected-among-multiple"}
            and (row.get("Molecular_Identity") or "").strip()
        }
    if not identities["kestrel"] or not identities["advntr"]:
        return "unresolved"
    if identities["kestrel"].intersection(identities["advntr"]):
        return "agreement"
    return "disagreement"


def _public_has_identity(public_rows: dict[str, tuple[dict[str, str], ...]], identity: str) -> bool:
    return any(
        (row.get("Molecular_Identity") or "").strip() == identity for rows in public_rows.values() for row in rows
    )


def _selected_projection(rows: tuple[dict[str, str], ...]) -> tuple[tuple[str, ...], ...]:
    projection: list[tuple[str, ...]] = []
    for row in rows:
        missing = tuple(column for column in SELECTED_PROJECTION_COLUMNS if column not in row)
        if missing:
            raise AssertionError(f"selected Kestrel row missing frozen projection columns: {missing}")
        projection.append(tuple(row[column] for column in SELECTED_PROJECTION_COLUMNS))
    return tuple(projection)


def _bam_replay_observation(sample_dir: Path, expected_identity: str) -> tuple[bool, bool]:
    path = sample_dir / "bam_identity_replay.v1.json"
    if not path.is_file():
        return False, False
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_json_constant,
        )
    except (OSError, json.JSONDecodeError) as error:
        raise AssertionError(f"invalid BAM replay artifact for golden oracle: {path}") from error
    root = _exact_json_object(value, {"schema_version", "loci"}, "BAM replay root", path)
    if root["schema_version"] != "bam-identity-replay-v1":
        raise AssertionError(f"unexpected BAM replay schema for golden oracle: {path}")
    loci = root["loci"]
    if not isinstance(loci, list):
        raise AssertionError(f"BAM replay loci are not a list: {path}")
    complete = False
    truth_match = False
    previous_group: tuple[int, ...] | None = None
    retained_ordinals: set[int] = set()
    for raw_locus in loci:
        locus = _exact_json_object(
            raw_locus,
            {"candidate_observation_ordinals", "state", "evidence"},
            "BAM replay locus",
            path,
        )
        group = _strict_ordinal_group(locus["candidate_observation_ordinals"], "candidate group", path)
        if previous_group is not None and group <= previous_group:
            raise AssertionError(f"BAM replay candidate groups are not deterministically increasing: {path}")
        if retained_ordinals.intersection(group):
            raise AssertionError(f"BAM replay candidate groups overlap: {path}")
        retained_ordinals.update(group)
        previous_group = group
        state = locus["state"]
        if state not in {"not-consulted", "unavailable", "observed"}:
            raise AssertionError(f"BAM replay state is invalid: {path}")
        if state != "observed":
            if locus["evidence"] is not None:
                raise AssertionError(f"unobserved BAM replay locus has evidence: {path}")
            continue
        if locus["evidence"] is None:
            raise AssertionError(f"observed BAM replay locus has no complete evidence: {path}")
        evidence = _exact_json_object(
            locus["evidence"],
            {"eligible_record_count", "records", "counts"},
            "BAM replay evidence",
            path,
        )
        complete = True
        derived_counts = _strict_bam_evidence(evidence, group, path)
        if derived_counts:
            maximum = max(derived_counts.values())
            winners = tuple(identity for identity, count in derived_counts.items() if count == maximum)
            truth_match |= winners == (expected_identity,)
    return complete, truth_match


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise AssertionError(f"BAM replay JSON contains duplicate key: {key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise AssertionError(f"BAM replay JSON contains non-standard constant: {value}")


def _exact_json_object(value: object, keys: set[str], name: str, path: Path) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != keys:
        raise AssertionError(f"{name} must contain exactly {sorted(keys)}: {path}")
    return value


def _strict_nonnegative_integer(value: object, name: str, path: Path) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise AssertionError(f"{name} must be a non-negative integer: {path}")
    return value


def _strict_ordinal_group(value: object, name: str, path: Path) -> tuple[int, ...]:
    if not isinstance(value, list) or not value:
        raise AssertionError(f"BAM replay {name} must be a nonempty list: {path}")
    group = tuple(_strict_nonnegative_integer(ordinal, f"BAM replay {name} ordinal", path) for ordinal in value)
    if group != tuple(sorted(group)) or len(group) != len(set(group)):
        raise AssertionError(f"BAM replay {name} must be unique and increasing: {path}")
    return group


def _strict_bam_evidence(evidence: dict[str, object], group: tuple[int, ...], path: Path) -> dict[str, int]:
    records = evidence["records"]
    raw_counts = evidence["counts"]
    eligible_record_count = _strict_nonnegative_integer(
        evidence["eligible_record_count"], "BAM eligible record count", path
    )
    if not isinstance(records, list):
        raise AssertionError(f"BAM replay records are not a list: {path}")
    if eligible_record_count != len(records):
        raise AssertionError(f"BAM eligible record count does not equal retained records: {path}")
    derived_counts: dict[str, int] = {}
    for raw_record in records:
        record = _exact_json_object(
            raw_record,
            {"identities", "candidate_observation_ordinals", "minimum_kmer_depth"},
            "BAM replay record",
            path,
        )
        identities = record["identities"]
        raw_bindings = record["candidate_observation_ordinals"]
        if not isinstance(identities, list) or not isinstance(raw_bindings, list):
            raise AssertionError(f"BAM replay identities and bindings must be lists: {path}")
        if len(identities) != len(raw_bindings):
            raise AssertionError(f"BAM replay identity/binding cardinality differs: {path}")
        parsed_identities: list[str] = []
        for identity in identities:
            _parse_canonical_identity(identity)
            parsed_identities.append(identity)  # type: ignore[arg-type]
        if len(parsed_identities) != len(set(parsed_identities)):
            raise AssertionError(f"BAM replay record identities are not unique: {path}")
        bindings = tuple(_strict_ordinal_group(binding, "record identity binding", path) for binding in raw_bindings)
        bound_ordinals = tuple(ordinal for binding in bindings for ordinal in binding)
        if len(bound_ordinals) != len(set(bound_ordinals)) or not set(bound_ordinals) <= set(group):
            raise AssertionError(f"BAM replay record bindings do not belong uniquely to the locus group: {path}")
        if bindings != tuple(sorted(bindings)):
            raise AssertionError(f"BAM replay identity/binding pairs are not canonically ordered: {path}")
        minimum_depth = record["minimum_kmer_depth"]
        if minimum_depth is not None:
            depth = _strict_nonnegative_integer(minimum_depth, "BAM minimum k-mer depth", path)
            if depth > 2_147_483_647:
                raise AssertionError(f"BAM minimum k-mer depth exceeds the signed Java range: {path}")
        for identity in parsed_identities:
            derived_counts[identity] = derived_counts.get(identity, 0) + 1

    if not isinstance(raw_counts, list):
        raise AssertionError(f"BAM replay counts are not a list: {path}")
    serialized_counts: dict[str, int] = {}
    serialized_order: list[str] = []
    for raw_count in raw_counts:
        entry = _exact_json_object(raw_count, {"identity", "record_count"}, "BAM replay count", path)
        identity = entry["identity"]
        _parse_canonical_identity(identity)
        if identity in serialized_counts:
            raise AssertionError(f"BAM replay counts contain duplicate identities: {path}")
        count = _strict_nonnegative_integer(entry["record_count"], "BAM identity record count", path)
        if count < 1:
            raise AssertionError(f"BAM identity record count must be positive: {path}")
        serialized_counts[identity] = count  # type: ignore[index]
        serialized_order.append(identity)  # type: ignore[arg-type]
    if serialized_order != list(derived_counts):
        raise AssertionError(f"BAM replay counts are not in deterministic derived order: {path}")
    if serialized_counts != derived_counts:
        raise AssertionError(f"BAM replay serialized counts do not equal derived record counts: {path}")
    return derived_counts
