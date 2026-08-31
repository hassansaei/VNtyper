"""Independent truth and public-artifact oracle for molecular identity.

This module deliberately uses only the Python standard library. It derives truth from
the simulator's mutation classes and a frozen canonical X-repeat sequence; it does not
import VNtyper normalization, naming, grouping, reconciliation, or decision code.
"""

from __future__ import annotations

import ast
import csv
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
class GoldenCorpus:
    """Independently loaded truth, historical projection, and public-row contract."""

    sim_root: Path
    advntr_root: Path
    expected_by_class: dict[str, ExpectedIdentity]
    expected_by_sample: dict[str, SampleExpectation]
    mutated_keys: frozenset[str]
    mutated_samples: int
    control_samples: int
    total: DisplayCounts
    by_tier: dict[str, DisplayCounts]
    tier_keys: dict[str, frozenset[str]]
    control_findings: int
    kestrel_positive_keys: frozenset[str]
    dupc_vcf_names: tuple[str, ...]
    identity_contract_violations: tuple[str, ...]


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
    kestrel_positive_keys: set[str] = set()
    dupc_vcf_names: list[str] = []
    violations: list[str] = []

    for key, expectation in sorted(expected_by_sample.items()):
        experiment, pair_id = key.split("/", maxsplit=1)
        public_rows = _public_rows(advntr, experiment, pair_id, "mutated")
        violations.extend(_identity_violations(key, "mutated", public_rows))
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

        raw_kestrel = _simulation_kestrel_rows(sim, experiment, pair_id, "mutated")
        if raw_kestrel:
            kestrel_positive_keys.add(key)
        if expectation.mutation == "dupC":
            for row in public_rows["kestrel"]:
                name = (row.get("Nomenclature_Kestrel") or "").strip()
                if name and name[0].isdigit():
                    dupc_vcf_names.append(name)

    for key in sorted(normal_keys):
        experiment, pair_id = key.split("/", maxsplit=1)
        public_rows = _public_rows(advntr, experiment, pair_id, "normal")
        violations.extend(_identity_violations(key, "normal", public_rows))
        control_findings += int(bool(_displayed_verdicts(public_rows)))
        _simulation_kestrel_rows(sim, experiment, pair_id, "normal")

    return GoldenCorpus(
        sim_root=sim,
        advntr_root=advntr,
        expected_by_class=expected_by_class,
        expected_by_sample=expected_by_sample,
        mutated_keys=frozenset(expected_by_sample),
        mutated_samples=len(expected_by_sample),
        control_samples=len(normal_keys),
        total=DisplayCounts(*total),
        by_tier={tier: DisplayCounts(*counts) for tier, counts in tier_counts.items()},
        tier_keys={tier: frozenset(keys) for tier, keys in tier_keys.items()},
        control_findings=control_findings,
        kestrel_positive_keys=frozenset(kestrel_positive_keys),
        dupc_vcf_names=tuple(dupc_vcf_names),
        identity_contract_violations=tuple(violations),
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


def _identity_violations(
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
) -> list[str]:
    violations: list[str] = []
    for caller, rows in public_rows.items():
        for index, row in enumerate(rows):
            row_key = f"{key}/{condition}/{caller}[{index}]"
            missing = tuple(column for column in IDENTITY_COLUMNS if column not in row)
            if missing:
                violations.append(f"{row_key}: missing {missing}")
                continue

            identity = (row["Molecular_Identity"] or "").strip()
            status = (row["Molecular_Identity_Status"] or "").strip()
            try:
                representations = int(row["Equivalent_Representation_Count"])
                hypotheses = int(row["Identity_Hypothesis_Count"])
            except (TypeError, ValueError):
                violations.append(f"{row_key}: identity counts are not integers")
                continue

            if status == "unresolved":
                if identity or representations != 0 or hypotheses < 0:
                    violations.append(f"{row_key}: inconsistent unresolved identity quartet")
            elif status in {"unique", "legacy-selected-among-multiple"}:
                if not identity or representations < 1 or hypotheses < 1:
                    violations.append(f"{row_key}: inconsistent resolved identity quartet")
            else:
                violations.append(f"{row_key}: unknown identity status {status!r}")
    return violations
