"""Identify the adVNTR model a run will use, and refuse the combinations that lie.

adVNTR derives its read-fetch window from the model's own content. Before #268 the
shipped hg38 model described 840 bp of an array GRCh38 carries at 3,525 bp, so adVNTR
fetched a quarter of the locus and every read outside it was invisible -- silently.
A v2 model records the array's genomic end explicitly, which is what corrects it.

That correction only holds if the model and the binary agree. adVNTR 2.0.3 selects the
legacy columns by name, so it would ignore a recorded end and reproduce the truncated
window with no error at all. Wrong answers are worse than a refusal, so this module
refuses, and it does so before the run rather than after.

It also produces the provenance record `pipeline_summary.json` carries, so a result can
be traced to the model that produced it.
"""

from __future__ import annotations

import hashlib
import re
import shlex
import sqlite3
import subprocess
from pathlib import Path
from typing import Any

#: adVNTR releases before this ignore a model's recorded genomic end.
SPAN_AWARE_ADVNTR = (2, 0, 4)

_VERSION = re.compile(r"(\d+)\.(\d+)\.(\d+)")

#: Every statement this module runs, keyed by the two table names it accepts. No SQL is
#: assembled from a variable, so there is no identifier for anything to be injected into.
#: SQLite cannot parameterise an identifier, so a fixed literal is the only safe form.
_PRAGMA = {
    "vntrs": "PRAGMA table_info(vntrs)",
    "vntrs_v2": "PRAGMA table_info(vntrs_v2)",
}
_SELECT = {
    "vntrs": {
        False: "SELECT id, chromosome, ref_start, repeats FROM vntrs",
        True: "SELECT id, chromosome, ref_start, repeats, ref_end FROM vntrs",
    },
    "vntrs_v2": {
        False: "SELECT id, chromosome, ref_start, repeats FROM vntrs_v2",
        True: "SELECT id, chromosome, ref_start, repeats, ref_end FROM vntrs_v2",
    },
}

#: The MUC1 VNTR. A database carrying more rows is selected from, not rejected.
MUC1_VID = 25561

_REINSTALL = "Reinstall the reference bundle with `vntyper install-references -d reference --references hg38`."


class AdvntrModelError(RuntimeError):
    """The resolved adVNTR model cannot be used with the available adVNTR."""


def parse_advntr_version(text: str | None) -> tuple[int, int, int] | None:
    """Extract a version tuple from `advntr --version` output.

    Returns None when no version can be read. adVNTR 2.0.3 has no `--version` flag, so
    an unreadable answer is the ordinary signal for "too old", not an anomaly -- and it
    is never guessed at.
    """
    if not text:
        return None
    match = _VERSION.search(text)
    if match is None:
        return None
    return tuple(int(part) for part in match.groups())  # type: ignore[return-value]


def detect_advntr_version(config: dict[str, Any]) -> tuple[int, int, int] | None:
    """Ask the configured adVNTR which version it is.

    Returns None when the question cannot be answered -- an absent binary, a failure, or
    a release with no `--version` flag, which is exactly what 2.0.3 is. Callers treat
    None as "too old" rather than "probably fine".
    """
    command = config.get("tools", {}).get("advntr")
    if not command:
        return None
    try:
        # Split rather than handing the string to a shell: the configured command is
        # multi-token ("mamba run -n envadvntr advntr"), which is the only reason a
        # shell was tempting, and shlex gives the same tokens without one.
        completed = subprocess.run(
            [*shlex.split(command), "--version"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, ValueError, subprocess.SubprocessError):
        return None
    return parse_advntr_version(f"{completed.stdout}\n{completed.stderr}")


def describe_model(db_path: str | Path) -> dict[str, Any]:
    """Return the provenance record for a model database.

    Deliberately reports `window_bp` -- the interval adVNTR will actually fetch -- and
    never a bare sum of segment lengths presented as a genomic span. Conflating those
    two quantities is the defect itself.
    """
    path = Path(db_path)
    if not path.is_file():
        raise AdvntrModelError(f"adVNTR model not found: {path}. {_REINSTALL}")

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    connection = sqlite3.connect(str(path))
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        # Statements are selected from a fixed map rather than built by interpolating a
        # name into SQL. SQLite cannot parameterise an identifier, so the only safe form
        # is a literal the caller never influences.
        if "vntrs_v2" in tables:
            table, schema_version = "vntrs_v2", "v2"
        elif "vntrs" in tables:
            table, schema_version = "vntrs", "v1"
        else:
            raise AdvntrModelError(
                f"{path} contains no adVNTR model table (expected `vntrs_v2` or `vntrs`). {_REINSTALL}"
            )

        columns = {row[1] for row in connection.execute(_PRAGMA[table])}
        rows = list(connection.execute(_SELECT[table]["ref_end" in columns]))
    finally:
        connection.close()

    if not rows:
        raise AdvntrModelError(f"{path} holds no VNTR rows. {_REINSTALL}")
    if len(rows) == 1:
        row = rows[0]
    else:
        # A general adVNTR database carries thousands of VNTRs. Select MUC1 rather than
        # assuming the file is the single-row bundle.
        muc1 = [candidate for candidate in rows if candidate[0] == MUC1_VID]
        if len(muc1) != 1:
            raise AdvntrModelError(
                f"{path} holds {len(rows)} VNTR rows and {len(muc1)} with VID "
                f"{MUC1_VID}; exactly one MUC1 model is required. {_REINSTALL}"
            )
        row = muc1[0]
    vid, chromosome, ref_start, repeats = row[0], row[1], int(row[2]), row[3]
    ref_end = int(row[4]) if len(row) > 4 and row[4] is not None else None
    segments = repeats.split(",")

    if ref_end is None:
        # No recorded end: adVNTR falls back to start + sum(len(units)), which is the
        # window it will really use. Report that, so the truncation is visible.
        ref_end = ref_start + sum(len(segment) for segment in segments)
    if ref_end <= ref_start:
        raise AdvntrModelError(f"{path} records ref_end={ref_end} at or before ref_start={ref_start}. {_REINSTALL}")

    return {
        "path": str(path),
        "sha256": digest,
        "schema_version": schema_version,
        "vid": vid,
        "genomic_interval": f"{chromosome}:{ref_start}-{ref_end}",
        "window_bp": ref_end - ref_start,
        "n_segments": len(segments),
        "n_distinct_segments": len(set(segments)),
        "max_segment_len": max(len(segment) for segment in segments),
    }


def require_compatible_advntr(model: dict[str, Any], advntr_version: tuple[int, int, int] | None) -> None:
    """Raise unless this model can be trusted with this adVNTR.

    Refuses two combinations, both of which would otherwise produce wrong results in
    silence: a v2 model read by an adVNTR that ignores its recorded end, and a v1
    model, which carries the truncated hg38 window regardless of the binary.
    """
    if model["schema_version"] == "v1":
        raise AdvntrModelError(
            f"{model['path']} is a v1 adVNTR model. On GRCh38 it describes "
            f"{model['window_bp']} bp of a 3,525 bp repeat array, so adVNTR would see "
            f"part of the locus and miss variants outside it "
            f"(hassansaei/VNtyper#268). {_REINSTALL}"
        )

    if advntr_version is None or advntr_version < SPAN_AWARE_ADVNTR:
        found = ".".join(str(part) for part in advntr_version) if advntr_version else "unknown"
        wanted = ".".join(str(part) for part in SPAN_AWARE_ADVNTR)
        raise AdvntrModelError(
            f"adVNTR {found} cannot read the recorded genomic end in "
            f"{model['path']}; it would fetch a truncated window and report "
            f"confident but incomplete results. Install adVNTR >= {wanted}."
        )
