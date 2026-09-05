# vntyper/modules/advntr/advntr_options.py

"""
Option definitions and validation rules for adVNTR command construction.

Separated from ``advntr_genotyping.py`` to keep options handling pure,
independently testable, and strictly within repository size guidelines.
"""

from __future__ import annotations

import re
import shlex
from collections.abc import Mapping
from typing import Any

#: Every spelling adVNTR declares for an option ``run_advntr`` sets itself
#: (``advntr/__main__.py`` at the pinned revision: ``-a/--alignment_file``, ``-o/--outfile``,
#: ``-fs/--frameshift``, ``-m/--models``, ``-t/--threads``, ``-vid/--vntr_id`` and
#: ``--working_directory``).
MANAGED_ADVNTR_OPTIONS: frozenset[str] = frozenset(
    {
        "-a",
        "--alignment_file",
        "-o",
        "--outfile",
        "-fs",
        "--frameshift",
        "-m",
        "--models",
        "-t",
        "--threads",
        "-vid",
        "--vntr_id",
        "--working_directory",
    }
)

#: What owns each managed option, so a refusal can say what to set instead of only what to
#: remove.
MANAGED_ADVNTR_OPTION_OWNERS: dict[str, str] = {
    "-t": "advntr_settings['threads']",
    "--threads": "advntr_settings['threads']",
    "-vid": "advntr_settings['vid']",
    "--vntr_id": "advntr_settings['vid']",
    "-o": "the output directory and sample name run_advntr is called with",
    "--outfile": "the output directory and sample name run_advntr is called with",
    "-m": "run_advntr's db_file argument",
    "--models": "run_advntr's db_file argument",
    "-a": "run_advntr's sorted_bam argument",
    "--alignment_file": "run_advntr's sorted_bam argument",
    "--working_directory": "run_advntr's output argument",
    "-fs": "the fixed `genotype -fs` mode this module runs",
    "--frameshift": "the fixed `genotype -fs` mode this module runs",
}

#: adVNTR options added in v2.1.0/v2.2.0 that require an upgraded adVNTR binary.
ADVNTR_V22_OPTIONS: frozenset[str] = frozenset(
    {
        "--prune-reverse",
        "--exact-frameshift-caller",
        "--frameshift-background",
        "--frameshift-calibration-out",
        "--rare-unit-coverage-guard",
    }
)

#: Every other option adVNTR's ``genotype`` subparser declares -- the ones this module does
#: not set, and which ``additional_commands`` therefore exists to carry.
#:
#: ``-h``/``--help`` are deliberately absent: adVNTR would print help and exit zero without
#: genotyping anything.
ADVNTR_EXTRA_OPTIONS: frozenset[str] = frozenset(
    {
        "-r",
        "--reference_filename",
        "-f",
        "--fasta",
        "-p",
        "--pacbio",
        "-n",
        "--nanopore",
        "--outfmt",
        "--vid_file",
        "--append",
        "--noref_aln",
        "--min_read_length",
        "-e",
        "--expansion",
        "-c",
        "--coverage",
        "--haploid",
        "-naive",
        "--naive",
        "-u",
        "--update",
        "--fullru",
        "-aln",
        "--aln",
        *ADVNTR_V22_OPTIONS,
    }
)

#: What argparse treats as a negative number rather than an option string.
_NEGATIVE_NUMBER = re.compile(r"^-\d+$|^-\d*\.\d+$")


def advntr_option_token(token: str) -> str | None:
    """Return the option spelling ``token`` is, or ``None`` if it is a value.

    Only the exact spelling is returned -- ``--opt=value`` yields ``--opt``, and nothing
    else is normalised. That is the point: see :func:`resolve_additional_commands` for why
    this is an allow-list rather than a pattern match.

    Args:
        token (str): One ``shlex``-split word of ``additional_commands``.

    Returns:
        str | None: The option spelling, or ``None`` if the word is a value.
    """
    if not token.startswith("-") or _NEGATIVE_NUMBER.match(token):
        return None
    return token.split("=", 1)[0] if token.startswith("--") else token


def resolve_additional_commands(
    settings: Mapping[str, Any],
    advntr_version: tuple[int, int, int] | None = None,
) -> str:
    """Return ``additional_commands``, refusing any option this module already sets.

    Args:
        settings: An ``advntr_settings`` mapping.
        advntr_version: Optional detected adVNTR version tuple (major, minor, patch).
            When provided, validates that requested options are supported by the installed
            binary.

    Returns:
        str: The validated command fragment string.

    Raises:
        ValueError: If the fragment does not parse as shell words, if any word is an option
            :func:`run_advntr` sets itself, if any word is not an option adVNTR declares,
            if an option requires a newer adVNTR version than installed, or if
            `--exact-frameshift-caller` is requested without `--frameshift-background`.
    """
    additional = str(settings.get("additional_commands", ""))

    try:
        tokens = shlex.split(additional)
    except ValueError as exc:
        raise ValueError(
            f"advntr_settings['additional_commands'] is not a parseable command fragment ({exc}): "
            f"{additional!r}. It is interpolated into a string run under `bash -c`, so an "
            "unbalanced quote would surface as a shell syntax error inside adVNTR's log."
        ) from exc

    seen_options: set[str] = set()

    for token in tokens:
        option = advntr_option_token(token)
        if option is None:
            continue
        seen_options.add(option)
        if option in MANAGED_ADVNTR_OPTIONS:
            owner = MANAGED_ADVNTR_OPTION_OWNERS[option]
            raise ValueError(
                f"advntr_settings['additional_commands'] contains {token!r}, which is adVNTR's "
                f"{option} -- an option run_advntr already sets. adVNTR parses with argparse, where "
                f"the last occurrence wins, so this would silently override {owner}. Set that instead, "
                "and keep additional_commands for flags adVNTR alone owns (such as -aln)."
            )
        if option not in ADVNTR_EXTRA_OPTIONS:
            raise ValueError(
                f"advntr_settings['additional_commands'] contains {token!r}, which adVNTR's `genotype` "
                "does not declare as an option. Abbreviations, attached values and single-dash flag "
                "groups are refused even when argparse would accept them, because each is a way to "
                "reach an option run_advntr sets -- `-pt3` sets the thread count and `-v` reaches "
                "-vid. Spell the option out in full, one word per option."
            )
        if advntr_version is not None and advntr_version < (2, 2, 0) and option in ADVNTR_V22_OPTIONS:
            raise ValueError(
                f"advntr_settings['additional_commands'] contains {token!r}, which requires adVNTR >= 2.2.0, "
                f"but installed adVNTR version is {advntr_version[0]}.{advntr_version[1]}.{advntr_version[2]}."
            )

    if "--exact-frameshift-caller" in seen_options and "--frameshift-background" not in seen_options:
        raise ValueError(
            "advntr_settings['additional_commands'] contains '--exact-frameshift-caller' without "
            "'--frameshift-background'. adVNTR's exact frameshift caller requires an explicit frozen "
            "background model file."
        )

    return additional
