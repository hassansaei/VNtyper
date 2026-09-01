#!/usr/bin/env python3
"""
Differential sweep for #192 -- summing every ``LEN`` token in an adVNTR ``State``.

Why this script exists
----------------------
#192 changes ``Insertion_len``, which feeds
``Net_indel_length = Insertion_len - Deletion_length``, which is the value the #182
pathogenic-frame filter tests. It therefore changes **reported genotypes**, and the
golden-cohort gate cannot see it: the only compound state in the cohort is
``example_dfc3``'s ``D17_2&D18_2&D19_2&D20_2&D21_2``, which carries no ``LEN`` token at
all, so ``Insertion_len`` is 0 under both semantics and a cohort PASS is not evidence for
this change. This sweep is the evidence.

The filter it models
--------------------
A row is reported exactly when its **signed** net change satisfies
``Net_indel_length % 3 == 1``: a net insertion of ``3n+1`` bases (kept by
``advntr_processing_ins``) or a net deletion of ``3n+2`` bases (kept by
``advntr_processing_del``). ``frame`` is ``|Net_indel_length|``, the magnitude the two
configurable series are matched against; each arm tests the sign separately.

Both arms previously tested ``frame`` alone, guarded only by "names at least one
deletion" / "names at least one insertion" -- guards a **mixed** state satisfies on both
sides, so a mixed state could be reported by the arm whose series its magnitude landed
in, whatever the direction of the change actually was. :func:`absolute_frame_survival`
keeps that retired model so the sweep can report exactly which state shapes the sign
repair moves; see the ``sign_fix_delta`` section of the result.

What it compares
----------------
For every generated ``State`` string it computes, on both sides:

* ``Insertion_len`` -- historic replica vs the production ``sum_insertion_lengths``;
* ``frame`` -- ``abs(Insertion_len - Deletion_length)``, with ``Deletion_length`` the
  unchanged ``Variant.str.count("D")``;
* the keep/drop outcome of each half of the filter, and whether the state is **reported**
  at all (``process_advntr_output`` concatenates the two halves, so either one keeping the
  row is enough).

Both sides are evaluated through the **corrected** filter, so ``reporting_delta`` is
"what #192's ``Insertion_len`` change does, measured under the signed #182 rule".
The separate ``sign_fix_delta`` isolates the other axis: same (summed) ``Insertion_len``,
retired absolute-value filter versus the signed one.

The new side's model of the filter is cross-checked against the real
``advntr_processing_ins`` / ``advntr_processing_del`` before anything is reported, so the
sweep cannot pass by comparing two of its own approximations.

The oracle, stated before the sweep is run
------------------------------------------
The historic derivation extracted greedily from the first literal ``LEN`` to the end of
the string, split that remainder once on ``LEN``, and coerced what was left with
``pd.to_numeric(errors="coerce")``. So it produced the intended length exactly when the
first literal ``LEN`` was followed by nothing but digits, and **zero** otherwise. The new
derivation sums every ``LEN<n>`` token. Hence:

    a state differs iff the text following its first literal "LEN" is not a plain
    digit run, AND the sum of its LEN<n> tokens is non-zero.

In the terms the plan uses: a state differs **iff material follows its first ``LEN``
token**. States with no ``LEN``, a single terminal ``LEN``, or a ``LEN`` only in the last
part must be byte-identical -- the generator includes all three classes deliberately,
because a sweep that omitted them would prove nothing about them. A difference in any of
those three classes is a regression and exits non-zero.

Two refinements the shorthand does not cover, reachable only by strings adVNTR cannot
emit and reported separately as an annex:

* ``LEN0`` -- material follows the first ``LEN`` yet both sides give 0, so the state does
  not differ (e.g. ``I9_2_A_LEN0&D50_2``);
* a stray ``LEN`` not followed by digits -- ``I9_2_A_LENX&I50_2_A_LEN3`` has its only real
  token in the last part, yet it differs, because the historic expression anchored on the
  malformed literal. The new value (3) is the correct one.

Usage:
    python scripts/advntr_len_differential.py [--out /tmp/advntr-len-diff.json]

Exit codes:
    0: every difference is predicted by the oracle, the three unchanged classes are
       byte-identical, and the sign repair moves no pure state.
    1: a regression -- a difference the oracle does not predict, a predicted difference
       that did not occur, a pure state whose verdict the sign repair changes, or a
       mismatch between the sweep's filter model and production.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.modules.advntr.advntr_decision_config import project_advntr_settings
from vntyper.scripts.run_configuration import resolve_run_configuration

#: A well-formed ``LEN<n>`` token, as the production parser defines one.
LEN_TOKEN = re.compile(r"LEN(\d+)")

#: Any literal ``LEN``, well-formed or not. Used only to detect the stray-literal class.
LEN_LITERAL = re.compile(r"LEN")

#: The part alphabet for the systematic product. Every member is a shape adVNTR emits:
#: insertion lengths that are 1, 2 and 0 mod 3 and one multi-digit length, a deletion
#: event, and an insertion-shaped part with no ``LEN`` token (which appears in #192's own
#: examples). Every ``LEN`` here is >= 1, which is what makes the shorthand oracle exact
#: on this set.
WELL_FORMED_PARTS = (
    "I9_2_A_LEN9",
    "I50_2_A_LEN3",
    "I22_2_G_LEN1",
    "I80_2_A_LEN2",
    "I14_2_G_LEN14",
    "D50_2",
    "I50_2",
)

#: Compounds are generated for 1..MAX_PARTS parts, each also with a trailing ``&``.
MAX_PARTS = 5

#: Parts carrying no ``LEN`` token, used to build the two unchanged classes on their own.
#: The systematic product above reaches those classes only incidentally -- a compound
#: lands in them only when *exactly one* of its parts carries a ``LEN`` and that part is
#: last -- so they are enumerated separately as well. Otherwise the "previously-parsing
#: inputs are identical" claim would rest on a few hundred probes while the changed class
#: had tens of thousands.
LEN_FREE_PARTS = ("D8_2", "D9_2", "D17_2", "I50_2", "I14_2")

#: Terminal ``LEN`` parts, including multi-digit and multiple-of-three lengths.
TERMINAL_LEN_PARTS = tuple(f"I{9 + value}_2_A_LEN{value}" for value in (1, 2, 3, 4, 5, 7, 9, 10, 12, 14, 20, 25))

#: Strings adVNTR does not emit, swept so the parser's behaviour on them is recorded
#: rather than assumed. These are where the two oracle refinements above live.
DEGENERATE_PROBES = (
    "",
    "&",
    "&&",
    "I9_2_A_LEN0",
    "I9_2_A_LEN0&D50_2",
    "I9_2_A_LEN0&I50_2_A_LEN3",
    "I9_2_A_LENX",
    "I9_2_A_LENX&I50_2_A_LEN3",
    "I9_2_A_LEN",
    "LEN",
    "LEN&LEN3",
    " I9_2_A_LEN3 ",
    "NOT_A_VARIANT",
)

#: Every state named in #192, in the task brief's measurement table, or in the golden
#: cohort. Swept explicitly so the report can quote them by name.
NAMED_PROBES = (
    "I22_2_G_LEN1",
    "I9_2_A_LEN9&I50_2_A_LEN3",
    "I9_2_A_LEN2&D50_2",
    "I9_2_A_LEN9&",
    "I50_2&I9_2_A_LEN3",
    "I50_2&D9_2&I80_2_A_LEN7",
    "D17_2&D18_2&D19_2&D20_2&D21_2",
    "I9_2_A_LEN9&I50_2_A_LEN1",
    "I9_2_A_LEN3&D50_2&D51_2",
    "I9_2_A_LEN2&D50_2&D51_2",
    "I9_2_A_LEN9&I50_2_A_LEN3&I80_2_A_LEN1",
    "D8_2&D9_2&I9_2_A_LEN9",
    "D2_2&I2_2_C_LEN5",
    "D49_2&I49_2_A_LEN12",
    "D58_2&D59_2",
    "D50_2",
)

#: The three classes the plan requires to be byte-identical.
UNCHANGED_CLASSES = ("no_len_token", "terminal_len_single_part", "terminal_len_compound")


def generate_states() -> list[str]:
    """Every ``State`` string the sweep runs, de-duplicated, in generation order.

    Returns:
        list[str]: single parts, 2--5-part compounds over :data:`WELL_FORMED_PARTS`
            (mixing insertions, deletions and ``LEN``-less insertions), each also with a
            trailing ``&``, followed by :data:`DEGENERATE_PROBES` and :data:`NAMED_PROBES`.
    """
    states: list[str] = []
    for size in range(1, MAX_PARTS + 1):
        for combination in itertools.product(WELL_FORMED_PARTS, repeat=size):
            joined = "&".join(combination)
            states.append(joined)
            states.append(f"{joined}&")
    states.extend(unchanged_class_states())
    states.extend(DEGENERATE_PROBES)
    states.extend(NAMED_PROBES)
    return list(dict.fromkeys(states))


def unchanged_class_states() -> list[str]:
    """States the oracle says must be byte-identical, enumerated in their own right.

    Two families, both of which the historic derivation handled correctly:

    * ``no_len_token`` -- every part is a deletion or a ``LEN``-less insertion;
    * ``terminal_len_*`` -- a single ``LEN`` part, last, ending the string.

    Returns:
        list[str]: the two families for 1..``MAX_PARTS`` parts.
    """
    states: list[str] = []
    for size in range(1, MAX_PARTS + 1):
        states.extend("&".join(prefix) for prefix in itertools.product(LEN_FREE_PARTS, repeat=size))
    for size in range(MAX_PARTS):
        for prefix in itertools.product(LEN_FREE_PARTS, repeat=size):
            states.extend("&".join((*prefix, terminal)) for terminal in TERMINAL_LEN_PARTS)
    return states


def historic_insertion_len(states: pd.Series) -> pd.Series:
    """The pre-#192 derivation, transcribed from ``advntr_genotyping.py`` before this change.

    Kept as a literal replica rather than a re-implementation: this is the only remaining
    copy of the behaviour the sweep measures against.

    Args:
        states (pd.Series): ``Variant`` strings.

    Returns:
        pd.Series: the historic ``Insertion_len``, as ``int``.
    """
    extracted = states.str.extract(r"(LEN.*)")[0].fillna("LEN")
    remainder = extracted.str.split("LEN", n=1, expand=True)[1]
    remainder = remainder.astype(str).replace("^$", "0", regex=True)
    return pd.to_numeric(remainder, errors="coerce").fillna(0).astype(int)


def frame_series(insertion_len: pd.Series, deletion_length: pd.Series) -> pd.Series:
    """``frame = abs(Net_indel_length)``, carried as a string as production does."""
    return (insertion_len - deletion_length).abs().astype(str)


def accepted_frames() -> tuple[set[str], set[str]]:
    """The two #182 frame series, built from the same settings the production filter reads.

    Returns:
        tuple[set[str], set[str]]: ``(ins_frame, del_frame)`` as sets of strings.
    """
    run_configuration = resolve_run_configuration()
    settings = project_advntr_settings(run_configuration.advntr, run_configuration.advntr_runtime)
    ins_frame = {str(step * settings.frameshift_multiplier + 1) for step in range(settings.max_frameshift)}
    del_frame = {str(step * settings.frameshift_multiplier + 2) for step in range(settings.max_frameshift)}
    return ins_frame, del_frame


def survival(insertion_len: pd.Series, deletion_length: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Which rows each half of the filter keeps, modelling ``advntr_processing_{ins,del}``.

    The signed rule: the insertion arm takes a **net** gain whose magnitude is in the
    ``3n+1`` series, the deletion arm a **net** loss whose magnitude is in the ``3n+2``
    series. Their union is exactly ``Net_indel_length % 3 == 1``, and they are disjoint.

    Args:
        insertion_len (pd.Series): summed (or historic) inserted length.
        deletion_length (pd.Series): ``Variant.str.count("D")``.

    Returns:
        tuple[pd.Series, pd.Series]: boolean ``(kept_by_ins, kept_by_del)``.
    """
    ins_frame, del_frame = accepted_frames()
    net = insertion_len - deletion_length
    frame = frame_series(insertion_len, deletion_length)
    kept_ins = (net > 0) & frame.isin(ins_frame)
    kept_del = (net < 0) & frame.isin(del_frame)
    return kept_ins, kept_del


def absolute_frame_survival(insertion_len: pd.Series, deletion_length: pd.Series) -> tuple[pd.Series, pd.Series]:
    """The **retired** filter: ``abs(Net_indel_length)`` against the series, presence guards only.

    Kept as a literal replica of the code before the sign repair, so the set of state
    shapes whose verdict the repair changes can be enumerated rather than argued. It is
    not a model of anything production still does.

    Args:
        insertion_len (pd.Series): summed inserted length.
        deletion_length (pd.Series): ``Variant.str.count("D")``.

    Returns:
        tuple[pd.Series, pd.Series]: boolean ``(kept_by_ins, kept_by_del)``.
    """
    ins_frame, del_frame = accepted_frames()
    frame = frame_series(insertion_len, deletion_length)
    kept_ins = (insertion_len >= 1) & frame.isin(ins_frame)
    kept_del = (deletion_length >= 1) & frame.isin(del_frame)
    return kept_ins, kept_del


def state_frame(states: list[str]) -> pd.DataFrame:
    """The batch in the shape ``pd.read_csv`` hands ``advntr_processing_{ins,del}``."""
    return pd.DataFrame(
        {
            "VID": 25561,
            "State": states,
            "NumberOfSupportingReads": 11,
            "MeanCoverage": 153.98,
            "Pvalue": 0.0001,
        }
    )


def classify(state: str) -> str:
    """Which grammar class ``state`` belongs to, by the position of its ``LEN`` tokens.

    Args:
        state (str): an adVNTR ``State`` string.

    Returns:
        str: one of ``stray_len_literal`` (a literal ``LEN`` not followed by digits, which
            adVNTR does not emit and which the historic expression anchored on),
            ``no_len_token``, ``terminal_len_single_part``, ``terminal_len_compound``
            (the only token is in the last ``&`` part and ends the string), or
            ``material_after_first_len``.
    """
    tokens = list(LEN_TOKEN.finditer(state))
    if len(LEN_LITERAL.findall(state)) != len(tokens):
        return "stray_len_literal"
    if not tokens:
        return "no_len_token"
    if tokens[0].end() != len(state):
        return "material_after_first_len"
    return "terminal_len_single_part" if "&" not in state else "terminal_len_compound"


def oracle_predicts_difference(state: str) -> bool:
    """Whether ``state`` is expected to differ, computed from the string alone.

    Stated before the sweep is run: the historic derivation yielded the intended length
    only when the text after the first literal ``LEN`` was a plain digit run, and zero
    otherwise; the new one sums every ``LEN<n>``.

    Args:
        state (str): an adVNTR ``State`` string.

    Returns:
        bool: ``True`` when the two derivations must disagree.
    """
    first = state.find("LEN")
    if first < 0:
        return False
    tail = state[first + len("LEN") :].strip()
    historic = int(tail) if tail.isdigit() else 0
    return historic != sum(int(value) for value in LEN_TOKEN.findall(state))


def cross_check_against_production(frame: pd.DataFrame, kept_ins: pd.Series, kept_del: pd.Series) -> list[str]:
    """Confirm the sweep's model of the new filter is the production filter.

    Args:
        frame (pd.DataFrame): the batch handed to the production functions.
        kept_ins (pd.Series): the sweep's model of ``advntr_processing_ins``.
        kept_del (pd.Series): the sweep's model of ``advntr_processing_del``.

    Returns:
        list[str]: human-readable mismatches; empty when the model is exact.
    """
    problems: list[str] = []
    for label, run, modelled in (
        ("advntr_processing_ins", advntr.advntr_processing_ins, kept_ins),
        ("advntr_processing_del", advntr.advntr_processing_del, kept_del),
    ):
        produced = set(run(frame.copy()).index)
        expected = set(modelled.index[modelled])
        if produced != expected:
            missing = sorted(expected - produced)[:5]
            extra = sorted(produced - expected)[:5]
            problems.append(f"{label}: model and production disagree (model-only {missing}, production-only {extra})")
    return problems


def sweep(max_examples: int = 200) -> dict[str, Any]:
    """Run both parsers over every generated state and collect the differences.

    Args:
        max_examples (int): how many differing states to record in full, per class. The
            counts, the oracle verdicts and the complete reported-call delta are always
            exhaustive; only the per-state ``before``/``after`` records are sampled, so
            the JSON stays readable when tens of thousands of states differ.

    Returns:
        dict[str, Any]: the full result, as written by ``--out``.
    """
    states = generate_states()
    frame = state_frame(states)
    variants = frame["State"]
    deletion_length = variants.str.count("D").astype(int)

    old_len = historic_insertion_len(variants)
    new_len = variants.map(advntr.sum_insertion_lengths).astype(int)
    old_frame = frame_series(old_len, deletion_length)
    new_frame = frame_series(new_len, deletion_length)
    old_ins, old_del = survival(old_len, deletion_length)
    new_ins, new_del = survival(new_len, deletion_length)

    model_problems = cross_check_against_production(frame, new_ins, new_del)

    differs = (old_len != new_len) | (old_frame != new_frame) | (old_ins != new_ins) | (old_del != new_del)
    classes = variants.map(classify)
    predicted = variants.map(oracle_predicts_difference)

    records: list[dict[str, Any]] = []
    sampled: Counter[str] = Counter()
    for position in variants.index[differs]:
        sampled[classes[position]] += 1
        # The named probes are always recorded in full: they are the states #192, the task
        # brief and the golden cohort talk about, so they must not fall off a sample.
        if sampled[classes[position]] > max_examples and variants[position] not in NAMED_PROBES:
            continue
        records.append(
            {
                "state": variants[position],
                "class": classes[position],
                "predicted_by_oracle": bool(predicted[position]),
                "deletion_length": int(deletion_length[position]),
                "before": {
                    "insertion_len": int(old_len[position]),
                    "frame": old_frame[position],
                    "kept_by_ins": bool(old_ins[position]),
                    "kept_by_del": bool(old_del[position]),
                    "reported": bool(old_ins[position] or old_del[position]),
                },
                "after": {
                    "insertion_len": int(new_len[position]),
                    "frame": new_frame[position],
                    "kept_by_ins": bool(new_ins[position]),
                    "kept_by_del": bool(new_del[position]),
                    "reported": bool(new_ins[position] or new_del[position]),
                },
            }
        )

    reported_before = old_ins | old_del
    reported_after = new_ins | new_del
    unchanged_mask = classes.isin(UNCHANGED_CLASSES)

    # The other axis: same (summed) Insertion_len, retired absolute-value filter versus the
    # signed one. This is the set of state shapes the sign repair moves, and the reason a
    # golden-cohort gate is needed at all.
    abs_ins, abs_del = absolute_frame_survival(new_len, deletion_length)
    reported_abs = abs_ins | abs_del
    net = new_len - deletion_length
    # A "pure" state changes in one direction only, so its sign is never in doubt and the
    # repair must not move it. Zero inserted bases OR zero deleted events is pure; the
    # empty state (neither) is pure by both counts and in frame under either model.
    pure = (new_len == 0) | (deletion_length == 0)
    sign_fix_moved = reported_abs != reported_after

    return {
        "probes": len(states),
        "differing": int(differs.sum()),
        "identical": int((~differs).sum()),
        "unchanged_classes": {
            "names": list(UNCHANGED_CLASSES),
            "probes": int(unchanged_mask.sum()),
            "identical": int((unchanged_mask & ~differs).sum()),
            "violations": sorted(variants[unchanged_mask & differs]),
        },
        "oracle": {
            "predicted": int(predicted.sum()),
            "differed_but_not_predicted": sorted(variants[differs & ~predicted]),
            "predicted_but_identical": sorted(variants[~differs & predicted]),
        },
        "by_class": {
            "probes": dict(Counter(classes)),
            "differing": dict(Counter(classes[differs])),
        },
        "trailing_ampersand": {
            "probes": int(variants.str.endswith("&").sum()),
            "differing": int((variants.str.endswith("&") & differs).sum()),
        },
        "example_cap_per_class": max_examples,
        "reporting_delta": {
            "gained": sorted(variants[~reported_before & reported_after]),
            "lost": sorted(variants[reported_before & ~reported_after]),
        },
        "sign_fix_delta": {
            "probes": len(states),
            "moved": int(sign_fix_moved.sum()),
            "lost": sorted(variants[reported_abs & ~reported_after]),
            "gained": sorted(variants[~reported_abs & reported_after]),
            "pure_states_moved": sorted(variants[pure & sign_fix_moved]),
            "residue_of_lost": dict(Counter(int(value) % 3 for value in net[reported_abs & ~reported_after])),
        },
        "model_problems": model_problems,
        "differences": records,
    }


def print_summary(result: dict[str, Any]) -> None:
    """Print the sweep result in the form the commit message and the report quote."""
    unchanged = result["unchanged_classes"]
    print(f"probes                     : {result['probes']}")
    print(f"identical                  : {result['identical']}")
    print(f"differing                  : {result['differing']}")
    print(f"previously-parsing probes  : {unchanged['identical']}/{unchanged['probes']} identical")
    print(f"                             (classes: {', '.join(unchanged['names'])})")
    print("\nby class:")
    for name, count in sorted(result["by_class"]["probes"].items()):
        print(f"  {name:26s} {count:6d} probes, {result['by_class']['differing'].get(name, 0):6d} differing")
    trailing = result["trailing_ampersand"]
    print(f"  {'(of which trailing &)':26s} {trailing['probes']:6d} probes, {trailing['differing']:6d} differing")
    gained = result["reporting_delta"]["gained"]
    lost = result["reporting_delta"]["lost"]
    print(f"\nreported-call delta        : +{len(gained)} gained, -{len(lost)} lost")
    for state in lost[:10]:
        print(f"  LOST  {state}")
    for state in gained[:10]:
        print(f"  GAIN  {state}")
    if len(gained) > 10:
        print(f"  ... and {len(gained) - 10} further gained states")
    sign_fix = result["sign_fix_delta"]
    print(f"\nsign-fix delta             : {sign_fix['moved']} states move verdict")
    print(f"  no longer reported         : {len(sign_fix['lost'])}")
    print(f"  newly reported             : {len(sign_fix['gained'])}")
    print(f"  net residue of the lost    : {sign_fix['residue_of_lost']}")
    # Truncated: under a regression this list can run to tens of thousands of states, and
    # the full set is always in the JSON written by --out.
    pure_moved = sign_fix["pure_states_moved"]
    print(f"  pure states moved          : {len(pure_moved)}{'' if not pure_moved else ' -- ' + str(pure_moved[:5])}")
    for state in sign_fix["lost"][:10]:
        print(f"  SIGN-LOST  {state}")
    if len(sign_fix["lost"]) > 10:
        print(f"  ... and {len(sign_fix['lost']) - 10} further states no longer reported")
    print("\noracle:")
    print(f"  predicted to differ        : {result['oracle']['predicted']}")
    print(f"  differed, not predicted    : {result['oracle']['differed_but_not_predicted']}")
    print(f"  predicted, did not differ  : {result['oracle']['predicted_but_identical']}")
    print(f"  unchanged-class violations : {unchanged['violations']}")


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv (list[str] | None): command-line arguments; ``sys.argv[1:]`` when ``None``.

    Returns:
        int: 0 when every difference is predicted and the unchanged classes are identical.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=None, help="write the full result as JSON to this path")
    parser.add_argument(
        "--max-examples",
        type=int,
        default=200,
        help="how many differing states to record in full per class (counts are always exhaustive)",
    )
    args = parser.parse_args(argv)

    result = sweep(max_examples=args.max_examples)
    print_summary(result)

    if args.out is not None:
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True))
        print(f"\nfull result written to {args.out}")

    failures = (
        result["model_problems"]
        + [f"unchanged class differed: {state!r}" for state in result["unchanged_classes"]["violations"]]
        + [f"differed but not predicted: {state!r}" for state in result["oracle"]["differed_but_not_predicted"]]
        + [f"predicted but identical: {state!r}" for state in result["oracle"]["predicted_but_identical"]]
        # The sign repair is a claim about mixed states only. A pure state changing verdict
        # would mean it altered the pure-state rule, which is not authorised.
        + [f"sign fix moved a pure state: {state!r}" for state in result["sign_fix_delta"]["pure_states_moved"]]
        # Every state the sign repair stops reporting must be one whose net change is not
        # the pathogenic frame. A residue of 1 in this bucket would mean a real call lost.
        + [
            f"sign fix dropped a pathogenic-frame state (residue {residue})"
            for residue in result["sign_fix_delta"]["residue_of_lost"]
            if residue == 1
        ]
    )
    if failures:
        print("\nREGRESSION -- the sweep did not match its oracle:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print(
        "\nOK -- every difference is predicted by the oracle; the unchanged classes are "
        "identical; the sign fix moves only mixed states, all outside the pathogenic frame."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
