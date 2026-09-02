# Mutation Testing

!!! warning "Advisory only - nothing gates on this number"

    This score is **not** a pass/fail threshold and CI does not bind against
    it. Read it as a map of which decisions are untested, not as a grade.

Line coverage answers *did a test execute this line?*. It does not answer
*would a test have noticed if this line were wrong?* - and for VNtyper the
second question is the one that matters, because the characteristic failure
is a silently wrong genotype call rather than a crash.

Mutation testing answers the second question directly: it introduces a
deliberate defect, runs the tests, and records whether anything failed. A
**surviving** mutant is a defect the suite cannot see.

## Result

**1253 of 1486 mutants killed - a raw mutation score of 84.3%.**

Of the 233 survivors, 0 are hand-classified as
*equivalent* (the mutation cannot change observable behaviour, so no test could
ever kill it) and 233 are genuine gaps. Excluding the equivalent
mutants the score is **84.3%** (1253/1486).

Both numbers are given because neither alone is honest: the raw score
understates the suite by counting unkillable mutants against it, and the
adjusted score depends on classifications that are a human judgement call.
Every classification is listed below with its reason so it can be checked.

| Module | Killed | Total | Raw score |
| --- | ---: | ---: | ---: |
| `vntyper/scripts/kestrel_genotyping.py` | 74 | 109 | 67.9% |
| `vntyper/scripts/motif_processing.py` | 44 | 61 | 72.1% |
| `vntyper/scripts/identity_candidates.py` | 140 | 184 | 76.1% |
| `vntyper/scripts/identity_reconciliation.py` | 230 | 295 | 78.0% |
| `vntyper/scripts/decision_profile.py` | 12 | 15 | 80.0% |
| `vntyper/scripts/nomenclature_dominance.py` | 79 | 96 | 82.3% |
| `vntyper/scripts/run_configuration.py` | 11 | 13 | 84.6% |
| `vntyper/scripts/variant_parsing.py` | 6 | 7 | 85.7% |
| `vntyper/scripts/decision_profile_schema.py` | 120 | 139 | 86.3% |
| `vntyper/scripts/motif_decisions.py` | 7 | 8 | 87.5% |
| `vntyper/scripts/profile_provenance.py` | 112 | 127 | 88.2% |
| `vntyper/scripts/confidence_assignment.py` | 10 | 11 | 90.9% |
| `vntyper/scripts/flagging.py` | 35 | 38 | 92.1% |
| `vntyper/scripts/calibration_statistics.py` | 168 | 178 | 94.4% |
| `vntyper/scripts/scoring.py` | 19 | 19 | 100.0% |
| `vntyper/scripts/calibration_objective.py` | 133 | 133 | 100.0% |
| `vntyper/scripts/decision_profile_semantics.py` | 53 | 53 | 100.0% |

## Surviving mutants

### Genuine gaps

Each of these is a change to the source that the **entire** unit tier does
not notice. They are the actionable output of this exercise: a test that
kills one is a test that would have caught a real defect of that shape.

| Module | Line | Mutation |
| --- | ---: | --- |
| `vntyper/scripts/calibration_statistics.py` | 85 | `or` &rarr; `and` |
| `vntyper/scripts/calibration_statistics.py` | 91 | `or` &rarr; `and` |
| `vntyper/scripts/calibration_statistics.py` | 92 | `or` &rarr; `and` |
| `vntyper/scripts/calibration_statistics.py` | 93 | `or` &rarr; `and` |
| `vntyper/scripts/calibration_statistics.py` | 94 | `or` &rarr; `and` |
| `vntyper/scripts/calibration_statistics.py` | 201 | `True` &rarr; `False` |
| `vntyper/scripts/calibration_statistics.py` | 211 | `0` &rarr; `1` |
| `vntyper/scripts/calibration_statistics.py` | 235 | `*` &rarr; `/` |
| `vntyper/scripts/calibration_statistics.py` | 237 | `1.0` &rarr; `2.0` |
| `vntyper/scripts/calibration_statistics.py` | 260 | `1` &rarr; `2` |
| `vntyper/scripts/confidence_assignment.py` | 138 | `-` &rarr; `+` |
| `vntyper/scripts/decision_profile.py` | 26 | `True` &rarr; `False` |
| `vntyper/scripts/decision_profile.py` | 60 | `*` &rarr; `/` |
| `vntyper/scripts/decision_profile.py` | 97 | `*` &rarr; `/` |
| `vntyper/scripts/decision_profile_schema.py` | 105 | `1.0` &rarr; `2.0` |
| `vntyper/scripts/decision_profile_schema.py` | 106 | `1.0` &rarr; `2.0` |
| `vntyper/scripts/decision_profile_schema.py` | 142 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 145 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 208 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 210 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 219 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 222 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 225 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 243 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 243 | `0` &rarr; `1` |
| `vntyper/scripts/decision_profile_schema.py` | 263 | `and` &rarr; `or` |
| `vntyper/scripts/decision_profile_schema.py` | 283 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 296 | `*` &rarr; `/` |
| `vntyper/scripts/decision_profile_schema.py` | 316 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 321 | `or` &rarr; `and` |
| `vntyper/scripts/decision_profile_schema.py` | 423 | `+` &rarr; `-` |
| `vntyper/scripts/decision_profile_schema.py` | 423 | `1` &rarr; `2` |
| `vntyper/scripts/decision_profile_schema.py` | 436 | `*` &rarr; `/` |
| `vntyper/scripts/flagging.py` | 148 | `True` &rarr; `False` |
| `vntyper/scripts/flagging.py` | 295 | `False` &rarr; `True` |
| `vntyper/scripts/flagging.py` | 478 | `False` &rarr; `True` |
| `vntyper/scripts/identity_candidates.py` | 53 | `True` &rarr; `False` |
| `vntyper/scripts/identity_candidates.py` | 68 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 101 | `True` &rarr; `False` |
| `vntyper/scripts/identity_candidates.py` | 121 | `True` &rarr; `False` |
| `vntyper/scripts/identity_candidates.py` | 132 | `False` &rarr; `True` |
| `vntyper/scripts/identity_candidates.py` | 156 | `True` &rarr; `False` |
| `vntyper/scripts/identity_candidates.py` | 167 | `True` &rarr; `False` |
| `vntyper/scripts/identity_candidates.py` | 311 | `-` &rarr; `+` |
| `vntyper/scripts/identity_candidates.py` | 316 | `3` &rarr; `4` |
| `vntyper/scripts/identity_candidates.py` | 482 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 495 | `*` &rarr; `/` |
| `vntyper/scripts/identity_candidates.py` | 495 | `False` &rarr; `True` |
| `vntyper/scripts/identity_candidates.py` | 498 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 498 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/identity_candidates.py` | 507 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 510 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 518 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 518 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 538 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 575 | `0` &rarr; `1` |
| `vntyper/scripts/identity_candidates.py` | 577 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 578 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 579 | `2` &rarr; `3` |
| `vntyper/scripts/identity_candidates.py` | 580 | `2` &rarr; `3` |
| `vntyper/scripts/identity_candidates.py` | 590 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 591 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 592 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 592 | `0` &rarr; `1` |
| `vntyper/scripts/identity_candidates.py` | 593 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 593 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 594 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 594 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 595 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 595 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 596 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 597 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 597 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 598 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 599 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 599 | `and` &rarr; `or` |
| `vntyper/scripts/identity_candidates.py` | 599 | `1` &rarr; `2` |
| `vntyper/scripts/identity_candidates.py` | 622 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 622 | `or` &rarr; `and` |
| `vntyper/scripts/identity_candidates.py` | 624 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 69 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 102 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 103 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 103 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 111 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 112 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 112 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 115 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 116 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 124 | `*` &rarr; `/` |
| `vntyper/scripts/identity_reconciliation.py` | 153 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 167 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 193 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 197 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 203 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 224 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 245 | `and` &rarr; `or` |
| `vntyper/scripts/identity_reconciliation.py` | 246 | `and` &rarr; `or` |
| `vntyper/scripts/identity_reconciliation.py` | 247 | `and` &rarr; `or` |
| `vntyper/scripts/identity_reconciliation.py` | 254 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 277 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 278 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 294 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 294 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 299 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 337 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 348 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/identity_reconciliation.py` | 382 | `==` &rarr; `!=` |
| `vntyper/scripts/identity_reconciliation.py` | 383 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 395 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/identity_reconciliation.py` | 406 | `*` &rarr; `/` |
| `vntyper/scripts/identity_reconciliation.py` | 445 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 456 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 474 | `continue` &rarr; `break` |
| `vntyper/scripts/identity_reconciliation.py` | 476 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/identity_reconciliation.py` | 476 | `and` &rarr; `or` |
| `vntyper/scripts/identity_reconciliation.py` | 476 | `!=` &rarr; `==` |
| `vntyper/scripts/identity_reconciliation.py` | 492 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 508 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 533 | `True` &rarr; `False` |
| `vntyper/scripts/identity_reconciliation.py` | 553 | `+=` &rarr; `-=` |
| `vntyper/scripts/identity_reconciliation.py` | 553 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 554 | `continue` &rarr; `break` |
| `vntyper/scripts/identity_reconciliation.py` | 587 | `*` &rarr; `/` |
| `vntyper/scripts/identity_reconciliation.py` | 614 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 615 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 616 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 622 | `==` &rarr; `!=` |
| `vntyper/scripts/identity_reconciliation.py` | 632 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/identity_reconciliation.py` | 672 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 778 | `*` &rarr; `/` |
| `vntyper/scripts/identity_reconciliation.py` | 818 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 827 | `and` &rarr; `or` |
| `vntyper/scripts/identity_reconciliation.py` | 827 | `>=` &rarr; `<` |
| `vntyper/scripts/identity_reconciliation.py` | 827 | `0` &rarr; `1` |
| `vntyper/scripts/identity_reconciliation.py` | 838 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 845 | `==` &rarr; `!=` |
| `vntyper/scripts/identity_reconciliation.py` | 845 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 861 | `1` &rarr; `2` |
| `vntyper/scripts/identity_reconciliation.py` | 866 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 866 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 872 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 872 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 872 | `or` &rarr; `and` |
| `vntyper/scripts/identity_reconciliation.py` | 872 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 210 | `and` &rarr; `or` |
| `vntyper/scripts/kestrel_genotyping.py` | 248 | `continue` &rarr; `break` |
| `vntyper/scripts/kestrel_genotyping.py` | 341 | `20` &rarr; `21` |
| `vntyper/scripts/kestrel_genotyping.py` | 342 | `30` &rarr; `31` |
| `vntyper/scripts/kestrel_genotyping.py` | 343 | `30` &rarr; `31` |
| `vntyper/scripts/kestrel_genotyping.py` | 352 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 433 | `break` &rarr; `continue` |
| `vntyper/scripts/kestrel_genotyping.py` | 557 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 586 | `*` &rarr; `/` |
| `vntyper/scripts/kestrel_genotyping.py` | 705 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 708 | `True` &rarr; `False` |
| `vntyper/scripts/kestrel_genotyping.py` | 729 | `True` &rarr; `False` |
| `vntyper/scripts/kestrel_genotyping.py` | 763 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/kestrel_genotyping.py` | 763 | `and` &rarr; `or` |
| `vntyper/scripts/kestrel_genotyping.py` | 763 | `or` &rarr; `and` |
| `vntyper/scripts/kestrel_genotyping.py` | 763 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 773 | `/` &rarr; `*` |
| `vntyper/scripts/kestrel_genotyping.py` | 774 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 775 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 775 | `5` &rarr; `6` |
| `vntyper/scripts/kestrel_genotyping.py` | 787 | `True` &rarr; `False` |
| `vntyper/scripts/kestrel_genotyping.py` | 790 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 957 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 1080 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1095 | `*` &rarr; `/` |
| `vntyper/scripts/kestrel_genotyping.py` | 1177 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1181 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1183 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1186 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1188 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1191 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1193 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 1200 | `1` &rarr; `2` |
| `vntyper/scripts/kestrel_genotyping.py` | 1325 | `==` &rarr; `!=` |
| `vntyper/scripts/kestrel_genotyping.py` | 1325 | `1` &rarr; `2` |
| `vntyper/scripts/motif_decisions.py` | 86 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 217 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 232 | `==` &rarr; `!=` |
| `vntyper/scripts/motif_processing.py` | 240 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 291 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 293 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 296 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 331 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 338 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 342 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 343 | `60` &rarr; `61` |
| `vntyper/scripts/motif_processing.py` | 368 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 494 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 513 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 513 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 526 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 528 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 537 | `False` &rarr; `True` |
| `vntyper/scripts/nomenclature_dominance.py` | 67 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 68 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 69 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 70 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 71 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 72 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 82 | `or` &rarr; `and` |
| `vntyper/scripts/nomenclature_dominance.py` | 90 | `True` &rarr; `False` |
| `vntyper/scripts/nomenclature_dominance.py` | 117 | `or` &rarr; `and` |
| `vntyper/scripts/nomenclature_dominance.py` | 122 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 158 | `or` &rarr; `and` |
| `vntyper/scripts/nomenclature_dominance.py` | 191 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 191 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 191 | `2` &rarr; `3` |
| `vntyper/scripts/nomenclature_dominance.py` | 194 | `0` &rarr; `1` |
| `vntyper/scripts/nomenclature_dominance.py` | 194 | `in` &rarr; `not in` |
| `vntyper/scripts/nomenclature_dominance.py` | 196 | `!=` &rarr; `==` |
| `vntyper/scripts/profile_provenance.py` | 105 | `True` &rarr; `False` |
| `vntyper/scripts/profile_provenance.py` | 122 | `True` &rarr; `False` |
| `vntyper/scripts/profile_provenance.py` | 128 | `and` &rarr; `or` |
| `vntyper/scripts/profile_provenance.py` | 130 | `True` &rarr; `False` |
| `vntyper/scripts/profile_provenance.py` | 135 | `*` &rarr; `/` |
| `vntyper/scripts/profile_provenance.py` | 168 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 213 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 236 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 237 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 243 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 255 | `continue` &rarr; `break` |
| `vntyper/scripts/profile_provenance.py` | 276 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 289 | `*` &rarr; `/` |
| `vntyper/scripts/profile_provenance.py` | 381 | `or` &rarr; `and` |
| `vntyper/scripts/profile_provenance.py` | 387 | `or` &rarr; `and` |
| `vntyper/scripts/run_configuration.py` | 30 | `True` &rarr; `False` |
| `vntyper/scripts/run_configuration.py` | 97 | `*` &rarr; `/` |
| `vntyper/scripts/variant_parsing.py` | 114 | `0.0` &rarr; `1.0` |

### Classified equivalent

None classified yet.

## How this compares to the 43.5% baseline

The experiment that motivated this work scored **43.5%** (27 of 62 mutants
killed) across the eight highest-coverage modules. That harness was never
committed, which is why this one exists.

!!! note "The two totals are not directly comparable"

    Different mutant population, different modules: 62 mutants over eight
    modules then, 1486 over 17 modules now, generated by a different operator set.
    A higher or lower headline number would not by itself mean the suite has
    improved or regressed. Only per-module figures on the same module carry
    across, and even those only loosely.

The one comparison that is meaningful is `confidence_assignment.py`, the
module that motivated the whole effort: it had **100% line coverage and a 21%
mutation score**, i.e. four of five deliberate defects in it went undetected by
a fully green build.

| `confidence_assignment.py` | Then | Now |
| --- | ---: | ---: |
| Line coverage | 100% | 100% |
| Mutation score | 21% | 90.9% raw, 90.9% adjusted |

## Reproducing this

```bash
make mutation
```

The harness is `scripts/mutation_test.py`. It mutates one token at a time,
runs the module's own tests first and escalates anything that survives to
the full unit tier before recording it as a survivor, so the score is not
biased by the scoping.

Before it mutates anything it runs those same pytest invocations against the
**unmutated** tree and aborts unless they pass, printing the failure. A mutant
is counted as killed whenever pytest exits non-zero, and pytest exits non-zero
for an unrelated failure, a collection error or a missing dependency just as
readily - so without that preflight a broken checkout scores 100% and
overwrites this page with the result.

!!! danger "No child may load bytecode built from a different revision"

    That is the invariant. CPython validates a cached `.pyc` against the
    source's `(mtime, size)` pair with one-second mtime granularity, so a
    mutant written in the same second as the file it replaces **and of the
    same byte length** (`==` to `!=`, `1` to `2`) is indistinguishable from
    it to the cache validator: the interpreter loads the stale `.pyc`, runs
    the **unmutated** code, every such mutant *survives* and the score is
    fiction. Two sweeps produced exactly that before it was found.

    Two defences hold it, and both are required. `run_pytest()` passes
    `python -B` and sets `PYTHONDONTWRITEBYTECODE=1` in the child, so no
    `.pyc` is written during the run; and every `__pycache__` under
    `vntyper/` is deleted before the sweep starts and again after each
    mutant is written, so none left by an earlier run or an earlier mutant
    can be loaded. The flags stop new caches, the deletion stops old ones.

    The `PYTHONDONTWRITEBYTECODE=1` on the `make mutation` recipe is defence
    in depth for the parent process, which never imports a target module -
    it is not what holds the invariant, and the harness is safe without it.

!!! note "Each measurement runs in an isolated workspace"

    The harness captures HEAD in a disposable detached worktree and overlays
    the current non-ignored working state, except selected mutation targets
    and requested output paths. Selected targets therefore come from the
    captured commit, while ordinary edits and new tests participate in the
    measurement without being written back.

    Import provenance is proved against the pinned worktree before testing.
    A green baseline and a known-killed canary must then pass before ordinary
    mutants are measured, and the post-overlay baseline is verified after
    the canary and after every target.

    Every mutant and bytecode-cache write is confined to that workspace;
    real production source is never mutated. Requested report artifacts are
    built completely and installed atomically in the real checkout.
    The cleanup is best effort: SIGINT, SIGTERM, SIGHUP and SIGQUIT attempt
    the common unwind path, while SIGKILL or a host crash can leave only an
    orphan disposable worktree for later inspection and removal.

## Related: branch coverage, now enabled

Mutation testing and branch coverage were investigated together, because both
ask a sharper question than line coverage and they agreed on which modules are
weakest. The branch-coverage half of that work is recorded here so it is not
re-derived from scratch.

`[tool.coverage.run]` now sets `branch = true`, so an `if` that is entered but
never taken no longer counts as fully covered. It was enabled in **#196**,
measured on `fix/issue-181-197-followups` at `5bb2463`:

| Measure | Value |
| --- | ---: |
| Line (statement) coverage | 76.60% |
| **Branch-inclusive total** | **74.22%** |
| Branch-only coverage | 66.00% |
| Branch exits never taken | 512 of 1506 |

`fail_under` was raised **70 &rarr; 74** in the same commit, to the figure
`scripts/coverage_gate.py` printed for that run. **The floor was raised to meet
the measurement; the measurement was not weakened to fit the floor.** That
distinction is the whole point of the ratchet.

!!! warning "74 is a branch-inclusive floor, and nothing else notices if that changes"

    Deleting `branch = true` does not fail any gate on its own - it *raises*
    the reported total, because statement-only coverage of the same suite is
    76.60% against the branch-inclusive 74.22%. The build would go green while
    measuring strictly less. `tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled`
    exists solely to make that edit fail, and it is the only thing that does.

### Correction: the previously recorded prerequisite was wrong

This section formerly recorded the opposite conclusion, and it is kept here
rather than deleted, because a document that quietly rewrites its own history
stops being worth trusting.

The earlier measurement was **63.80%** branch-inclusive (76.60% is the current
line figure; it was 66.82% then) against a `fail_under` of **66**. Enabling
branch coverage at that point really would have failed CI on the enabling
commit, and the decision not to enable it - and specifically not to lower the
floor to admit it - was correct.

What was wrong was the stated route out. The old text identified
`cohort_summary.py` and `install_references.py` as holding 275 of the 685
untaken exits, both on the oversized-file list in `AGENTS.md`, and concluded:
*"Splitting them is the prerequisite, not writing more tests against them."*

**It was not a prerequisite.** Branch coverage cleared the floor with both
files still unsplit and still untested. The gap was closed instead by testing
five small, already-testable modules to 100% and a sixth to 98%:
`cross_match.py`, `utils.py`, `file_processing.py`,
`extract_unmapped_from_offset.py`, `variant_parsing.py`, and `docker/app/tasks.py`.

The generalisable mistake was reading "these two files hold the most untaken
exits" as "these two files are the ones that must be fixed". Concentration of
missing coverage is not the same as cheapness of covering it: the two oversized
modules are expensive precisely because they fuse I/O with logic, while the
same number of units was available across several modules that could be called
directly. Splitting `cohort_summary.py` and `install_references.py` remains
worth doing for the reasons `AGENTS.md` gives - it was simply never a blocker
for this.

## Raw output

```text
VNtyper mutation testing - advisory score
============================================================

Command:  make mutation
Total:    1486 mutants, 1253 killed, 233 survived
Score:    84.3%
Duration: 190.9 min

Per module
------------------------------------------------------------
  67.9%   74/109  vntyper/scripts/kestrel_genotyping.py
  72.1%   44/ 61  vntyper/scripts/motif_processing.py
  76.1%  140/184  vntyper/scripts/identity_candidates.py
  78.0%  230/295  vntyper/scripts/identity_reconciliation.py
  80.0%   12/ 15  vntyper/scripts/decision_profile.py
  82.3%   79/ 96  vntyper/scripts/nomenclature_dominance.py
  84.6%   11/ 13  vntyper/scripts/run_configuration.py
  85.7%    6/  7  vntyper/scripts/variant_parsing.py
  86.3%  120/139  vntyper/scripts/decision_profile_schema.py
  87.5%    7/  8  vntyper/scripts/motif_decisions.py
  88.2%  112/127  vntyper/scripts/profile_provenance.py
  90.9%   10/ 11  vntyper/scripts/confidence_assignment.py
  92.1%   35/ 38  vntyper/scripts/flagging.py
  94.4%  168/178  vntyper/scripts/calibration_statistics.py
 100.0%   19/ 19  vntyper/scripts/scoring.py
 100.0%  133/133  vntyper/scripts/calibration_objective.py
 100.0%   53/ 53  vntyper/scripts/decision_profile_semantics.py

Surviving mutants  [E] = hand-classified equivalent, [ ] = genuine gap
------------------------------------------------------------
vntyper/scripts/calibration_statistics.py
  [ ] line   85  'or' -> 'and'
  [ ] line   91  'or' -> 'and'
  [ ] line   92  'or' -> 'and'
  [ ] line   93  'or' -> 'and'
  [ ] line   94  'or' -> 'and'
  [ ] line  201  'True' -> 'False'
  [ ] line  211  '0' -> '1'
  [ ] line  235  '*' -> '/'
  [ ] line  237  '1.0' -> '2.0'
  [ ] line  260  '1' -> '2'

vntyper/scripts/confidence_assignment.py
  [ ] line  138  '-' -> '+'

vntyper/scripts/decision_profile.py
  [ ] line   26  'True' -> 'False'
  [ ] line   60  '*' -> '/'
  [ ] line   97  '*' -> '/'

vntyper/scripts/decision_profile_schema.py
  [ ] line  105  '1.0' -> '2.0'
  [ ] line  106  '1.0' -> '2.0'
  [ ] line  142  'or' -> 'and'
  [ ] line  145  'or' -> 'and'
  [ ] line  208  'or' -> 'and'
  [ ] line  210  'or' -> 'and'
  [ ] line  219  'or' -> 'and'
  [ ] line  222  'or' -> 'and'
  [ ] line  225  'or' -> 'and'
  [ ] line  243  'or' -> 'and'
  [ ] line  243  '0' -> '1'
  [ ] line  263  'and' -> 'or'
  [ ] line  283  'or' -> 'and'
  [ ] line  296  '*' -> '/'
  [ ] line  316  'or' -> 'and'
  [ ] line  321  'or' -> 'and'
  [ ] line  423  '+' -> '-'
  [ ] line  423  '1' -> '2'
  [ ] line  436  '*' -> '/'

vntyper/scripts/flagging.py
  [ ] line  148  'True' -> 'False'
  [ ] line  295  'False' -> 'True'
  [ ] line  478  'False' -> 'True'

vntyper/scripts/identity_candidates.py
  [ ] line   53  'True' -> 'False'
  [ ] line   68  '1' -> '2'
  [ ] line  101  'True' -> 'False'
  [ ] line  121  'True' -> 'False'
  [ ] line  132  'False' -> 'True'
  [ ] line  156  'True' -> 'False'
  [ ] line  167  'True' -> 'False'
  [ ] line  311  '-' -> '+'
  [ ] line  316  '3' -> '4'
  [ ] line  482  'or' -> 'and'
  [ ] line  495  '*' -> '/'
  [ ] line  495  'False' -> 'True'
  [ ] line  498  'or' -> 'and'
  [ ] line  498  'not' -> ''
  [ ] line  507  '1' -> '2'
  [ ] line  510  'or' -> 'and'
  [ ] line  518  'or' -> 'and'
  [ ] line  518  'or' -> 'and'
  [ ] line  538  '1' -> '2'
  [ ] line  575  '0' -> '1'
  [ ] line  577  '1' -> '2'
  [ ] line  578  '1' -> '2'
  [ ] line  579  '2' -> '3'
  [ ] line  580  '2' -> '3'
  [ ] line  590  'and' -> 'or'
  [ ] line  591  'and' -> 'or'
  [ ] line  592  'and' -> 'or'
  [ ] line  592  '0' -> '1'
  [ ] line  593  'and' -> 'or'
  [ ] line  593  '1' -> '2'
  [ ] line  594  'and' -> 'or'
  [ ] line  594  '1' -> '2'
  [ ] line  595  'and' -> 'or'
  [ ] line  595  'and' -> 'or'
  [ ] line  596  'and' -> 'or'
  [ ] line  597  'and' -> 'or'
  [ ] line  597  '1' -> '2'
  [ ] line  598  'and' -> 'or'
  [ ] line  599  'and' -> 'or'
  [ ] line  599  'and' -> 'or'
  [ ] line  599  '1' -> '2'
  [ ] line  622  'or' -> 'and'
  [ ] line  622  'or' -> 'and'
  [ ] line  624  '1' -> '2'

vntyper/scripts/identity_reconciliation.py
  [ ] line   69  'True' -> 'False'
  [ ] line  102  'or' -> 'and'
  [ ] line  103  'or' -> 'and'
  [ ] line  103  'or' -> 'and'
  [ ] line  111  'or' -> 'and'
  [ ] line  112  'or' -> 'and'
  [ ] line  112  'or' -> 'and'
  [ ] line  115  'or' -> 'and'
  [ ] line  116  'or' -> 'and'
  [ ] line  124  '*' -> '/'
  [ ] line  153  'or' -> 'and'
  [ ] line  167  'True' -> 'False'
  [ ] line  193  'or' -> 'and'
  [ ] line  197  'or' -> 'and'
  [ ] line  203  'or' -> 'and'
  [ ] line  224  'or' -> 'and'
  [ ] line  245  'and' -> 'or'
  [ ] line  246  'and' -> 'or'
  [ ] line  247  'and' -> 'or'
  [ ] line  254  'True' -> 'False'
  [ ] line  277  'or' -> 'and'
  [ ] line  278  'or' -> 'and'
  [ ] line  294  'or' -> 'and'
  [ ] line  294  'or' -> 'and'
  [ ] line  299  'or' -> 'and'
  [ ] line  337  'or' -> 'and'
  [ ] line  348  'not' -> ''
  [ ] line  382  '==' -> '!='
  [ ] line  383  '1' -> '2'
  [ ] line  395  'not' -> ''
  [ ] line  406  '*' -> '/'
  [ ] line  445  'or' -> 'and'
  [ ] line  456  'True' -> 'False'
  [ ] line  474  'continue' -> 'break'
  [ ] line  476  'not' -> ''
  [ ] line  476  'and' -> 'or'
  [ ] line  476  '!=' -> '=='
  [ ] line  492  'True' -> 'False'
  [ ] line  508  'True' -> 'False'
  [ ] line  533  'True' -> 'False'
  [ ] line  553  '+=' -> '-='
  [ ] line  553  '1' -> '2'
  [ ] line  554  'continue' -> 'break'
  [ ] line  587  '*' -> '/'
  [ ] line  614  'or' -> 'and'
  [ ] line  615  'or' -> 'and'
  [ ] line  616  'or' -> 'and'
  [ ] line  622  '==' -> '!='
  [ ] line  632  'not' -> ''
  [ ] line  672  '1' -> '2'
  [ ] line  778  '*' -> '/'
  [ ] line  818  '1' -> '2'
  [ ] line  827  'and' -> 'or'
  [ ] line  827  '>=' -> '<'
  [ ] line  827  '0' -> '1'
  [ ] line  838  'or' -> 'and'
  [ ] line  845  '==' -> '!='
  [ ] line  845  '1' -> '2'
  [ ] line  861  '1' -> '2'
  [ ] line  866  'or' -> 'and'
  [ ] line  866  'or' -> 'and'
  [ ] line  872  'or' -> 'and'
  [ ] line  872  'or' -> 'and'
  [ ] line  872  'or' -> 'and'
  [ ] line  872  '0' -> '1'

vntyper/scripts/kestrel_genotyping.py
  [ ] line  210  'and' -> 'or'
  [ ] line  248  'continue' -> 'break'
  [ ] line  341  '20' -> '21'
  [ ] line  342  '30' -> '31'
  [ ] line  343  '30' -> '31'
  [ ] line  352  'False' -> 'True'
  [ ] line  433  'break' -> 'continue'
  [ ] line  557  'False' -> 'True'
  [ ] line  586  '*' -> '/'
  [ ] line  705  '0' -> '1'
  [ ] line  708  'True' -> 'False'
  [ ] line  729  'True' -> 'False'
  [ ] line  763  'not' -> ''
  [ ] line  763  'and' -> 'or'
  [ ] line  763  'or' -> 'and'
  [ ] line  763  'False' -> 'True'
  [ ] line  773  '/' -> '*'
  [ ] line  774  'False' -> 'True'
  [ ] line  775  '0' -> '1'
  [ ] line  775  '5' -> '6'
  [ ] line  787  'True' -> 'False'
  [ ] line  790  'False' -> 'True'
  [ ] line  957  'False' -> 'True'
  [ ] line 1080  '0' -> '1'
  [ ] line 1095  '*' -> '/'
  [ ] line 1177  '0' -> '1'
  [ ] line 1181  '0' -> '1'
  [ ] line 1183  '0' -> '1'
  [ ] line 1186  '0' -> '1'
  [ ] line 1188  '0' -> '1'
  [ ] line 1191  '0' -> '1'
  [ ] line 1193  '0' -> '1'
  [ ] line 1200  '1' -> '2'
  [ ] line 1325  '==' -> '!='
  [ ] line 1325  '1' -> '2'

vntyper/scripts/motif_decisions.py
  [ ] line   86  'True' -> 'False'

vntyper/scripts/motif_processing.py
  [ ] line  217  'False' -> 'True'
  [ ] line  232  '==' -> '!='
  [ ] line  240  'True' -> 'False'
  [ ] line  291  'False' -> 'True'
  [ ] line  293  'False' -> 'True'
  [ ] line  296  'False' -> 'True'
  [ ] line  331  'False' -> 'True'
  [ ] line  338  'True' -> 'False'
  [ ] line  342  'False' -> 'True'
  [ ] line  343  '60' -> '61'
  [ ] line  368  'True' -> 'False'
  [ ] line  494  'True' -> 'False'
  [ ] line  513  '-' -> '+'
  [ ] line  513  '1' -> '2'
  [ ] line  526  'not' -> ''
  [ ] line  528  'not' -> ''
  [ ] line  537  'False' -> 'True'

vntyper/scripts/nomenclature_dominance.py
  [ ] line   67  '0' -> '1'
  [ ] line   68  '0' -> '1'
  [ ] line   69  '0' -> '1'
  [ ] line   70  '0' -> '1'
  [ ] line   71  '0' -> '1'
  [ ] line   72  '0' -> '1'
  [ ] line   82  'or' -> 'and'
  [ ] line   90  'True' -> 'False'
  [ ] line  117  'or' -> 'and'
  [ ] line  122  '0' -> '1'
  [ ] line  158  'or' -> 'and'
  [ ] line  191  '0' -> '1'
  [ ] line  191  '0' -> '1'
  [ ] line  191  '2' -> '3'
  [ ] line  194  '0' -> '1'
  [ ] line  194  'in' -> 'not in'
  [ ] line  196  '!=' -> '=='

vntyper/scripts/profile_provenance.py
  [ ] line  105  'True' -> 'False'
  [ ] line  122  'True' -> 'False'
  [ ] line  128  'and' -> 'or'
  [ ] line  130  'True' -> 'False'
  [ ] line  135  '*' -> '/'
  [ ] line  168  'or' -> 'and'
  [ ] line  213  'or' -> 'and'
  [ ] line  236  'or' -> 'and'
  [ ] line  237  'or' -> 'and'
  [ ] line  243  'or' -> 'and'
  [ ] line  255  'continue' -> 'break'
  [ ] line  276  'or' -> 'and'
  [ ] line  289  '*' -> '/'
  [ ] line  381  'or' -> 'and'
  [ ] line  387  'or' -> 'and'

vntyper/scripts/run_configuration.py
  [ ] line   30  'True' -> 'False'
  [ ] line   97  '*' -> '/'

vntyper/scripts/variant_parsing.py
  [ ] line  114  '0.0' -> '1.0'

```
