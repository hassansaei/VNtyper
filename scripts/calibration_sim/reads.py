"""Bounded deterministic paired-read simulation with independent origin truth."""

from __future__ import annotations

import logging
import math
import random
from collections.abc import Iterator
from dataclasses import dataclass
from typing import NoReturn

from .haplotypes import SimulatedHaplotype, diploid_truth

logger = logging.getLogger(__name__)
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_ALTERNATES = {base: tuple(other for other in "ACGT" if other != base) for base in "ACGT"}


@dataclass(frozen=True)
class SimulatedReadPair:
    """Forward FASTQ sequences with an independent anonymous fragment origin."""

    name: str
    allele_index: int
    fragment_start: int
    fragment_end: int
    read1: str
    read2: str
    qualities: str
    read1_substitutions: int
    read2_substitutions: int


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _integer(value: object, name: str, minimum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(f"read simulation {name} must be an integer >= {minimum}")


def _nonnegative(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"read simulation {name} must be finite and nonnegative")
    try:
        number = float(value)
    except OverflowError:
        _fail(f"read simulation {name} must be finite and nonnegative")
    if not math.isfinite(number) or number < 0:
        _fail(f"read simulation {name} must be finite and nonnegative")
    return number


def _substitute(sequence: str, rate: float, generator: random.Random) -> tuple[str, int]:
    if rate == 0:
        return sequence, 0
    bases = []
    count = 0
    for base in sequence:
        draw = generator.random()
        alternate = generator.choice(_ALTERNATES[base])
        if draw < rate:
            bases.append(alternate)
            count += 1
        else:
            bases.append(base)
    return "".join(bases), count


def _pairs(
    haplotypes: tuple[SimulatedHaplotype, SimulatedHaplotype],
    pair_count: int,
    read_length: int,
    fragment_length: int,
    seed: int,
    rate: float,
    quality_score: int,
    starts: tuple[int, int],
    first_allele_probability: float,
) -> Iterator[SimulatedReadPair]:
    generator = random.Random(seed)
    errors = random.Random(f"calibration-read-errors-v1:{seed}")
    qualities = chr(33 + quality_score) * read_length
    for ordinal in range(pair_count):
        allele = 0 if generator.random() < first_allele_probability else 1
        start = generator.randrange(starts[allele])
        end = start + fragment_length
        source = haplotypes[allele].sequence
        first, first_errors = _substitute(source[start : start + read_length], rate, errors)
        second, second_errors = _substitute(source[end - read_length : end].translate(_COMPLEMENT)[::-1], rate, errors)
        yield SimulatedReadPair(
            f"synthetic-read-{ordinal:012d}", allele, start, end, first, second, qualities, first_errors, second_errors
        )


def generate_read_pairs(
    *,
    haplotypes: tuple[SimulatedHaplotype, SimulatedHaplotype],
    pair_count: int,
    read_length: int,
    fragment_length: int,
    seed: int,
    substitution_rate: float,
    quality_score: int,
    allele_copy_weights: tuple[float, float],
    maximum_generated_bases: int,
) -> Iterator[SimulatedReadPair]:
    """Stream paired reads sampled uniformly over weighted valid fragment starts.

    Equal copy weights give every possible fragment start the same probability;
    a longer allele therefore contributes more pairs. Zero copy weight explicitly
    simulates dropout. The same seed and arguments reproduce reads within the
    recorded Python random-generator implementation. Fragment origins and errors
    use separate streams, so error-rate comparisons retain identical origins and
    nested substitutions. No process-global RNG is used.

    Args:
        haplotypes: Two validated known-truth sequence constructions.
        pair_count: Exact number of paired fragments; zero models absent coverage.
        read_length: Bases emitted per mate, at most the fragment length.
        fragment_length: Fixed physical fragment length, supported by both alleles.
        seed: Nonnegative deterministic seed recorded by the simulation protocol.
        substitution_rate: Independent per-base substitution probability in [0, 1].
        quality_score: Explicit constant Phred+33 score in [0, 93], independently
            controlled so a protocol can test mismatched quality/error assumptions.
        allele_copy_weights: Relative nonnegative copy/recovery weights, not read fractions.
        maximum_generated_bases: Explicit budget for pair_count * 2 * read_length.

    Returns:
        A bounded-memory iterator of sequences and anonymous fragment truth.

    Raises:
        ValueError: Before iteration if geometry, parameters or budget are invalid.
    """
    if not isinstance(haplotypes, tuple) or len(haplotypes) != 2:
        _fail("read simulation requires exactly two immutable haplotypes")
    diploid_truth(*haplotypes)
    for value, name, minimum in (
        (pair_count, "pair count", 0),
        (read_length, "read length", 1),
        (fragment_length, "fragment length", 1),
        (seed, "seed", 0),
        (quality_score, "quality score", 0),
        (maximum_generated_bases, "output base budget", 0),
    ):
        _integer(value, name, minimum)
    if quality_score > 93:
        _fail("read simulation quality score exceeds Phred+33 range")
    if read_length > fragment_length or any(fragment_length > len(item.sequence) for item in haplotypes):
        _fail("read simulation fragment geometry is incompatible with reads or haplotypes")
    if pair_count * 2 * read_length > maximum_generated_bases:
        _fail("read simulation exceeds its declared output base budget")
    rate = _nonnegative(substitution_rate, "substitution rate")
    if rate > 1:
        _fail("read simulation substitution rate must be <= 1")
    if not isinstance(allele_copy_weights, tuple) or len(allele_copy_weights) != 2:
        _fail("read simulation copy weights must be an immutable pair")
    first_weight, second_weight = (_nonnegative(value, "copy weight") for value in allele_copy_weights)
    maximum_weight = max(first_weight, second_weight)
    if maximum_weight == 0:
        _fail("read simulation requires at least one nonzero copy weight")
    starts = tuple(len(item.sequence) - fragment_length + 1 for item in haplotypes)
    first_mass = starts[0] * (first_weight / maximum_weight)
    second_mass = starts[1] * (second_weight / maximum_weight)
    return _pairs(
        haplotypes,
        pair_count,
        read_length,
        fragment_length,
        seed,
        rate,
        quality_score,
        (starts[0], starts[1]),
        first_mass / (first_mass + second_mass),
    )
