"""The cross-assembly coverage claim, measured against the paired fixtures.

`tests/data/remapped/bwa/` holds seven samples aligned to both GRCh37 and GRCh38 -
same reads, same molecules, two references. That makes the comparison paired, with
no sample-level confounding, and it is the only place in this repository where the
claim behind #222 can be checked rather than asserted.

**What this test is, and is not.** It is a regression guard on this pipeline's
behaviour under `DEPTH_COUNTING_POLICY`. It is *not* evidence that the flank figure
is a general biological invariant, and it must not be cited as validation:

* All seven fixtures began as hg19 MUC1-region subsets and every GRCh37/GRCh38 pair
  was remapped from the same `example_*_hg19_subset_R{1,2}.fastq.gz`. The read pool
  was therefore ascertained by hg19 alignment *before* the comparison.
* Duplicate flags were lost in that remapping - `example_6449`'s original BAM carries
  65,867 duplicate-marked records and its remapped counterpart none - so these files
  do not represent what the duplicate-marking path produces either.
* One aligner, one capture design, no sample with an independently established repeat
  count.

Establishing generality needs whole-library alignments to full references from data
of both assemblies' own origin. Until then the honest claim is the narrow one this
test pins.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from vntyper.scripts.coverage_stats import VNTR_FLANK_BASES, vntr_geometry

pytestmark = pytest.mark.integration

REPOSITORY = Path(__file__).resolve().parents[2]
REMAPPED = REPOSITORY / "tests/data/remapped/bwa"
SAMPLES = ("6449", "66bf", "6c28", "7a61", "a5c1", "b178", "dfc3")

#: Both builds' fixtures use UCSC naming, so the contig is the same string.
CONTIG = "chr1"

#: Measured spread of the flank ratio across these fixtures is 0.976-1.001. 5% absorbs
#: that with room to spare and still fails on an asymmetric array-bound error, which is
#: the failure mode that matters: a one-sided 140 bp error moves the comparison by 17.5%.
TOLERANCE = 0.05

ASSEMBLIES = {
    "hg19": {
        "chromosome": 1,
        "bam_region_coords": "155158000-155163000",
        "vntr_region_coords": "155160500-155162000",
        "vntr_array_coords": "155161000-155161810",
    },
    "hg38": {
        "chromosome": 1,
        "bam_region_coords": "155184000-155194000",
        "vntr_region_coords": "155188000-155192500",
        "vntr_array_coords": "155188530-155192010",
    },
}


def _depths(bam: Path, start: int, end: int) -> dict[int, int]:
    """Per-base depth under exactly the policy the pipeline runs.

    No `-q`/`-Q`: the array figures are only conserved because every primary alignment
    is counted whatever its mapping quality, so a filtering flag here would silently
    measure a different quantity from the one the pipeline reports.
    """
    completed = subprocess.run(
        ["samtools", "depth", "-a", "-r", f"{CONTIG}:{start}-{end}", str(bam)],
        capture_output=True,
        text=True,
        check=True,
    )
    depths = {}
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        _, position, depth = line.split("\t")
        depths[int(position)] = int(depth)
    return depths


def _figures(sample: str, build: str) -> tuple[float, float]:
    """Return (window mean, flank mean) for one sample on one build."""
    bam = REMAPPED / build / f"example_{sample}_{build}_bwa.bam"
    if not bam.is_file():
        pytest.skip(f"paired fixture {bam} is not present")

    geometry = vntr_geometry(ASSEMBLIES[build])
    assert geometry is not None, f"{build} fixture config must carry vntr_array_coords"

    depths = _depths(bam, *geometry.window)
    window = [depths.get(position, 0) for position in range(geometry.window[0], geometry.window[1] + 1)]
    flank = [depths.get(position, 0) for start, end in geometry.flank for position in range(start, end + 1)]
    return sum(window) / len(window), sum(flank) / len(flank)


@pytest.fixture(scope="module", autouse=True)
def _requires_samtools() -> None:
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not on PATH")


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_flank_depth_is_the_same_on_both_assemblies(sample: str) -> None:
    """The figure the report labels as comparable actually is."""
    _, flank37 = _figures(sample, "hg19")
    _, flank38 = _figures(sample, "hg38")

    ratio = flank37 / flank38

    assert abs(ratio - 1.0) <= TOLERANCE, (
        f"{sample}: flank depth differs between assemblies by {abs(ratio - 1.0) * 100:.1f}% "
        f"({flank37:.1f} on GRCh37 against {flank38:.1f} on GRCh38). The flank is unique "
        f"sequence and should measure the sample, not the reference - so this most likely "
        f"means the {VNTR_FLANK_BASES} bp flanks are no longer symmetric about the array."
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_window_mean_is_not_comparable_and_that_is_the_defect(sample: str) -> None:
    """The other half of the claim, pinned so a later change cannot quietly erase it.

    Without this, someone "fixing" #222 by harmonising the two configured windows
    would see the flank test still pass and conclude the problem was solved. It would
    not be: the difference is in how much array each assembly represents, and no
    window can change that.
    """
    mean37, _ = _figures(sample, "hg19")
    mean38, _ = _figures(sample, "hg38")

    ratio = mean37 / mean38

    assert ratio > 1.0 + TOLERANCE, (
        f"{sample}: window mean is now comparable between assemblies (ratio {ratio:.3f}). "
        f"That contradicts the measurement behind #222, where GRCh37 ran about 2.7x higher "
        f"over the window. Either the fixtures changed or the configured windows moved - "
        f"and if the windows moved, the build-comparable columns need re-deriving together."
    )
