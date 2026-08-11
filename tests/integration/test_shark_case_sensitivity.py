"""Integration coverage of SHARK's case sensitivity against soft-masked references.

Moved out of ``tests/unit/test_shark_filtering.py``: that module applies
``pytestmark = pytest.mark.unit`` at module scope, so an integration test added there
would carry *both* markers and still be selected by ``pytest -m unit``. This test needs
the real ``shark`` binary on ``PATH`` (the ``shark_env`` conda environment), so it belongs
in the integration tier and is skipped outright when that binary is absent.
"""

import pytest

pytestmark = pytest.mark.integration


class TestSharkIsCaseInsensitive:
    """Verified against shark 1.2.0 (h077b44d_5): its to_int table maps a/c/g/t and
    A/C/G/T identically, which is why the soft masking in the derived references is
    harmless. A future base rebuild could change that, and ordinary cohort reads need
    not contain a lowercase-only discriminating k-mer to catch it.
    """

    def test_a_lowercase_only_match_is_retained_exactly_as_its_uppercase_copy(self, tmp_path):
        import shutil
        import subprocess

        shark = shutil.which("shark")
        if shark is None:
            pytest.skip("shark binary not on PATH")

        # Every 17-mer overlapping the read is lowercase in `lower`, uppercase in `upper`.
        sequence = "acgt" * 40
        lower = tmp_path / "lower.fa"
        upper = tmp_path / "upper.fa"
        lower.write_text(f">r\n{sequence}\n")
        upper.write_text(f">r\n{sequence.upper()}\n")

        read = sequence[20:100].upper()  # reads are always uppercase
        quality = "I" * len(read)
        for name in ("r1", "r2"):
            (tmp_path / f"{name}.fastq").write_text(f"@read1\n{read}\n+\n{quality}\n")

        def retained(reference):
            out1 = tmp_path / f"{reference.stem}_1.fq"
            out2 = tmp_path / f"{reference.stem}_2.fq"
            subprocess.run(
                [
                    shark,
                    "-r",
                    str(reference),
                    "-1",
                    str(tmp_path / "r1.fastq"),
                    "-2",
                    str(tmp_path / "r2.fastq"),
                    "-o",
                    str(out1),
                    "-p",
                    str(out2),
                ],
                check=True,
                capture_output=True,
            )
            return out1.read_text()

        assert retained(lower) == retained(upper), (
            "SHARK has become case-sensitive; the derived references are soft-masked, "
            "so both derivations now need an upper-casing step and new expected digests"
        )
