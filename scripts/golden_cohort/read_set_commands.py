"""Bounded-memory command orchestration for golden-cohort CRAM read evidence."""

from __future__ import annotations

import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FutureTimeoutError
from contextlib import suppress
from pathlib import Path

from golden_cohort.read_sets import ReadSetEvidence, read_name_from_sam_record, summarize_sorted_read_names


def collect_read_set_evidence(
    unmapped_bam: Path,
    samtools_path: str,
    *,
    cwd: Path,
    temporary_parent: Path,
    timeout: int,
) -> ReadSetEvidence:
    """Count and hash one unmapped BAM without retaining its SAM text in memory.

    The record view is consumed on a worker so ``Future.result(timeout=...)`` bounds the
    entire streaming phase, including a child that never closes stdout. QNAMEs alone are
    written to managed temporary storage, externally sorted under the C locale and then
    hashed incrementally.

    Args:
        unmapped_bam: Intermediate BAM produced by a successful non-fast CRAM run.
        samtools_path: Exact samtools executable from the case's complete config.
        cwd: Source tree used as the command working directory.
        temporary_parent: Existing case-log directory for managed temporary files.
        timeout: Wall-clock bound in seconds for each evidence command and the stream.

    Returns:
        Validated count and sorted-QNAME SHA-256 evidence.

    Raises:
        ValueError: If the promised BAM is missing or command outputs disagree.
        RuntimeError: If count, view, sort or streaming fails or times out.
    """
    if not unmapped_bam.is_file():
        raise ValueError(f"unmapped BAM does not exist: {unmapped_bam}")

    completed_count = subprocess.run(
        [samtools_path, "view", "-c", str(unmapped_bam)],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if completed_count.returncode != 0:
        detail = completed_count.stderr.strip() or completed_count.stdout.strip() or "no diagnostic"
        raise RuntimeError(f"samtools read-set command exited {completed_count.returncode}: {detail}")

    with tempfile.TemporaryDirectory(prefix="vntyper-read-set-", dir=temporary_parent) as temporary_dir:
        temporary_root = Path(temporary_dir)
        names_path = temporary_root / "read-names.txt"
        sorted_names_path = temporary_root / "read-names.sorted.txt"
        stderr_path = temporary_root / "samtools-view.stderr"
        view_command = [samtools_path, "view", str(unmapped_bam)]
        with stderr_path.open("w+", encoding="utf-8") as view_stderr:
            view_process = subprocess.Popen(
                view_command,
                cwd=str(cwd),
                stdout=subprocess.PIPE,
                stderr=view_stderr,
                text=True,
            )

            def kill_and_reap() -> None:
                with suppress(ProcessLookupError):
                    view_process.kill()
                view_process.wait()

            view_stdout = view_process.stdout
            if view_stdout is None:
                kill_and_reap()
                raise RuntimeError("samtools read-set command did not expose stdout")

            def consume_view() -> int:
                try:
                    with names_path.open("w", encoding="utf-8") as names_file:
                        for line_number, sam_record in enumerate(view_stdout, start=1):
                            names_file.write(f"{read_name_from_sam_record(sam_record, line_number)}\n")
                finally:
                    view_stdout.close()
                return view_process.wait()

            executor = ThreadPoolExecutor(max_workers=1)
            view_future = executor.submit(consume_view)
            try:
                view_returncode = view_future.result(timeout=timeout)
            except (FutureTimeoutError, subprocess.TimeoutExpired) as exc:
                kill_and_reap()
                raise RuntimeError(f"samtools read-set command timed out after {timeout} seconds") from exc
            except Exception:
                kill_and_reap()
                raise
            finally:
                executor.shutdown(wait=True, cancel_futures=True)

            view_stderr.seek(0)
            view_detail = view_stderr.read().strip() or "no diagnostic"
        if view_returncode != 0:
            raise RuntimeError(f"samtools read-set command exited {view_returncode}: {view_detail}")

        sort_env = dict(os.environ)
        sort_env["LC_ALL"] = "C"
        completed_sort = subprocess.run(
            ["sort", "-o", str(sorted_names_path), str(names_path)],
            cwd=str(cwd),
            env=sort_env,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        if completed_sort.returncode != 0:
            detail = completed_sort.stderr.strip() or completed_sort.stdout.strip() or "no diagnostic"
            raise RuntimeError(f"read-name sort command exited {completed_sort.returncode}: {detail}")
        with sorted_names_path.open(encoding="utf-8") as sorted_names:
            return summarize_sorted_read_names(completed_count.stdout, sorted_names)
