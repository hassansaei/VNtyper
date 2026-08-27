"""Shared helpers for ordinary pipeline-job handoff tests."""

import hashlib
from pathlib import Path


def _identity(path: Path) -> list[int]:
    metadata = path.stat()
    return [metadata.st_dev, metadata.st_ino]


def _digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _invoke_job(tasks, **kwargs) -> None:
    identity = kwargs.get("workspace_identity")
    bam_path = Path(kwargs["bam_path"])
    index_path_value = kwargs.get("index_path")
    if isinstance(identity, dict) and "alignment_sha256" not in identity:
        identity["alignment_sha256"] = _digest(bam_path.read_bytes())
        identity["index_sha256"] = None if index_path_value is None else _digest(Path(index_path_value).read_bytes())
    tasks.settings.DEFAULT_HANDOFF_SPOOL_DIR = str(bam_path.parent.parent)
    tasks.settings.DEFAULT_INPUT_DIR = str(bam_path.parent.parent.with_name("shared-input"))
    tasks.settings.DEFAULT_OUTPUT_DIR = str(Path(kwargs["output_dir"]).parent)
    Path(tasks.settings.DEFAULT_INPUT_DIR).mkdir(exist_ok=True)
    tasks.run_vntyper_job.push_request(id="task-inode")
    try:
        tasks.run_vntyper_job.run(**kwargs)
    finally:
        tasks.run_vntyper_job.pop_request()
