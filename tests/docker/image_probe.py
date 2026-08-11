"""
In-container probe for the VNtyper image.

Executed *inside* the image by tests/docker/test_image_structure.py, which pipes this
file to the image's own interpreter and parses the single JSON document it writes to
stdout. It runs under Python 3.9 with the standard library only: it must not import
pytest, and it must not be collected by the host test run (hence no `test_` prefix).

One container start covers every runtime assertion. Twenty `docker run` invocations
would each re-pay interpreter startup on a 4.8 GB image; this pays it once and returns
data that pytest turns into individually-named tests.

Resolving config.json the way vntyper does - from the INSTALLED package, then relative
to WORKDIR - is the whole point of this file. A hardcoded list of reference paths would
drift from config.json, and that drift is precisely the bug this tier exists to catch.
"""

import json
import os
import shutil
import subprocess
import sys

PREFIX = "/opt/vntyper"
CONDA_ENVS = ("vntyper", "envadvntr", "shark_env")
REQUIRED_BINARIES = ("bwa", "samtools", "fastp", "bcftools", "java")
# The six physical reference ids install_references_config.json declares BWA entries
# for (ucsc_references + ncbi_references + ensembl_references). config.json ships all
# six bwa_reference_* keys as real paths, so the sidecar check below must cover all six
# too, not just the two UCSC ones - otherwise a missing GRCh37/GRCh38/*_ensembl index
# would slip past this tier and only surface as a family-fallback warning at runtime.
BWA_PHYSICAL_IDS = ("hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl")
# `git` is deliberately absent: condaforge/miniforge3 ships it in its own layer, so
# asserting on it would be a permanently red test we cannot fix from this Dockerfile.
BUILD_TOOLS = ("gcc", "g++", "cc", "make", "ld")


def _run(argv):
    """Run argv and return a JSON-safe record rather than raising.

    Args:
        argv (list): Command to execute.

    Returns:
        dict: Keys ``rc``, ``out`` and ``err``.
    """
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, timeout=60, check=False)
        return {
            "rc": proc.returncode,
            "out": proc.stdout.strip(),
            "err": proc.stderr.strip()[:300],
        }
    except Exception as exc:
        return {"rc": 127, "out": "", "err": repr(exc)}


def _entry(rel):
    """Resolve a config-declared path the way the pipeline does, and stat it.

    Args:
        rel (str): Path as written in a shipped config file.

    Returns:
        dict: Declared path, resolved path, existence and size.
    """
    path = rel if os.path.isabs(rel) else os.path.join(PREFIX, rel)
    exists = os.path.isfile(path)
    return {
        "declared": rel,
        "resolved": path,
        "exists": exists,
        "size": os.path.getsize(path) if exists else 0,
    }


def _collect_references(pkg_dir, config):
    """Collect every reference path the shipped configs declare, plus BWA sidecars.

    Args:
        pkg_dir (str): Directory of the installed ``vntyper`` package.
        config (dict): Parsed ``config.json``.

    Returns:
        dict: Mapping of config key to a path record.
    """
    refs = {}
    for key, value in sorted(config.get("reference_data", {}).items()):
        if value is None:
            continue
        if not isinstance(value, str):
            raise ValueError(f"reference_data.{key} must be a path string or null")
        refs[key] = _entry(value)

    shark_config = os.path.join(pkg_dir, "modules", "shark", "shark_config.json")
    if os.path.isfile(shark_config):
        with open(shark_config) as handle:
            for key, rel in json.load(handle).get("shark_settings", {}).items():
                if isinstance(rel, str) and rel.endswith((".fa", ".fasta")):
                    refs["shark." + key] = _entry(rel)

    # bwa_index_extensions is itself config-declared, so the sidecar set stays in sync
    # with config.json rather than with a list maintained here.
    extensions = config.get("tool_params", {}).get("bwa_index_extensions", [])
    for assembly in BWA_PHYSICAL_IDS:
        base = config.get("reference_data", {}).get("bwa_reference_" + assembly)
        if base:
            for ext in extensions:
                refs[f"bwa_index_{assembly}{ext}"] = _entry(base + ext)

    return refs


def main():
    """Emit a single JSON report describing the image's structure."""
    import importlib.metadata as metadata

    import vntyper

    pkg_dir = os.path.dirname(vntyper.__file__)
    with open(os.path.join(pkg_dir, "config.json")) as handle:
        config = json.load(handle)

    json.dump(
        {
            "uid": os.getuid(),
            "gid": os.getgid(),
            "package_version": metadata.version("vntyper"),
            "cli_version": _run(["vntyper", "--version"]),
            "conda_env_dirs": {env: os.path.isdir(f"/opt/conda/envs/{env}") for env in CONDA_ENVS},
            # shark_env holds only bin/shark and has no interpreter, so only the two
            # Python environments are probed for a version.
            "interpreters": {
                env: _run(
                    [
                        f"/opt/conda/envs/{env}/bin/python",
                        "-c",
                        "import sys;print('%d.%d' % sys.version_info[:2])",
                    ]
                )
                for env in ("vntyper", "envadvntr")
            },
            "shark_binary": os.path.isfile("/opt/conda/envs/shark_env/bin/shark"),
            "advntr_import": _run(
                [
                    "/opt/conda/envs/envadvntr/bin/python",
                    "-c",
                    "import advntr; print(advntr.__file__)",
                ]
            ),
            "binaries": {name: shutil.which(name) for name in REQUIRED_BINARIES},
            "leaked_build_tools": {name: shutil.which(name) for name in BUILD_TOOLS if shutil.which(name)},
            "references": _collect_references(pkg_dir, config),
            "kestrel_jar": bool(
                [name for name in os.listdir(os.path.join(pkg_dir, "dependencies", "kestrel")) if name.endswith(".jar")]
            )
            if os.path.isdir(os.path.join(pkg_dir, "dependencies", "kestrel"))
            else False,
        },
        sys.stdout,
    )


if __name__ == "__main__":
    main()
