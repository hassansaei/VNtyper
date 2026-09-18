#!/usr/bin/env python3
"""Generate deterministic calibration inputs and verify their external commitment.

These commands establish input identity only. They do not run callers, evaluate
performance, or authorize a calibrated model. Keep generated studies outside Git.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from calibration_sim.generation import generate_simulation_bundle, verify_generated_bundle
from calibration_sim.protocol import decode_simulation_protocol

from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Run generation or verification without interpreting inputs as results.

    Args:
        argv: Explicit argument list, or the process arguments when omitted.

    Returns:
        Zero for a completed input operation, one for a refused or failed operation.

    Raises:
        SystemExit: Argparse's exit two for invalid command usage.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    generate = commands.add_parser("generate", help="generate private deterministic input files")
    generate.add_argument("--protocol", type=Path, required=True)
    generate.add_argument("--output", type=Path, required=True)
    verify = commands.add_parser("verify-inputs", help="verify exact input inventory against an external digest")
    verify.add_argument("--input", type=Path, required=True)
    verify.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args(argv)
    try:
        if args.operation == "generate":
            protocol = decode_simulation_protocol(load_strict_json_object(read_regular_path(args.protocol)))
            generated = generate_simulation_bundle(protocol, args.output)
            manifest_sha256 = generated.manifest_sha256
            status = "generated-inputs-only"
        else:
            protocol = verify_generated_bundle(args.input, expected_manifest_sha256=args.manifest_sha256)
            manifest_sha256 = args.manifest_sha256
            status = "verified-inputs-only"
    except (OSError, ValueError, RuntimeError) as error:
        logger.error(f"calibration simulation input operation failed: {error}")
        return 1
    print(
        json.dumps(
            {
                "schema_version": "calibration-simulation-input-operation-v1",
                "evidence_status": status,
                "protocol_sha256": protocol.sha256,
                "manifest_sha256": manifest_sha256,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
