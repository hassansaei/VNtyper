"""Generate and extract the local IGV report used by the summary renderer."""

from __future__ import annotations

import logging
import os
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from vntyper.scripts import report_assets
from vntyper.scripts.report_formatting import extract_igv_fragments

logger = logging.getLogger(__name__)


def run_igv_report(
    bed_file: str | Path,
    bam_file: str | Path | None,
    fasta_file: str | Path | None,
    output_html: str | Path,
    flanking: int = 50,
    vcf_file: str | Path | None = None,
    config: Mapping[str, Any] | None = None,
    *,
    report_igv: str,
) -> None:
    """Run igv-reports against VNtyper's controlled, network-free template.

    igv-reports' ``--standalone`` mode follows the CDN URL in its own template at
    generation time. VNtyper instead passes a shipped local template with no external
    resources. In sidecar mode, the verified vendored igv.js source is inserted only
    after igv-reports has populated the table and session literals. Embedded mode uses
    the small generated page only as a temporary extraction artifact.

    Args:
        bed_file: BED file defining the displayed loci.
        bam_file: Optional BAM/CRAM alignment track.
        fasta_file: Reference FASTA used to embed sequence data.
        output_html: Destination for igv-reports' generated page.
        flanking: Bases to include around each locus.
        vcf_file: Optional VCF track.
        config: Optional configuration carrying the default flanking value.
        report_igv: ``embedded`` for a temporary extraction page or ``sidecar`` for
            the final, self-contained alignment browser.

    Raises:
        ValueError: If the mode cannot produce a report, the controlled template is
            missing, or its library marker is not present exactly once in a sidecar.
        subprocess.CalledProcessError: If igv-reports fails.
    """
    if report_igv not in (report_assets.REPORT_IGV_EMBEDDED, report_assets.REPORT_IGV_SIDECAR):
        msg = f"--report-igv {report_igv!r} cannot generate an IGV report."
        logger.error(msg)
        raise ValueError(msg)
    if not report_assets.IGV_REPORT_TEMPLATE_PATH.is_file():
        msg = f"Controlled IGV report template is missing: {report_assets.IGV_REPORT_TEMPLATE_PATH}."
        logger.error(msg)
        raise ValueError(msg)

    logger.debug(
        "run_igv_report called with bed_file=%s, bam_file=%s, fasta_file=%s, output_html=%s, "
        "vcf_file=%s, flanking=%s, report_igv=%s",
        bed_file,
        bam_file,
        fasta_file,
        output_html,
        vcf_file,
        flanking,
        report_igv,
    )

    if config is not None and flanking == 50:
        flanking = config.get("default_values", {}).get("flanking", 50)
        logger.debug("Flanking region set to %s based on config.", flanking)

    command = [
        "create_report",
        str(bed_file),
        "--template",
        str(report_assets.IGV_REPORT_TEMPLATE_PATH),
        "--flanking",
        str(flanking),
        "--fasta",
        str(fasta_file) if fasta_file else "",
    ]
    tracks: list[str] = []
    for track_name, track_path in (("VCF", vcf_file), ("BAM", bam_file)):
        if track_path and os.path.exists(track_path):
            tracks.append(str(track_path))
            logger.debug("Adding %s track: %s", track_name, track_path)
        elif track_path:
            logger.warning(
                "%s file specified but not found: %s. Skipping %s track.", track_name, track_path, track_name
            )
        else:
            logger.info("No %s file provided. IGV report will not include %s track.", track_name, track_name)
    if tracks:
        command.extend(["--tracks", *tracks])
    else:
        logger.warning("No valid tracks are available; the IGV report will show embedded reference sequence only.")
    command.extend(["--output", str(output_html)])

    logger.info("Running IGV report: %s", " ".join(command))
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError:
        logger.error("Error generating IGV report: %s", " ".join(command))
        raise
    except Exception as error:
        logger.error("Unexpected error generating IGV report: %s", error)
        raise

    if report_igv == report_assets.REPORT_IGV_SIDECAR:
        output_path = Path(output_html)
        generated = output_path.read_text(encoding="utf-8")
        marker_count = generated.count(report_assets.IGV_REPORT_LIBRARY_MARKER)
        if marker_count != 1:
            msg = f"Generated IGV sidecar does not contain exactly one verified library marker; found {marker_count}."
            logger.error(msg)
            raise ValueError(msg)
        source = report_assets.igv_library_source().decode("utf-8")
        output_path.write_text(generated.replace(report_assets.IGV_REPORT_LIBRARY_MARKER, source), encoding="utf-8")
    logger.info("IGV report successfully generated at %s", output_html)


def extract_igv_content(igv_report_html: str | Path) -> tuple[str, str, str]:
    """Extract container markup and the two data literals from an IGV report.

    Args:
        igv_report_html: Page written by :func:`run_igv_report`.

    Returns:
        tuple[str, str, str]: Container markup, table JSON and session dictionary,
        or three empty strings when the optional artifact cannot be read.
    """
    logger.debug("extract_igv_content called with igv_report_html=%s", igv_report_html)
    try:
        content = Path(igv_report_html).read_text(encoding="utf-8")
    except FileNotFoundError:
        logger.error("IGV report file not found: %s", igv_report_html)
        return "", "", ""
    except Exception as error:
        logger.error("Unexpected error reading IGV report: %s", error)
        return "", "", ""
    return extract_igv_fragments(content)
