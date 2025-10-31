"""
Shared test data utilities for downloading and validating test data.

This module provides functions to download test data from Zenodo and validate
MD5 checksums. It's used by both integration and Docker tests.

Environment Variables:
    VNTYPER_TEST_DATA_SKIP_DOWNLOAD: If set to "1" or "true", skip automatic downloads
        and fail fast if data is missing. Used in CI environments where data should
        be pre-downloaded by the CI workflow.
"""

import hashlib
import logging
import os
import tempfile
import zipfile
from pathlib import Path

import pytest
import requests

logger = logging.getLogger(__name__)


def compute_md5(file_path: Path) -> str:
    """
    Compute the MD5 hash of a file, returning the hex digest string.

    Args:
        file_path: Path to the file whose MD5 should be computed.

    Returns:
        Hex digest string of the file's MD5 sum.
    """
    hasher = hashlib.md5()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def download_file(url: str, dest_path: Path, timeout: int = 60) -> None:
    """
    Download a file from URL and save it to dest_path.

    Args:
        url: The remote URL to download from.
        dest_path: The local path where the downloaded file is saved.
        timeout: Timeout in seconds for the download request.

    Raises:
        requests.HTTPError: If the download fails (non-200 status code).
    """
    logger.info(f"Downloading {url} to {dest_path}")
    # Use tuple timeout: (connection timeout, read timeout)
    # This ensures the request doesn't hang indefinitely
    resp = requests.get(url, stream=True, timeout=(timeout, timeout))
    resp.raise_for_status()

    # Get file size if available for progress reporting
    total_size = int(resp.headers.get("content-length", 0))
    if total_size > 0:
        logger.info(f"File size: {total_size / (1024 * 1024):.2f} MB")

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    chunk_size = 65536  # 64KB chunks

    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            if chunk:  # filter out keep-alive chunks
                f.write(chunk)
                downloaded += len(chunk)

                # Log progress every 10MB for large files
                if total_size > 10 * 1024 * 1024 and downloaded % (10 * 1024 * 1024) < chunk_size:
                    progress = (downloaded / total_size * 100) if total_size > 0 else 0
                    logger.info(
                        f"Download progress: {downloaded / (1024 * 1024):.1f} MB / {total_size / (1024 * 1024):.1f} MB ({progress:.1f}%)"
                    )

    logger.info(f"Download complete: {dest_path} ({downloaded / (1024 * 1024):.2f} MB)")


def ensure_test_data_downloaded(test_config: dict) -> None:
    """
    Ensure all test data is present and valid (by MD5).

    This function supports both:
    1. Archive-based downloads (data.zip from Zenodo)
    2. Individual file downloads (legacy support)

    If 'archive_file' is specified in config, it will download and extract the archive
    when any file is missing or has MD5 mismatch. Otherwise, it falls back to individual
    file downloads.

    In CI environments (when VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1), this function will
    fail immediately if data is missing instead of attempting to download. This ensures
    fast failure and prevents pytest fixtures from consuming test execution time.

    Args:
        test_config: The loaded JSON config with:
            - archive_file (optional): dict with url, extract_to
            - file_resources: list of file dicts with local_path, filename, md5sum
            - server_base_url (optional): base URL for individual file downloads

    Raises:
        SystemExit: If download fails or files don't match expected MD5 sums.
            In CI mode, exits immediately if data is missing without attempting download.
    """
    file_resources = test_config.get("file_resources", [])
    archive_config = test_config.get("archive_file")

    # Check if running in CI mode (download disabled)
    skip_download = os.environ.get("VNTYPER_TEST_DATA_SKIP_DOWNLOAD", "").lower() in ("1", "true", "yes")

    if skip_download:
        logger.info("Running in CI mode - automatic downloads disabled (VNTYPER_TEST_DATA_SKIP_DOWNLOAD is set)")

    # Check if we need to download anything
    need_download = False
    missing_files = []
    mismatched_files = []

    for resource in file_resources:
        local_dir = Path(resource["local_path"])
        filename = resource["filename"]
        local_path = local_dir / filename
        expected_md5 = resource.get("md5sum")

        # Skip if no MD5 check required
        if not expected_md5:
            continue

        if not local_path.exists():
            logger.info(f"File {local_path} is missing.")
            missing_files.append(str(local_path))
            need_download = True
            continue

        current_md5 = compute_md5(local_path)
        if current_md5.lower() != expected_md5.lower():
            logger.warning(f"File {local_path} has MD5 mismatch. Expected={expected_md5}, Got={current_md5}")
            mismatched_files.append(f"{local_path} (MD5: expected={expected_md5}, got={current_md5})")
            need_download = True
            continue

    if not need_download:
        logger.info("All test data files verified. No download needed.")
        return

    # If in CI mode and data is missing, fail immediately with helpful message
    if skip_download:
        error_msg = ["Test data is missing or invalid, but automatic downloads are disabled in CI mode."]
        error_msg.append("\nThis indicates a problem with the CI cache or download step.")
        error_msg.append("\n\nTo fix this:")
        error_msg.append("  1. Check that the cache key in .github/workflows/*.yml is correct")
        error_msg.append("  2. Ensure the 'Download test data' step completed successfully")
        error_msg.append("  3. Verify the extraction path matches what tests expect (tests/data)")
        error_msg.append("\n\nFor local development, download test data by running:")
        error_msg.append("  make download-test-data")

        if missing_files:
            error_msg.append(f"\n\nMissing files ({len(missing_files)}):")
            for f in missing_files[:10]:  # Show first 10
                error_msg.append(f"  - {f}")
            if len(missing_files) > 10:
                error_msg.append(f"  ... and {len(missing_files) - 10} more")

        if mismatched_files:
            error_msg.append(f"\n\nFiles with MD5 mismatch ({len(mismatched_files)}):")
            for f in mismatched_files[:10]:  # Show first 10
                error_msg.append(f"  - {f}")
            if len(mismatched_files) > 10:
                error_msg.append(f"  ... and {len(mismatched_files) - 10} more")

        pytest.exit("\n".join(error_msg), returncode=1)

    # Need to download - use archive if configured
    if archive_config:
        archive_url = archive_config["url"]
        extract_to = Path(archive_config["extract_to"])

        logger.info(f"Downloading test data archive from {archive_url}")

        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)

        try:
            # Download archive
            download_file(archive_url, tmp_path, timeout=300)  # 5 min timeout for large archive
            logger.info(f"Archive downloaded to {tmp_path}")

            # Extract archive
            logger.info(f"Extracting archive to {extract_to}")
            extract_to.mkdir(parents=True, exist_ok=True)

            with zipfile.ZipFile(tmp_path, "r") as zip_ref:
                # Get all file names in the archive
                all_files = zip_ref.namelist()
                logger.info(f"Archive contains {len(all_files)} entries")

                # Robust detection of dominant top-level directory
                # This handles cases where the zip has a nested "data/" directory
                # even if there are other files at root (like README.md)
                dir_counts: dict[str, int] = {}
                files_at_root = 0

                for name in all_files:
                    # Skip directory entries (they end with /)
                    if name.endswith("/"):
                        continue

                    # Check if file is in a subdirectory or at root
                    if "/" in name:
                        top_dir = name.split("/", 1)[0] + "/"
                        dir_counts[top_dir] = dir_counts.get(top_dir, 0) + 1
                    else:
                        files_at_root += 1

                # Determine if there's a dominant top-level directory to strip
                common_prefix = ""
                if dir_counts:
                    # Find the directory with the most files
                    dominant_dir, file_count = max(dir_counts.items(), key=lambda x: x[1])
                    total_files = sum(dir_counts.values()) + files_at_root

                    logger.info(f"Top-level directories: {dict(dir_counts)}")
                    logger.info(f"Files at root: {files_at_root}")
                    logger.info(f"Dominant directory: '{dominant_dir}' with {file_count}/{total_files} files")

                    # Only strip prefix if:
                    # 1. Dominant directory has >90% of files (very dominant)
                    # 2. OR: Has >80% AND fewer than 5 files at root (just metadata files like README)
                    #
                    # This handles two cases:
                    # - Case 1: data/ with 28 files + README.md at root (96.5% → strip prefix)
                    # - Case 2: remapped/ with 87 files + 36 test files at root (70.7% → DON'T strip)
                    ratio = file_count / total_files
                    if ratio > 0.9 or (ratio > 0.8 and files_at_root < 5):
                        common_prefix = dominant_dir
                        logger.info(f"Will strip '{common_prefix}' from extraction paths (ratio: {ratio:.1%}, root files: {files_at_root})")
                    else:
                        logger.info(f"Mixed archive structure detected (ratio: {ratio:.1%}, root files: {files_at_root})")
                        logger.info("Will extract all files normally (no prefix stripping)")

                if common_prefix:
                    logger.info("Extracting files while stripping root directory...")

                    # Extract each file, stripping the common prefix
                    for member in zip_ref.infolist():
                        # Skip directory entries
                        if member.filename.endswith("/"):
                            continue

                        # Remove common prefix from the path
                        member_path = member.filename
                        if member_path.startswith(common_prefix):
                            member_path = member_path[len(common_prefix) :]

                        # Skip files that don't have the prefix (e.g., README.md at root)
                        if not member_path or member_path == member.filename:
                            logger.info(f"Skipping file not in dominant directory: {member.filename}")
                            continue

                        # Extract to target directory
                        target_path = extract_to / member_path
                        target_path.parent.mkdir(parents=True, exist_ok=True)

                        with zip_ref.open(member) as source, open(target_path, "wb") as target:
                            target.write(source.read())
                else:
                    # No dominant directory, extract normally
                    logger.info("No dominant root directory detected, extracting all files normally...")
                    zip_ref.extractall(extract_to)

            logger.info("Archive extracted successfully")

            # Verify all files after extraction
            for resource in file_resources:
                local_dir = Path(resource["local_path"])
                filename = resource["filename"]
                local_path = local_dir / filename
                expected_md5 = resource.get("md5sum")

                if not local_path.exists():
                    pytest.exit(
                        f"File {local_path} not found after archive extraction!",
                        returncode=1,
                    )

                if expected_md5:
                    current_md5 = compute_md5(local_path)
                    if current_md5.lower() != expected_md5.lower():
                        pytest.exit(
                            f"File {local_path} MD5 mismatch after extraction.\n"
                            f"Expected={expected_md5}, Got={current_md5}",
                            returncode=1,
                        )
                logger.info(f"Verified {local_path}")
        finally:
            # Clean up temp file
            if tmp_path.exists():
                tmp_path.unlink()
    else:
        # Legacy: download individual files
        logger.info("Using legacy individual file download mode")
        server_base_url = test_config.get("server_base_url", "")

        for resource in file_resources:
            local_dir = Path(resource["local_path"])
            filename = resource["filename"]
            local_path = local_dir / filename
            expected_md5 = resource.get("md5sum")

            # Skip if file exists and MD5 matches
            if local_path.exists() and expected_md5:
                current_md5 = compute_md5(local_path)
                if current_md5.lower() == expected_md5.lower():
                    logger.info(f"File {local_path} already exists and is valid")
                    continue

            # Download file
            url_suffix = resource.get("url_suffix", "")
            if not url_suffix:
                logger.warning(f"No url_suffix for {filename}, skipping download")
                continue

            file_url = server_base_url + url_suffix
            logger.info(f"Downloading {filename} from {file_url}")

            try:
                download_file(file_url, local_path)

                # Verify MD5 after download
                if expected_md5:
                    current_md5 = compute_md5(local_path)
                    if current_md5.lower() != expected_md5.lower():
                        pytest.exit(
                            f"Downloaded file {local_path} MD5 mismatch.\nExpected={expected_md5}, Got={current_md5}",
                            returncode=1,
                        )
                logger.info(f"Successfully downloaded and verified {local_path}")
            except Exception as e:
                pytest.exit(f"Failed to download {filename}: {e}", returncode=1)
