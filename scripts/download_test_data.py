#!/usr/bin/env python3
"""
Standalone script to download and validate VNtyper test data.

This script downloads the test data archive from Zenodo and extracts it to the
correct location. It uses the same extraction logic as the test utilities to ensure
consistent behavior between local development and CI environments.

Usage:
    python scripts/download_test_data.py [--force] [--verify-only]

Options:
    --force         Force re-download even if files exist
    --verify-only   Only verify existing files, don't download
    --quiet         Minimal output
    --verbose       Detailed output

Environment Variables:
    VNTYPER_TEST_DATA_DIR: Override test data directory (default: tests/data)

This script is designed to be run:
1. Locally by developers: `make download-test-data` or `python scripts/download_test_data.py`
2. In CI workflows: Called before running tests to ensure data is present
"""

import argparse
import hashlib
import json
import logging
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)8s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def compute_md5(file_path: Path) -> str:
    """Compute MD5 hash of a file."""
    hasher = hashlib.md5()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def download_file_requests(url: str, dest_path: Path, timeout: int = 300) -> None:
    """
    Download a file using Python requests library.

    Args:
        url: Remote URL to download from
        dest_path: Local path to save file
        timeout: Timeout in seconds for connection and read operations
    """
    try:
        import requests
    except ImportError:
        logger.error("requests library not installed. Install with: pip install requests")
        sys.exit(1)

    logger.info(f"Downloading {url}")
    logger.info(f"Destination: {dest_path}")

    resp = requests.get(url, stream=True, timeout=(timeout, timeout))
    resp.raise_for_status()

    total_size = int(resp.headers.get("content-length", 0))
    if total_size > 0:
        logger.info(f"File size: {total_size / (1024 * 1024):.2f} MB")

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    chunk_size = 65536  # 64KB chunks

    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)

                # Log progress every 50MB for large files
                if total_size > 10 * 1024 * 1024 and downloaded % (50 * 1024 * 1024) < chunk_size:
                    progress = (downloaded / total_size * 100) if total_size > 0 else 0
                    logger.info(
                        f"Download progress: {downloaded / (1024 * 1024):.1f} MB / "
                        f"{total_size / (1024 * 1024):.1f} MB ({progress:.1f}%)"
                    )

    logger.info(f"Download complete: {downloaded / (1024 * 1024):.2f} MB")

    # Verify download completed successfully
    actual_size = dest_path.stat().st_size
    if total_size > 0 and actual_size != total_size:
        logger.error(f"Download size mismatch: expected {total_size} bytes, got {actual_size} bytes")
        raise RuntimeError(f"Incomplete download: {actual_size}/{total_size} bytes")

    logger.info(f"✓ Download verified: {actual_size / (1024 * 1024):.2f} MB")


def extract_archive(archive_path: Path, extract_to: Path) -> None:
    """
    Extract zip archive with intelligent handling of nested directories.

    This uses the SAME logic as tests/test_data_utils.py to ensure
    consistent extraction behavior.

    Args:
        archive_path: Path to zip archive
        extract_to: Directory to extract files to
    """
    logger.info(f"Extracting archive to {extract_to}")
    extract_to.mkdir(parents=True, exist_ok=True)

    # Verify zip file integrity first
    logger.info("Verifying archive integrity...")
    try:
        with zipfile.ZipFile(archive_path, "r") as zip_test:
            bad_file = zip_test.testzip()
            if bad_file:
                raise RuntimeError(f"Corrupted file in archive: {bad_file}")
        logger.info("✓ Archive integrity verified")
    except zipfile.BadZipFile as e:
        logger.error(f"Archive is corrupted: {e}")
        raise

    with zipfile.ZipFile(archive_path, "r") as zip_ref:
        all_files = zip_ref.namelist()
        logger.info(f"Archive contains {len(all_files)} entries")

        # List first few files for debugging
        logger.info("First 10 entries in archive:")
        for i, name in enumerate(all_files[:10]):
            logger.info(f"  [{i+1}] {name}")
        if len(all_files) > 10:
            logger.info(f"  ... and {len(all_files) - 10} more entries")

        # Robust detection of dominant top-level directory
        dir_counts = {}
        files_at_root = 0

        for name in all_files:
            if name.endswith("/"):  # Skip directory entries
                continue

            if "/" in name:
                top_dir = name.split("/", 1)[0] + "/"
                dir_counts[top_dir] = dir_counts.get(top_dir, 0) + 1
            else:
                files_at_root += 1

        # Determine if there's a dominant top-level directory to strip
        common_prefix = ""
        if dir_counts:
            dominant_dir, file_count = max(dir_counts.items(), key=lambda x: x[1])
            total_files = sum(dir_counts.values()) + files_at_root

            logger.info(f"Top-level directories: {dict(dir_counts)}")
            logger.info(f"Files at root: {files_at_root}")
            logger.info(f"Dominant directory: '{dominant_dir}' with {file_count}/{total_files} files")

            # If dominant directory contains >50% of files, strip it
            # Lowered from 0.9 to 0.5 to handle archives with extra root files (README, etc)
            if file_count / total_files > 0.5:
                common_prefix = dominant_dir
                logger.info(f"Will strip '{common_prefix}' from extraction paths (ratio: {file_count}/{total_files})")

        if common_prefix:
            logger.info(f"Extracting files while stripping '{common_prefix}' prefix...")
            logger.info(f"Target directory: {extract_to.absolute()}")
            extracted_count = 0
            total_members = len([m for m in zip_ref.infolist() if not m.filename.endswith("/")])

            for member in zip_ref.infolist():
                if member.filename.endswith("/"):  # Skip directory entries
                    continue

                member_path = member.filename
                original_path = member_path  # Keep for logging

                if member_path.startswith(common_prefix):
                    member_path = member_path[len(common_prefix) :]

                # Skip files not in dominant directory
                if not member_path or member_path == member.filename:
                    logger.debug(f"Skipping: {member.filename} (not in dominant directory)")
                    continue

                target_path = extract_to / member_path
                target_path.parent.mkdir(parents=True, exist_ok=True)

                # Log progress for large extractions
                extracted_count += 1
                if extracted_count % 10 == 0 or extracted_count == total_members:
                    logger.info(f"Extracting [{extracted_count}/{total_members}]: {original_path} → {target_path}")

                # Use shutil.copyfileobj for robust large file extraction
                try:
                    with zip_ref.open(member) as source, open(target_path, "wb") as target:
                        shutil.copyfileobj(source, target, length=65536)

                    # Verify file was written
                    if not target_path.exists():
                        raise RuntimeError(f"File not found after extraction: {target_path}")

                    file_size = target_path.stat().st_size
                    if file_size == 0:
                        raise RuntimeError(f"Extracted file is empty: {target_path}")

                    # Log details for first few files
                    if extracted_count <= 3:
                        logger.info(f"  ✓ Verified: {target_path.name} ({file_size / (1024*1024):.2f} MB)")

                except Exception as e:
                    logger.error(f"Failed to extract {member.filename} to {target_path}: {e}")
                    raise

            logger.info(f"✓ Successfully extracted {extracted_count} files to {extract_to.absolute()}")
        else:
            logger.info("Extracting all files normally...")
            zip_ref.extractall(extract_to)

    logger.info("Archive extracted successfully")


def verify_test_data(config_path: Path, data_dir: Path) -> tuple[bool, list[str]]:
    """
    Verify test data files exist and have correct MD5 checksums.

    Args:
        config_path: Path to test_data_config.json
        data_dir: Base directory for test data

    Returns:
        Tuple of (all_valid: bool, error_messages: list[str])
    """
    with config_path.open("r") as f:
        config = json.load(f)

    file_resources = config.get("file_resources", [])
    errors = []
    verified_count = 0

    logger.info(f"Verifying {len(file_resources)} test data files...")

    for resource in file_resources:
        # Compute path relative to data_dir
        local_dir = Path(resource["local_path"])
        filename = resource["filename"]

        # Handle both absolute and relative paths
        if local_dir.is_absolute():
            local_path = local_dir / filename
        else:
            # If local_path starts with "tests/data", use data_dir as base
            local_path_str = str(local_dir / filename)
            if local_path_str.startswith("tests/data/"):
                local_path = data_dir / local_path_str.replace("tests/data/", "")
            else:
                local_path = data_dir / filename

        expected_md5 = resource.get("md5sum")

        if not local_path.exists():
            errors.append(f"Missing: {local_path}")
            continue

        if expected_md5:
            current_md5 = compute_md5(local_path)
            if current_md5.lower() != expected_md5.lower():
                errors.append(f"MD5 mismatch: {local_path} (expected={expected_md5}, got={current_md5})")
                continue

        verified_count += 1
        logger.debug(f"✓ Verified: {local_path}")

    if errors:
        logger.error(f"Verification failed: {len(errors)} issues found")
        for err in errors[:10]:  # Show first 10 errors
            logger.error(f"  - {err}")
        if len(errors) > 10:
            logger.error(f"  ... and {len(errors) - 10} more")
        return False, errors
    else:
        logger.info(f"✓ All {verified_count} test data files verified successfully")
        return True, []


def main():
    parser = argparse.ArgumentParser(
        description="Download and verify VNtyper test data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Download test data:
    python scripts/download_test_data.py

  Force re-download:
    python scripts/download_test_data.py --force

  Only verify existing data:
    python scripts/download_test_data.py --verify-only

  Quiet mode for CI:
    python scripts/download_test_data.py --quiet
        """,
    )
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    parser.add_argument("--verify-only", action="store_true", help="Only verify existing files, don't download")
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    parser.add_argument("--verbose", action="store_true", help="Detailed output")

    args = parser.parse_args()

    # Configure logging level
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    elif args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Determine paths
    repo_root = Path(__file__).parent.parent
    config_path = repo_root / "tests" / "test_data_config.json"
    data_dir = Path(repo_root / "tests" / "data")

    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        sys.exit(1)

    # Load config
    with config_path.open("r") as f:
        config = json.load(f)

    archive_config = config.get("archive_file")
    if not archive_config:
        logger.error("No archive_file configuration found in test_data_config.json")
        sys.exit(1)

    # Verify only mode
    if args.verify_only:
        logger.info("Verify-only mode: checking existing files")
        success, _ = verify_test_data(config_path, data_dir)
        sys.exit(0 if success else 1)

    # Check if download needed
    need_download = args.force
    if not need_download:
        logger.info("Checking if test data exists and is valid...")
        success, _ = verify_test_data(config_path, data_dir)
        if success:
            logger.info("Test data already exists and is valid. Use --force to re-download.")
            sys.exit(0)
        logger.info("Test data is missing or invalid. Starting download...")
        need_download = True

    # Download archive
    archive_url = archive_config["url"]
    extract_to = data_dir

    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)

    try:
        logger.info("=" * 80)
        logger.info("DOWNLOADING TEST DATA")
        logger.info("=" * 80)
        logger.info(f"Source: {archive_url}")
        logger.info("This may take 10-30 minutes depending on network speed...")
        logger.info("")

        download_file_requests(archive_url, tmp_path)

        logger.info("")
        logger.info("=" * 80)
        logger.info("EXTRACTING ARCHIVE")
        logger.info("=" * 80)

        extract_archive(tmp_path, extract_to)

        logger.info("")
        logger.info("=" * 80)
        logger.info("VERIFYING DATA")
        logger.info("=" * 80)

        success, errors = verify_test_data(config_path, data_dir)

        if not success:
            logger.error("Verification failed after extraction!")
            sys.exit(1)

        logger.info("")
        logger.info("=" * 80)
        logger.info("✓ TEST DATA DOWNLOAD COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Test data installed to: {data_dir}")
        logger.info("You can now run tests with: make test")

    except Exception as e:
        logger.error(f"Download failed: {e}")
        sys.exit(1)
    finally:
        # Cleanup temp file
        if tmp_path.exists():
            tmp_path.unlink()


if __name__ == "__main__":
    main()
