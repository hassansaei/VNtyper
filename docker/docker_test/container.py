"""
Docker container operations with robust error handling.

Implements research-based best practices for subprocess Docker interactions:
- subprocess.run() with timeout parameter
- Comprehensive exception handling (TimeoutExpired, CalledProcessError, FileNotFoundError)
- BuildKit integration for faster builds
- Clear, actionable error messages
"""

import os
import subprocess
from pathlib import Path
from typing import Optional

from .utils.logging import log


# Custom exception classes for clear error handling
class DockerError(Exception):
    """Base exception for Docker operations."""

    pass


class DockerNotFoundError(DockerError):
    """Docker command not found on system."""

    pass


class DockerBuildError(DockerError):
    """Docker build failed."""

    pass


class DockerExecutionError(DockerError):
    """Docker run command failed."""

    pass


class DockerTimeoutError(DockerError):
    """Docker operation timed out."""

    pass


class DockerHealthCheckError(DockerError):
    """Container health check failed."""

    pass


class DockerImageNotFoundError(DockerError):
    """Docker image not found."""

    pass


class DockerOperations:
    """
    Docker container operations with robust error handling.

    Follows SOLID principles:
    - Single Responsibility: Docker interactions only
    - Dependency Inversion: Depends on pathlib.Path, not concrete paths

    Best practices from 2024 research:
    - Use subprocess.run() with timeout
    - Handle 3 critical exceptions: TimeoutExpired, CalledProcessError, FileNotFoundError
    - Enable BuildKit for 30-50% faster builds
    """

    def __init__(self, image_name: str = "vntyper:latest", default_timeout: int = 1800):
        """
        Initialize Docker operations.

        Args:
            image_name: Docker image name
            default_timeout: Default timeout for operations (seconds)
                           Default is 1800s (30 minutes) to accommodate long-running
                           tests like adVNTR which can take ~25 minutes

        Raises:
            DockerNotFoundError: If Docker not installed
        """
        self.image_name = image_name
        self.default_timeout = default_timeout
        self._verify_docker_available()

    def _verify_docker_available(self) -> None:
        """
        Verify docker command is available.

        Raises:
            DockerNotFoundError: If docker not installed
        """
        try:
            subprocess.run(
                ["docker", "--version"], capture_output=True, text=True, check=True, timeout=5
            )
            log.info(f"Docker available for image: {self.image_name}")
        except FileNotFoundError as e:
            raise DockerNotFoundError(
                "Docker not found. Install from: https://docs.docker.com/get-docker/"
            ) from e
        except subprocess.TimeoutExpired:
            log.warning("Docker version check slow (may indicate Docker daemon issues)")
        except subprocess.CalledProcessError as e:
            log.warning(f"Docker version check failed: {e.stderr}")

    def build_image(
        self, dockerfile: Path, context: Path = Path("."), timeout: Optional[int] = None
    ) -> bool:
        """
        Build Docker image with BuildKit.

        Args:
            dockerfile: Path to Dockerfile
            context: Build context directory
            timeout: Build timeout (default: 600s for builds)

        Returns:
            True if build succeeds

        Raises:
            DockerBuildError: If build fails
            DockerTimeoutError: If build exceeds timeout
        """
        timeout = timeout or 600  # Builds can take longer

        if not dockerfile.exists():
            raise DockerBuildError(f"Dockerfile not found: {dockerfile}")

        log.info(f"Building Docker image: {self.image_name}")
        log.info(f"  Dockerfile: {dockerfile}")
        log.info(f"  Context: {context}")
        log.info(f"  Timeout: {timeout}s")

        cmd = ["docker", "build", "-f", str(dockerfile), "-t", self.image_name, str(context)]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=True,
                env={**os.environ, "DOCKER_BUILDKIT": "1"},  # Enable BuildKit
            )

            # Log build output (last 10 lines for debugging)
            output_lines = result.stdout.strip().split("\n")
            log.info("Build output (last 10 lines):")
            for line in output_lines[-10:]:
                if line.strip():
                    log.info(f"  {line}")

            log.success(f"Docker image built: {self.image_name}")
            return True

        except subprocess.TimeoutExpired as e:
            raise DockerTimeoutError(
                f"Docker build timed out after {timeout}s. "
                "Try increasing timeout or check network/disk speed."
            ) from e
        except subprocess.CalledProcessError as e:
            # Extract useful error info from stderr
            error_lines = e.stderr.strip().split("\n")[-10:]
            raise DockerBuildError("Docker build failed:\n" + "\n".join(error_lines)) from e

    def run_command(
        self,
        command: list[str],
        volumes: Optional[dict[Path, Path]] = None,
        timeout: Optional[int] = None,
        capture_output: bool = True,
        entrypoint: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        """
        Run command in Docker container.

        Args:
            command: Command to run in container (e.g., ["vntyper", "--version"])
            volumes: Volume mounts {host_path: container_path}
            timeout: Execution timeout (default: self.default_timeout)
            capture_output: Capture stdout/stderr (False for interactive)
            entrypoint: Override container entrypoint (e.g., "bash")

        Returns:
            subprocess.CompletedProcess with stdout/stderr

        Raises:
            DockerExecutionError: If command fails
            DockerTimeoutError: If execution exceeds timeout
        """
        timeout = timeout or self.default_timeout

        # Build docker run command
        docker_cmd = ["docker", "run", "--rm"]

        # Add volume mounts (paths should already be resolved by caller)
        if volumes:
            for host_path, container_path in volumes.items():
                docker_cmd.extend(["-v", f"{host_path}:{container_path}"])

        # Add custom entrypoint if specified
        if entrypoint:
            docker_cmd.extend(["--entrypoint", entrypoint])

        docker_cmd.append(self.image_name)
        docker_cmd.extend(command)

        # Log command (truncated for readability)
        cmd_str = " ".join(str(x) for x in docker_cmd[:7])
        log.info(f"Running: {cmd_str}... (truncated)")

        try:
            result = subprocess.run(
                docker_cmd, capture_output=capture_output, text=True, timeout=timeout, check=True
            )

            return result

        except subprocess.TimeoutExpired as e:
            raise DockerTimeoutError(
                f"Container execution timed out after {timeout}s. "
                "Pipeline may be hanging or data too large."
            ) from e
        except subprocess.CalledProcessError as e:
            # Provide helpful error message
            stderr_preview = e.stderr[:500] if e.stderr else "No error output"
            raise DockerExecutionError(
                f"Container command failed (exit code {e.returncode}):\n"
                f"Command: {' '.join(command)}\n"
                f"Error output:\n{stderr_preview}"
            ) from e

    def run_pipeline(
        self,
        bam_file: Path,
        output_dir: Path,
        reference_assembly: str = "hg19",
        threads: int = 4,
        extra_modules: Optional[list[str]] = None,
        timeout: Optional[int] = None,
    ) -> subprocess.CompletedProcess:
        """
        Run VNtyper pipeline in container.

        Args:
            bam_file: Input BAM file path
            output_dir: Output directory (host path)
            reference_assembly: Reference assembly (hg19, hg38, GRCh37, GRCh38)
            threads: Number of threads
            extra_modules: Extra modules (e.g., ["advntr"])
            timeout: Pipeline timeout (default: 1800s / 30 minutes)
                    Note: adVNTR tests can take ~25 minutes

        Returns:
            subprocess.CompletedProcess

        Raises:
            DockerExecutionError: If pipeline fails
            DockerTimeoutError: If pipeline times out
        """
        # Setup volumes (paths should already be resolved by caller)
        input_dir = bam_file.parent
        output_dir.mkdir(parents=True, exist_ok=True)

        volumes = {
            input_dir: Path("/opt/vntyper/input"),
            output_dir: Path("/opt/vntyper/output"),
        }

        # Build vntyper command
        vntyper_cmd = [
            "vntyper",
            "pipeline",
            "--bam",
            f"/opt/vntyper/input/{bam_file.name}",
            "-o",
            "/opt/vntyper/output/result",
            "--reference-assembly",
            reference_assembly,
            "--threads",
            str(threads),
            "--fast-mode",
            "--keep-intermediates",
        ]

        if extra_modules:
            vntyper_cmd.extend(["--extra-modules", ",".join(extra_modules)])
            if "advntr" in extra_modules:
                vntyper_cmd.extend(["--advntr-max-coverage", "300"])

        log.info(f"Pipeline: BAM={bam_file.name}, ref={reference_assembly}, threads={threads}")
        if extra_modules:
            log.info(f"  Extra modules: {', '.join(extra_modules)}")

        return self.run_command(
            command=vntyper_cmd, volumes=volumes, timeout=timeout or self.default_timeout
        )

    def check_health(self) -> bool:
        """
        Run container health checks.

        Checks:
        1. Container starts and vntyper accessible
        2. Java runtime available (for Kestrel)
        3. Bioinformatics tools available

        Returns:
            True if all health checks pass

        Raises:
            DockerHealthCheckError: If any health check fails
        """
        log.info("Running container health checks...")

        # Check 1: VNtyper version
        try:
            result = self.run_command(["vntyper", "--version"], timeout=10)
            log.success("✓ VNtyper CLI accessible")
        except Exception as e:
            raise DockerHealthCheckError(f"VNtyper not accessible: {e}") from e

        # Check 2: Java runtime
        try:
            result = self.run_command(
                ["-c", "conda run -n vntyper java -version"],
                entrypoint="bash",
                timeout=10
            )
            # Java version goes to stderr
            if "openjdk" in result.stderr.lower() or "openjdk" in result.stdout.lower():
                log.success("✓ Java runtime available")
            else:
                raise DockerHealthCheckError("Java runtime not openjdk")
        except DockerExecutionError as e:
            # Java might still work, check stderr
            if "openjdk" in str(e).lower():
                log.success("✓ Java runtime available")
            else:
                raise DockerHealthCheckError(f"Java runtime check failed: {e}") from e

        # Check 3: Bioinformatics tools
        tools = ["samtools", "bwa", "fastp", "bcftools"]
        for tool in tools:
            try:
                self.run_command(
                    ["-c", f"conda run -n vntyper which {tool}"],
                    entrypoint="bash",
                    timeout=10
                )
                log.success(f"✓ {tool} available")
            except Exception as e:
                raise DockerHealthCheckError(f"{tool} check failed: {e}") from e

        log.success("All health checks passed")
        return True

    def get_image_info(self) -> dict[str, str]:
        """
        Get Docker image information.

        Returns:
            Dictionary with image metadata (size, created)

        Raises:
            DockerImageNotFoundError: If image not found
        """
        try:
            result = subprocess.run(
                [
                    "docker",
                    "images",
                    "--format",
                    "{{.Size}}\t{{.CreatedAt}}",
                    self.image_name,
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=5,
            )

            if result.stdout.strip():
                parts = result.stdout.strip().split("\t", 1)
                size = parts[0] if len(parts) > 0 else "unknown"
                created = parts[1] if len(parts) > 1 else "unknown"
                return {"size": size, "created": created}
            else:
                raise DockerImageNotFoundError(f"Image not found: {self.image_name}")

        except subprocess.CalledProcessError as e:
            raise DockerImageNotFoundError(f"Failed to get image info: {self.image_name}") from e

    def image_exists(self) -> bool:
        """
        Check if Docker image exists locally.

        Returns:
            True if image exists
        """
        try:
            subprocess.run(
                ["docker", "image", "inspect", self.image_name],
                capture_output=True,
                check=True,
                timeout=5,
            )
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return False
