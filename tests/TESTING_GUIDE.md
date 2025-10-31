# VNtyper Testing Guide

This guide explains how to run VNtyper tests with enhanced logging and progress visibility.

## Quick Start

```bash
# Install dev dependencies (includes pytest-timeout)
pip install -e .[dev]

# Run all tests with live logging (default now)
pytest

# Run integration tests with live logging
pytest -m integration

# Run unit tests (fast)
pytest -m unit
```

## Enhanced Logging Features

The test suite now includes several enhancements to improve visibility during slow integration tests:

### 1. **Live Logging** (Enabled by Default)

All test logs are now displayed in real-time as tests execute. You'll see:
- Test start/stop markers
- Pipeline execution progress
- Individual test durations
- Clear visual separators

**What you'll see:**
```
10:25:30 [    INFO] ================================================================================
10:25:30 [    INFO] STARTING: tests/integration/test_pipeline_integration.py::test_bam_input...
10:25:30 [    INFO] ================================================================================
10:25:31 [    INFO] EXECUTING: tests/integration/test_pipeline_integration.py::test_bam_input...
10:28:45 [    INFO] DURATION: 194.32s for tests/integration/test_pipeline_integration.py::test_bam_input...
10:28:45 [    INFO] ✓ PASSED (194.32s): tests/integration/test_pipeline_integration.py::test_bam_input...
```

### 2. **Test Durations Report**

The 10 slowest tests are automatically reported at the end:

```
============================== slowest 10 durations ================================
194.32s call     tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf]
142.67s call     tests/integration/test_pipeline_integration.py::test_advntr_input[example_a5c1]
...
```

### 3. **Persistent Log File**

All test output is saved to `tests/pytest.log` for later analysis:

```bash
# View the log file
cat tests/pytest.log

# Follow the log in real-time (in another terminal)
tail -f tests/pytest.log
```

### 4. **Progress Indicators**

At the start of each test session, you'll see:
```
10:20:15 [    INFO] ================================================================================
10:20:15 [    INFO] PYTEST SESSION STARTED
10:20:15 [    INFO] ================================================================================
10:20:16 [    INFO] Collected 15 test(s)
10:20:16 [    INFO]   - Integration tests: 10
10:20:16 [    INFO]   - Unit tests: 5
```

## Test Commands

### Basic Commands

```bash
# Run all tests with default configuration
pytest

# Run tests verbosely (shows test names as they run)
pytest -v

# Run tests with even more detail (captures disabled)
pytest -v -s

# Run tests quietly (minimal output)
pytest -q
```

### Filtered Test Execution

```bash
# Run only integration tests
pytest -m integration

# Run only unit tests
pytest -m unit

# Run specific test file
pytest tests/integration/test_pipeline_integration.py

# Run specific test
pytest tests/integration/test_pipeline_integration.py::test_bam_input_with_kestrel_checks
```

### Coverage Reports

```bash
# Run tests with coverage
pytest --cov=vntyper

# Run tests with coverage HTML report
pytest --cov=vntyper --cov-report=html

# View coverage report
open htmlcov/index.html
```

### Advanced Options

```bash
# Run tests with timeout protection (requires pytest-timeout)
pytest --timeout=600  # 10 minute timeout per test

# Run tests in parallel (requires pytest-xdist)
pytest -n auto

# Run only tests that failed last time
pytest --lf

# Run failed tests first, then others
pytest --ff

# Stop on first failure
pytest -x

# Stop after N failures
pytest --maxfail=3
```

## Logging Levels

You can control the logging level for different scenarios:

### During Development (More Detail)

```bash
# Show DEBUG level logs in console
pytest --log-cli-level=DEBUG

# Or edit pytest.ini temporarily:
# log_cli_level = DEBUG
```

### During CI/CD (Less Noise)

```bash
# Show only WARNING and above
pytest --log-cli-level=WARNING

# Disable live logging entirely
pytest --log-cli=false
```

## Troubleshooting Slow/Stuck Tests

### Real-time Monitoring

If a test seems stuck, in another terminal:

```bash
# Monitor the log file in real-time
tail -f tests/pytest.log

# Monitor just the key events
tail -f tests/pytest.log | grep -E "(STARTING|EXECUTING|PASSED|FAILED|DURATION)"
```

### Enable Test Timeout

Uncomment these lines in `pytest.ini` to enable automatic timeout:

```ini
# Timeout for integration tests (10 minutes)
timeout = 600
timeout_method = thread
```

This will automatically fail tests that run longer than 10 minutes.

### Check Pipeline Logs

Integration tests run the full VNtyper pipeline, which has its own logging:

```bash
# The pipeline logs to stderr, which pytest captures
# Use -s to see pipeline logs in real-time:
pytest -v -s -m integration
```

## Best Practices

### 1. **Use Live Logging for Integration Tests**

Integration tests are slow, so live logging helps you track progress:

```bash
pytest -m integration  # Live logging enabled by default
```

### 2. **Use Quiet Mode for Unit Tests**

Unit tests are fast, so you can use quieter output:

```bash
pytest -m unit -q  # Quiet mode
```

### 3. **Check Log File After Failures**

The log file contains full details:

```bash
# Run tests
pytest -m integration

# If tests fail, check the log
less tests/pytest.log

# Search for errors
grep -A 10 "FAILED" tests/pytest.log
```

### 4. **Use Markers for Selective Testing**

During development, run only the tests you need:

```bash
# Run specific subset
pytest -k "example_66bf"

# Run tests matching pattern
pytest -k "kestrel"

# Exclude tests matching pattern
pytest -k "not advntr"
```

## Performance Tips

### Speed Up Test Runs

1. **Run unit tests first** (they're fast):
   ```bash
   pytest -m unit
   ```

2. **Use parallel execution** (requires pytest-xdist):
   ```bash
   pip install pytest-xdist
   pytest -n auto -m unit
   ```

3. **Skip slow tests during development**:
   ```bash
   pytest -k "not slow"
   ```

### Reduce Log Verbosity for Speed

Excessive logging can slow down tests. For CI/CD:

```bash
# Minimal logging
pytest --log-cli-level=WARNING -q
```

## Configuration Files

### pytest.ini

Main pytest configuration with:
- Live logging settings
- Log file settings
- Default options
- Timeout configuration

### tests/conftest.py

Pytest hooks that provide:
- Test progress tracking
- Timing information
- Status indicators (✓/✗/⊘)
- Session summaries

## Example Test Run Output

Here's what you'll see when running integration tests:

```
================================ test session starts =================================
platform linux -- Python 3.10.14, pytest-8.4.2, pluggy-1.6.0
rootdir: /mnt/c/development/hassansaei/VNtyper
configfile: pytest.ini
plugins: cov-6.2.1, timeout-2.3.1
collected 10 items

10:20:16 [    INFO] ================================================================================
10:20:16 [    INFO] PYTEST SESSION STARTED
10:20:16 [    INFO] ================================================================================
10:20:16 [    INFO] Collected 10 test(s)
10:20:16 [    INFO]   - Integration tests: 10

10:20:17 [    INFO] ================================================================================
10:20:17 [    INFO] STARTING: tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf] [integration]
10:20:17 [    INFO] ================================================================================
10:20:18 [    INFO] EXECUTING: tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf]
10:23:32 [    INFO] DURATION: 194.32s for tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf]
10:23:32 [    INFO] ✓ PASSED (194.32s): tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf]

[... more tests ...]

10:35:45 [    INFO] ================================================================================
10:35:45 [    INFO] PYTEST SESSION FINISHED
10:35:45 [    INFO] Exit status: 0
10:35:45 [    INFO] ================================================================================

============================== slowest 10 durations =================================
194.32s call     tests/integration/test_pipeline_integration.py::test_bam_input[example_66bf]
142.67s call     tests/integration/test_pipeline_integration.py::test_advntr_input[example_a5c1]
...

============================= 10 passed in 925.47s (15:25) ==========================
```

## CI/CD Integration

For continuous integration, use these options:

```bash
# Recommended CI command
pytest -v -m integration --log-cli-level=INFO --maxfail=3 --timeout=600

# Minimal output for logs
pytest -q -m integration --log-cli-level=WARNING

# With coverage
pytest --cov=vntyper --cov-report=xml --cov-report=term -m integration
```

## Getting Help

```bash
# Show all pytest options
pytest --help

# Show available fixtures
pytest --fixtures

# Show available markers
pytest --markers

# Collect tests without running
pytest --collect-only
```

## Further Reading

- [Pytest Documentation](https://docs.pytest.org/)
- [Pytest Logging](https://docs.pytest.org/en/stable/how-to/logging.html)
- [Pytest Fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html)
- [Pytest Markers](https://docs.pytest.org/en/stable/example/markers.html)
