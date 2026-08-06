"""The upload-filename boundary for the web service.

`UploadFile.filename` is supplied by the client and is never sanitised by
Starlette. It is used to build the destination path of the uploaded file and is
carried onwards into the analysis pipeline, so it has to be constrained here, at
the boundary, before it reaches any filesystem or tooling call. Downstream
layers assume a plain, well-formed name and must not be relied upon to repair a
hostile one.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.uploads` is importable here.
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from app.uploads import ALIGNMENT_EXTENSIONS, INDEX_EXTENSIONS, safe_upload_path  # noqa: E402

# Rejected names, each paired with the reason it must not be accepted.
UNSAFE_NAMES = [
    ("../../../../etc/cron.d/x.bam", "relative traversal out of the job directory"),
    ("../x.bam", "a single traversal step is no more acceptable than four"),
    ("/etc/cron.d/x.bam", "absolute path: os.path.join discards the base entirely"),
    ("a.bam; touch marker", "command separator"),
    ("a.bam && touch marker", "command chaining"),
    ("a.bam | tee marker", "pipe"),
    ("$(touch marker).bam", "command substitution"),
    ("`touch marker`.bam", "backtick command substitution"),
    ("a b.bam", "a space breaks unquoted interpolation even without malice"),
    ("a'b.bam", "single quote"),
    ('a"b.bam', "double quote"),
    ("a*.bam", "glob wildcard"),
    ("~/a.bam", "tilde is expanded by the shell"),
    ("-a.bam", "a leading dash can be read as an option by a downstream tool"),
    ("+a.bam", "a plus is permitted inside a name but not as its first character"),
    ("", "an empty name yields the directory itself, not a file"),
    (None, "no filename part at all"),
    ("x.txt", "not an alignment file"),
    ("x.bam.exe", "the allowed extension must be the final one, not merely present"),
    (".bam", "extension with no stem"),
    (".", "the job directory itself"),
    ("..", "the parent of the job directory"),
    ("a\x00.bam", "NUL byte: truncated at the syscall boundary, so the tail is invisible"),
    ("a.bam\n", "trailing newline: '$' would accept this, the anchor must not"),
    ("a\nb.bam", "embedded newline"),
    ("a.bam\r", "trailing carriage return"),
    ("a.bam ", "trailing space"),
    ("..\\..\\x.bam", "backslash separators: os.path.basename is a no-op under POSIX"),
    ("‮emag.bam", "right-to-left override disguises how the name reads"),
    ("a" * 252 + ".bam", "256 characters, past the 255-byte filesystem limit"),
    ("x.bai", "an index extension is not an alignment file"),
    # The allowlist is ASCII-only. Matching an extension case-insensitively is
    # exactly where that can be lost: Python's IGNORECASE folds across Unicode
    # unless it is told not to, so `[A-Za-z]` starts matching these four.
    ("Müller.bam", "non-ASCII: normalisation and encoding questions with no clinical benefit"),
    ("KKlein.bam", "U+212A KELVIN SIGN case-folds to 'k' under a Unicode-aware match"),
    ("ſample.bam", "U+017F LATIN SMALL LETTER LONG S case-folds to 's'"),
    ("İ.bam", "U+0130 LATIN CAPITAL LETTER I WITH DOT ABOVE case-folds to 'i'"),
    ("ı.bam", "U+0131 LATIN SMALL LETTER DOTLESS I case-folds to 'i'"),
]


@pytest.mark.parametrize("filename,reason", UNSAFE_NAMES, ids=[repr(name) for name, _ in UNSAFE_NAMES])
def test_rejects_unsafe_filenames(tmp_path: Path, filename, reason: str) -> None:
    """Every unacceptable name raises instead of producing a path.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: The client-supplied name under test.
        reason: Why this name must be rejected; reported on failure.
    """
    with pytest.raises(ValueError):
        safe_upload_path(str(tmp_path), filename), reason


@pytest.mark.parametrize(
    "filename",
    [
        "sample.bam",
        "S1_2.cram",
        "a-b.c_1.bam",
        "a..bam",  # consecutive dots inside the stem are fine; they cannot traverse
        "0.bam",  # a digit is a valid first character
        "a" * 251 + ".bam",  # exactly 255 characters: the limit itself is allowed
        # A plus is neither a path separator nor a shell metacharacter, and it is
        # how tumour/normal pairs are conventionally named.
        "sample+tumor.bam",
        # Sequencers and LIMS exports routinely upper-case the extension. The
        # name is stored exactly as sent; only the *match* ignores case.
        "SAMPLE.BAM",
        "X.BAM",
        "x.Bam",
        "S1.CRAM",
        "S1.Cram",
    ],
)
def test_accepts_and_contains_safe_filenames(tmp_path: Path, filename: str) -> None:
    """Ordinary names are accepted unchanged and land directly in the job directory.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: A name that must survive the boundary untouched.
    """
    result = Path(safe_upload_path(str(tmp_path), filename))

    assert result.parent == tmp_path
    assert result.name == filename


@pytest.mark.parametrize("filename", ["SAMPLE.BAM", "Sample.Bam", "sample+tumor.bam", "S1.CRAM.CRAI"])
def test_an_accepted_name_is_never_case_folded_or_rewritten(tmp_path: Path, filename: str) -> None:
    """Matching an extension without regard to case does not change the name.

    Storing `SAMPLE.BAM` as `sample.bam` would hand the worker a path the client
    never named, and two submissions differing only in case would collide.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: A name that must be stored exactly as sent.
    """
    extensions = INDEX_EXTENSIONS if filename.lower().endswith(("bai", "crai")) else ALIGNMENT_EXTENSIONS

    assert Path(safe_upload_path(str(tmp_path), filename, extensions)).name == filename


@pytest.mark.parametrize("filename", ["sample.BAI", "S1.cram.CRAI", "sample+tumor.bam.bai"])
def test_the_index_slot_matches_its_extensions_the_same_way(tmp_path: Path, filename: str) -> None:
    """Both upload slots apply one rule, so neither can drift from the other.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: An index name that must be accepted.
    """
    assert Path(safe_upload_path(str(tmp_path), filename, INDEX_EXTENSIONS)).name == filename


def test_result_is_always_inside_the_job_directory(tmp_path: Path) -> None:
    """The containment property, stated independently of the allowlist.

    Args:
        tmp_path: Stand-in for the per-job input directory.
    """
    for name in ("ok.bam", "ok.cram"):
        resolved = Path(safe_upload_path(str(tmp_path), name)).resolve()
        assert tmp_path.resolve() in resolved.parents


def test_accepted_name_is_returned_verbatim_never_repaired(tmp_path: Path) -> None:
    """A hostile name is refused outright rather than quietly rewritten.

    Stripping the directory part of `../../evil.bam` and storing `evil.bam`
    would turn a rejected request into a silently accepted one, so the caller
    would never learn the client sent something it should not have.

    Args:
        tmp_path: Stand-in for the per-job input directory.
    """
    with pytest.raises(ValueError):
        safe_upload_path(str(tmp_path), "../../evil.bam")

    assert not (tmp_path / "evil.bam").exists()
    assert list(tmp_path.iterdir()) == []


def test_error_message_names_the_accepted_extensions(tmp_path: Path) -> None:
    """The rejection explains what is acceptable, without echoing the input back.

    Args:
        tmp_path: Stand-in for the per-job input directory.
    """
    with pytest.raises(ValueError) as excinfo:
        safe_upload_path(str(tmp_path), "totally-not-a-bam.txt")

    message = str(excinfo.value)
    assert ".bam" in message
    assert ".cram" in message
    assert "totally-not-a-bam" not in message


def test_error_message_says_which_characters_are_acceptable(tmp_path: Path) -> None:
    """A caller whose name is refused is told what to send instead.

    The ASCII-only rule is the one a legitimate name is most likely to fall foul
    of, and it is not guessable from the extension list alone, so the message
    states it rather than leaving the caller to experiment.

    Args:
        tmp_path: Stand-in for the per-job input directory.
    """
    with pytest.raises(ValueError) as excinfo:
        safe_upload_path(str(tmp_path), "Müller.bam")

    message = str(excinfo.value)
    assert "ASCII" in message
    assert "Müller" not in message


@pytest.mark.parametrize("filename", ["sample.bam.bai", "sample.bai", "S1.cram.crai", "S1.crai"])
def test_index_extensions_are_accepted_for_the_index_slot(tmp_path: Path, filename: str) -> None:
    """The index upload accepts the index extensions the pipeline looks for.

    `tasks.run_vntyper_job` looks for `<bam>.bai` beside the alignment, and the
    input-format docs tell users to bring `.bam.bai` or `.bai`, so refusing
    those would break a working feature.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: An index name that must be accepted.
    """
    result = Path(safe_upload_path(str(tmp_path), filename, INDEX_EXTENSIONS))

    assert result.parent == tmp_path
    assert result.name == filename


@pytest.mark.parametrize("filename", ["sample.bam", "sample.cram", "../x.bai", "x.bai.txt"])
def test_index_slot_rejects_everything_else(tmp_path: Path, filename: str) -> None:
    """The two slots do not accept each other's extensions, and the rest still fails.

    Args:
        tmp_path: Stand-in for the per-job input directory.
        filename: A name the index slot must reject.
    """
    with pytest.raises(ValueError):
        safe_upload_path(str(tmp_path), filename, INDEX_EXTENSIONS)


def test_extension_sets_are_the_ones_the_endpoint_relies_on() -> None:
    """Pin the two allowlists so a future edit to either is a deliberate act."""
    assert ALIGNMENT_EXTENSIONS == ("bam", "cram")
    assert INDEX_EXTENSIONS == ("bai", "crai")


# ---------------------------------------------------------------------------
# Endpoint level: the constraint has to be wired in, not merely available.
# ---------------------------------------------------------------------------


def _job_files(input_dir: Path):
    """List every file written under the job input tree.

    Args:
        input_dir: The service's configured input directory.

    Returns:
        list[Path]: Files (not directories) found anywhere beneath it.
    """
    return [path for path in input_dir.rglob("*") if path.is_file()]


def test_traversal_filename_is_rejected_at_the_endpoint(client, web_app, tmp_path: Path) -> None:
    """A traversing upload name gets a 400 and writes nothing.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, sharing this test's `tmp_path`.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("../../evil.bam", b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 400
    assert ".bam" in response.json()["detail"]
    assert _job_files(tmp_path / "input") == []
    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_not_called()


@pytest.mark.parametrize("filename", ["a.bam; touch marker", "/etc/cron.d/x.bam", "x.txt"])
def test_other_unsafe_upload_names_are_rejected_at_the_endpoint(client, tmp_path: Path, filename: str) -> None:
    """The endpoint refuses the whole rejected class, not just traversal.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The same directory the fixture configured as the input tree.
        filename: The client-supplied name under test.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": (filename, b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 400
    assert _job_files(tmp_path / "input") == []


def test_unsafe_index_filename_is_rejected_at_the_endpoint(client, tmp_path: Path) -> None:
    """A good BAM with an unacceptable index name is refused before anything is written.

    Both names are checked before the first byte is stored, so a rejected
    submission leaves no half-written job behind.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        files={
            "bam_file": ("ok.bam", b"x", "application/octet-stream"),
            "bai_file": ("../../evil.bai", b"x", "application/octet-stream"),
        },
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 400
    assert not (tmp_path / "input" / "evil.bai").exists()
    assert _job_files(tmp_path / "input") == []


def test_valid_upload_still_succeeds_and_lands_in_the_job_directory(client, web_app, tmp_path: Path) -> None:
    """The happy path is unchanged: the file is stored under its own name.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, sharing this test's `tmp_path`.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        files={
            "bam_file": ("sample.bam", b"bamdata", "application/octet-stream"),
            "bai_file": ("sample.bam.bai", b"baidata", "application/octet-stream"),
        },
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    job_id = response.json()["job_id"]
    job_input_dir = tmp_path / "input" / job_id
    assert (job_input_dir / "sample.bam").read_bytes() == b"bamdata"
    assert (job_input_dir / "sample.bam.bai").read_bytes() == b"baidata"
    web_app.run_vntyper_job.delay.assert_called_once()
    assert web_app.run_vntyper_job.delay.call_args.kwargs["bam_path"] == str(job_input_dir / "sample.bam")


_BOUNDARY = "vntyperuploadboundary"


def _multipart_with_empty_index_filename() -> bytes:
    """Build a request body whose index part carries an empty filename.

    `httpx` omits the `filename` attribute entirely when it is empty, which
    Starlette then parses as an ordinary form field rather than a file, so the
    high-level `files=` helper cannot reach this branch at all. The body is
    written out by hand so the test exercises what a real client can send.

    Returns:
        bytes: A complete `multipart/form-data` body.
    """
    lines = [
        f"--{_BOUNDARY}",
        'Content-Disposition: form-data; name="thread"',
        "",
        "1",
        f"--{_BOUNDARY}",
        'Content-Disposition: form-data; name="reference_assembly"',
        "",
        "hg19",
        f"--{_BOUNDARY}",
        'Content-Disposition: form-data; name="bam_file"; filename="ok.bam"',
        "Content-Type: application/octet-stream",
        "",
        "x",
        f"--{_BOUNDARY}",
        'Content-Disposition: form-data; name="bai_file"; filename=""',
        "Content-Type: application/octet-stream",
        "",
        "",
        f"--{_BOUNDARY}--",
        "",
    ]
    return "\r\n".join(lines).encode()


def test_empty_bai_part_does_not_500(client, tmp_path: Path) -> None:
    """An index part with an empty filename is ignored rather than crashing.

    Without the guard this reaches `open()` on the job directory itself and the
    request fails with a 500.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart_with_empty_index_filename(),
        headers={"Content-Type": f"multipart/form-data; boundary={_BOUNDARY}"},
    )

    assert response.status_code == 200
    assert [path.name for path in _job_files(tmp_path / "input")] == ["ok.bam"]


def test_absent_bai_part_still_submits_the_job(client, tmp_path: Path) -> None:
    """Omitting the index part entirely remains valid; the pipeline builds one.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("ok.bam", b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    assert [path.name for path in _job_files(tmp_path / "input")] == ["ok.bam"]
