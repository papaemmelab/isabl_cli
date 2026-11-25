# Changelog

All notable changes to isabl_cli are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.3] - Unreleased

### Added

- Every API request now sends an `X-Isabl-CLI-Version: <version>` header and a
  `User-Agent` of `isabl_cli/<version> python-requests/<x.y>`, so the API can tell
  which CLI version made a request. This includes the login request and
  `api_request(authenticate=False)`. Until now every client looked like plain
  `python-requests`. This lets the API refuse writes from installs older than its
  `MIN_CLI_VERSION`.
- While pytest is loaded, isabl_cli refuses `POST`, `PUT`, `PATCH` and `DELETE`
  requests to any API that isn't on localhost, and raises the new
  `TestSessionWriteError` before anything is sent. This applies to every test suite
  that uses isabl_cli. Reads are unaffected. Set `ISABL_ALLOW_TEST_WRITES_TO` to the
  URL of a disposable remote test server to allow writes to it. (#89)

### Changed

- Analyses marked `FINISHED` now get their results and storage usage patched, just as
  `SUCCEEDED` ones do. Before, analyses run by non-admin users had no results in the
  web UI until an admin ran `process-finished`. (#88)

### Testing

- `pytest.ini` sets `ISABL_API_URL` to `http://localhost:8000/api/v1/`, overriding
  whatever the shell exports. (#90)
- `tests/conftest.py` refuses to run unless `ISABL_API_URL` is local. It checks on
  import, at session start and before each test. (#90)
- The README explains how to run the tests against a local API, and against a remote
  test server. (#90)

## [2.0.2] - 2026-04-24

### Added

- `delete_instance` can be imported from `isabl_cli`, like the other API helpers. (#75)
- A result specification can set `store_as_bam` to save that result as the target
  experiment's BAM file. (#78)
- A dependency in `dependencies_results` can set `app_assembly` (for example `GRCh37`)
  to take results only from that assembly's version of the app. (#79)
- A pattern result specification can set `optional`. When no file matches, the result
  is empty instead of raising `FileNotFoundError`.
- Merged analyses record how many analyses were merged in `data.merged_analyses`.

### Changed

- User settings, including the login token, are stored separately for each
  `ISABL_API_URL`, so one machine can stay logged in to several APIs. Tokens saved by
  earlier versions aren't read, so users have to `isabl login` again after upgrading.
  (#82)
- When several files match a result pattern, the newest is used.
- "Checking for existing analyses..." and section titles go to stdout again instead of
  stderr.

### Fixed

- The latest dependency analysis is picked by numeric primary key. Before, keys were
  compared as strings, so analysis 9 was treated as newer than analysis 10.
- Project-level merges are skipped while any of the project's analyses are `FINISHED`
  but not yet processed, as they already were for `STARTED` and `SUBMITTED` ones.
  (fixes #46)

### Development

- CI moved from Travis to GitHub Actions. (#77, #80)
- The workflow that publishes tags to PyPI was fixed. (#73)
- The README covers the development and test setup (#83) and links to the Isabl MCP
  server.

## [2.0.1] - 2025-04-18

### Changed

- `setup.py` declares the package metadata itself instead of reading `setup.json`, and
  uses the README as the PyPI description.
- The README was updated.

## [2.0.0] - 2025-04-18

This is the first release published to PyPI from GitHub Actions. It covers five years
of changes since 1.0.0.

### Added

- Data import:
  - The data importer is explicit. (#30 in the old `isabl-io/cli` repository)
  - It can copy local files instead of symlinking them, including from YAML imports.
  - It follows symlinked directories.
  - It checks that files are readable before importing them.
- `EXTRA_RAW_DATA_FORMATS` setting to register more raw data formats. (#16)
- Applications can depend on other applications' results:
  - Declare these in `dependencies_results`, matching the app by class or by
    `app_name`, with `app_version` either a version or `latest`. (#22, #51)
  - A dependency can set `status`, with several statuses joined by `/`.
  - Set `linked` to `False` to use a dependency's results without linking the
    analysis as a dependency.
  - `get_dependency_analyses_option` accepts extra filters.
  - `filter_analyses` also accepts extra filters.
- Result specifications can set a `pattern`, and an `exclude`, to take a result from
  the first matching file in the analysis directory. (#71)
- `REJECTED` analysis status and a hidden `reject-analysis` command that rejects an
  analysis and stores the reason in its notes. (#58)
- Validators:
  - `validate_same_platform` and `validate_individuals`, for unmatched pipelines that
    set `IS_UNMATCHED`;
  - `validate_source`, to check the samples' source material.
- `get-metadata --all` prints every field, and `--pretty` lines up the columns.
  (experimental, #43)
- `process-finished --force` patches analyses again even if they were already
  processed, and the command reports how many analyses it skipped. (#36, #37)
- `--commit` shows a summary of what will be run. (#53)
- Applications email the analysts of a project when it has errors, through the new
  `send_error_email` API helper.
- Batch systems:
  - Slurm and SGE are supported alongside LSF.
  - Slurm accepts `extra_args` for every job, including its exit and clean jobs. (#61)
  - Before submitting, isabl_cli checks that LSF and Slurm are installed.
- Usage events are sent to analytics for each command and app run.

### Changed

- Error messages go to stderr. (#44)
- Analyses that already exist are validated again in case their data changed. (#14)
- If one app fails to import, the others still work. (#28)
- Species is only validated when the app has an assembly. (#54)
- Batch submission:
  - Local and Slurm jobs run with `/bin/bash`. (#60, #62, #63)
  - The Slurm `requirements` setting takes precedence over `extra_args`.
- Execution logs can be shown as ANSI in the web UI.
- The import warning threshold went up to 800 files.
- Dependencies were relaxed, and Python 3.10 is supported.

### Fixed

- Restarting an analysis started by a different user is refused unless `--force` is
  passed. (#50)
- `process-finished`:
  - no longer races with itself;
  - restores the analysis's tags when patching it fails, so the analysis isn't left
    marked as in progress.
- `submit_local` no longer sets `SUCCEEDED` regardless of `ADMIN_USER`.
- LSF job arrays with more than 10,000 jobs are split.
- Fixed:
  - the rsync command (#49);
  - querying dependency analyses;
  - file permissions when patching results (#65);
  - choosing the latest dependency analysis (#64).
- App notifications are skipped when a project has no analysts.

## [1.0.0] - 2020-01-27

This is the first tagged release, covering the work since May 2018. It provides:

- The `isabl` command line tool, with these commands:
  - `login`;
  - `get-*` commands for metadata, data, BAMs, BED files, paths, outdirs, results,
    references and counts;
  - `patch-status`, `patch-results` and `process-finished`;
  - project and individual merges;
  - signals.
- A Python client for the Isabl API, with helpers such as `get_instances`,
  `create_instance` and `patch_instance`.
- `AbstractApplication`, the framework for writing Isabl applications, along with
  validators.
- Running analyses locally or on LSF.
- Data import into Isabl.
- Test factories for Isabl models.

[2.0.3]: https://github.com/papaemmelab/isabl_cli/compare/v2.0.2...develop
[2.0.2]: https://github.com/papaemmelab/isabl_cli/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/papaemmelab/isabl_cli/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/papaemmelab/isabl_cli/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/papaemmelab/isabl_cli/releases/tag/v1.0.0
