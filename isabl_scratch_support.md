# Scratch Storage Implementation Plan

## Overview

Add support for running analyses on fast Lustre scratch storage (`/scratch`) while maintaining final results on gcsfuse-mounted datalake. Analyses will run entirely in scratch, then rsync/export results to the final location upon completion.

## Architecture

### Current Flow
- Analysis dir: `/isabl/data/analysis/01/01/1` (gcsfuse, slow)
- Pipeline runs directly in that location
- Results stay there permanently

### New Flow
1. **Creation**: Dual directories created
   - Final: `/isabl/data/analysis/01/01/1` (gcsfuse)
   - Scratch: `/scratch/isabl/data/analysis/01/01/1` (Lustre)
   - Symlinks in final dir → scratch for log files
2. **Execution**: Pipeline runs in scratch directory
3. **Completion**: rsync or lustre-export from scratch → final
4. **Cleanup**: Remove scratch directory after successful copy

## Configuration

### System-Level Settings

Add to [settings.py](isabl_cli/settings.py):

```python
"SCRATCH_STORAGE_DIRECTORY": None,  # e.g., "/scratch" to enable
```

- Default `None` = backward compatible (scratch disabled)
- Add to `_PATH_STRINGS` set for automatic path expansion
- **Note:** `SCRATCH_USE_LUSTRE_EXPORT` has been removed - now configured per-application

### Application-Level Settings

Each application can configure its scratch copy strategy:

```python
class MyApplication(AbstractApplication):
    # Scratch copy strategy: "rsync" (default) or "lustre-export"
    scratch_copy_strategy = "rsync"  # or "lustre-export"
```

- **"rsync"** (default): Uses rsync for copying (works anywhere)
- **"lustre-export"**: Uses GCP Lustre export (faster, requires GCP_CONFIGURATION.lustre_export_enabled=True)

## Critical Files & Changes

### 1. [data.py](isabl_cli/data.py) - Storage Path Generation

**Add two helper functions** (after line 258):

```python
def get_scratch_storage_url(endpoint, identifier, use_hash=False):
    """Get scratch storage path if configured, else None."""
    scratch_root = system_settings.SCRATCH_STORAGE_DIRECTORY
    if not scratch_root:
        return None
    return system_settings.MAKE_STORAGE_DIRECTORY(
        root=scratch_root, base=endpoint, identifier=identifier, use_hash=use_hash
    )

def get_final_storage_url(endpoint, identifier, use_hash=False):
    """Get final (permanent) storage path."""
    return system_settings.MAKE_STORAGE_DIRECTORY(
        root=system_settings.BASE_STORAGE_DIRECTORY,
        base=endpoint, identifier=identifier, use_hash=use_hash
    )
```

**Update `get_storage_url`** (line 250) to call `get_final_storage_url` for backward compatibility.

---

### 2. [app.py](isabl_cli/app.py) - Main Application Logic

#### A. Application-Level Attribute

**Add class attribute** (after line 105, near `gcp_lustre_export`):

```python
# Scratch copy strategy configuration. Set to one of:
# - "rsync" (default): Use rsync to copy from scratch to final storage
# - "lustre-export": Use GCP Lustre export for high-performance copying
# Only relevant when SCRATCH_STORAGE_DIRECTORY is configured at system level.
scratch_copy_strategy = "rsync"
```

**Rationale**: Follows the same pattern as `gcp_lustre_export` - per-application configuration.

#### B. Directory Creation & Symlinks

**Modify `_patch_analysis`** (lines 1759-1769):

```python
def _patch_analysis(self, analysis):
    # Get both storage URLs
    final_url = data.get_final_storage_url("analyses", analysis["pk"], use_hash=True)
    scratch_url = data.get_scratch_storage_url("analyses", analysis["pk"], use_hash=True)

    analysis["storage_url"] = final_url
    analysis["_scratch_storage_url"] = scratch_url

    # Create directories
    utils.makedirs(final_url)
    if scratch_url:
        utils.makedirs(scratch_url)
        # Create symlinks for logs: final → scratch
        for filename in ["head_job.sh", "head_job.log", "head_job.err"]:
            final_path = join(final_url, filename)
            scratch_path = join(scratch_url, filename)
            if os.path.exists(final_path) or os.path.islink(final_path):
                os.remove(final_path)
            os.symlink(scratch_path, final_path)

    return api.patch_instance(
        "analyses", analysis["pk"],
        results=self._get_analysis_results(analysis, created=True),
        storage_url=analysis["storage_url"],
    )
```

**Rationale**: Creates both dirs upfront, symlinks enable real-time log monitoring.

#### C. Execution Directory Helpers

**Add helper methods** (before line 1117):

```python
def _get_execution_directory(self, analysis):
    """Return scratch dir if configured, else final dir."""
    scratch_url = analysis.get("_scratch_storage_url")
    if scratch_url and os.path.isdir(scratch_url):
        return scratch_url
    return analysis["storage_url"]

def _should_use_lustre_export_for_scratch(self):
    """Determine if lustre-export should be used for scratch copy.

    Returns True if:
    1. Application requests lustre-export via scratch_copy_strategy, AND
    2. System allows lustre-export (GCP_CONFIGURATION.lustre_export_enabled=True)

    Follows the same pattern as _should_export_to_gcs().
    """
    # Check application preference
    strategy = getattr(self, "scratch_copy_strategy", "rsync")
    if strategy != "lustre-export":
        return False

    # Check if system allows lustre-export
    gcp_config = getattr(system_settings, "GCP_CONFIGURATION", None) or {}
    return gcp_config.get("lustre_export_enabled", False)

def _get_scratch_to_final_copy_command(self, analysis):
    """Generate rsync or lustre-export command based on app preference."""
    scratch_url = analysis.get("_scratch_storage_url")
    if not scratch_url:
        return ""

    final_url = analysis["storage_url"]

    # Use application-level setting with system validation
    use_lustre_export = self._should_use_lustre_export_for_scratch()

    if use_lustre_export:
        return f"isabl lustre-export --lustre-path {scratch_url} --analysis-pk {analysis['pk']} --no-delete-after"
    else:
        return (
            f"mkdir -p {final_url} && "
            f"rsync -av --delete --partial --partial-dir=.rsync-partial {scratch_url}/ {final_url}/"
        )

def _get_scratch_cleanup_command(self, analysis):
    """Generate cleanup command."""
    scratch_url = analysis.get("_scratch_storage_url")
    return f"rm -rf {scratch_url}" if scratch_url else ""
```

**Rationale**:
- `_should_use_lustre_export_for_scratch()` follows the same pattern as `_should_export_to_gcs()` (lines 1271-1281)
- Application chooses strategy, system validates it's available
- Falls back to rsync if system doesn't support lustre-export

#### D. Command Script Generation

**Modify `write_command_script`** (lines 1117-1182):

Change `outdir` determination (line 1128):
```python
outdir = self._get_execution_directory(analysis)  # scratch if available
```

Add copy/cleanup commands (after line 1150):
```python
copy_command = self._get_scratch_to_final_copy_command(analysis)
cleanup_command = self._get_scratch_cleanup_command(analysis)
```

Update command chains (lines 1152-1178) to add `copy_command` and `cleanup_command` after `export_command` and before `finished`:

**Example with all features**:
```python
if import_command and export_command and copy_command:
    command = (
        f"umask g+wrx && date && cd {outdir} && {tmpdir} && "
        f"( {import_command} ) && {started} && {command} && "
        f"( {export_command} ) && ( {copy_command} ) && {finished} && "
        f"( {cleanup_command} ) && date"
    )
```

**Order**: `import → started → pipeline → export → copy → finished → cleanup`

#### E. Path Conversion for get_command

**Modify `run_analyses`** (before line 1029):

```python
# Pass scratch path to application via storage_url
analysis_for_command = dict(i)
scratch_url = i.get("_scratch_storage_url")
if scratch_url:
    analysis_for_command["storage_url"] = scratch_url

command = self.get_command(analysis_for_command, inputs, self.settings)
```

**Rationale**: Applications receive scratch path transparently, no code changes needed.

#### F. Update Path Helper Methods

**Modify these methods** (lines 1105-1115) to use `_get_execution_directory`:

```python
def get_command_script_path(self, analysis):
    return join(self._get_execution_directory(analysis), "head_job.sh")

def get_command_log_path(self, analysis):
    return join(self._get_execution_directory(analysis), "head_job.log")

def get_command_err_path(self, analysis):
    return join(self._get_execution_directory(analysis), "head_job.err")
```

---

### 3. [gcp_lustre.py](isabl_cli/gcp_lustre.py) - Retry Logic

**Add retry helper** (before line 119):

```python
def _retry_with_exponential_backoff(func, max_retries=5, base_delay=10, operation_name="operation"):
    """Retry function with exponential backoff: 10s, 20s, 40s, 80s, 160s."""
    import time

    for attempt in range(1, max_retries + 1):
        try:
            return func()
        except Exception as e:
            if attempt == max_retries:
                raise GCPLustreExportError(
                    f"{operation_name} failed after {max_retries} attempts: {e}"
                )

            delay = base_delay * (2 ** (attempt - 1))
            click.echo(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"{operation_name} attempt {attempt} failed: {e}. "
                f"Retrying in {delay}s..."
            )
            time.sleep(delay)
```

**Modify `initiate_export`** (lines 119-180):

Wrap the subprocess.run call:

```python
def _run_export_command():
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    output = json.loads(result.stdout)
    operation_name = output.get("name")
    if not operation_name:
        raise GCPLustreExportError(f"No operation name: {result.stdout}")
    return operation_name

operation_name = _retry_with_exponential_backoff(
    _run_export_command,
    max_retries=5,
    base_delay=10,
    operation_name="Lustre export initiation"
)
```

**Rationale**: Handles transient GCP API failures, user-specified retry pattern.

---

### 4. [api.py](isabl_cli/api.py) - Permission Finalization

**Modify `_set_analysis_permissions`** (line 673):

Add check at beginning:

```python
def _set_analysis_permissions(analysis):
    """Set permissions on analysis directory.

    If scratch storage used, permissions handled by rsync --chmod.
    """
    scratch_url = analysis.get("_scratch_storage_url")
    if scratch_url:
        click.echo("Skipping permission changes (handled by rsync from scratch)")
        return

    # Existing permission logic...
```

**Rationale**: rsync sets permissions atomically via `--chmod`, no double-handling needed.

---

## Backward Compatibility

- **Default**: `SCRATCH_STORAGE_DIRECTORY=None` → scratch disabled, existing behavior
- **Migration**: Set config → new analyses use scratch, old analyses unchanged
- **Force/Restart**: If old analysis restarted after scratch enabled, automatically add scratch path

**Handling in `run_analyses`** (after line 1001):
```python
if force and i["status"] not in {"SUCCEEDED", "FINISHED", "REJECTED"}:
    system_settings.TRASH_ANALYSIS_STORAGE(i)
    # Add scratch if now configured but not in analysis
    if not i.get("_scratch_storage_url") and system_settings.SCRATCH_STORAGE_DIRECTORY:
        i["_scratch_storage_url"] = data.get_scratch_storage_url(
            "analyses", i["pk"], use_hash=True
        )
    utils.makedirs(i["storage_url"])
    if i.get("_scratch_storage_url"):
        utils.makedirs(i["_scratch_storage_url"])
```

---

## Error Handling

### Copy Failure
- **Behavior**: Analysis marked FAILED, scratch preserved
- **Mechanism**: Existing `{ command } || { failed }` wrapper catches copy errors
- **Recovery**: Admin can manually rsync and update status

### Symlink Issues
- **Prevention**: Symlinks created synchronously before job submission
- **Contingency**: Users can access scratch logs directly

### Disk Space
- **Mitigation**: Immediate cleanup after successful copy
- **Monitoring**: Document scratch sizing requirements

---

## Verification & Testing

### End-to-End Test

1. **Configure scratch** (system settings):
   ```python
   SCRATCH_STORAGE_DIRECTORY = "/scratch"
   GCP_CONFIGURATION = {
       "lustre_export_enabled": True,  # Required if any app uses lustre-export
       # ... other GCP settings
   }
   ```

   **And in your application**:
   ```python
   class TestApp(AbstractApplication):
       scratch_copy_strategy = "rsync"  # or "lustre-export"
   ```

2. **Run analysis**:
   ```bash
   isabl run-analyses <app> --targets <ids> --commit
   ```

3. **Verify during execution**:
   - Check both directories exist:
     ```bash
     ls -la /isabl/data/analysis/01/01/1/
     ls -la /scratch/isabl/data/analysis/01/01/1/
     ```
   - Verify symlinks work:
     ```bash
     tail -f /isabl/data/analysis/01/01/1/head_job.log
     ```
   - Should show real-time output from scratch

4. **Verify after completion**:
   - Check final location has results:
     ```bash
     ls -la /isabl/data/analysis/01/01/1/
     ```
   - Check scratch cleaned up:
     ```bash
     ls -la /scratch/isabl/data/analysis/01/01/1/  # Should not exist
     ```
   - Check status:
     ```bash
     isabl get-metadata --analyses <id>
     ```

5. **Test failure scenario**:
   - Manually kill analysis mid-execution
   - Verify scratch preserved for debugging
   - Verify analysis marked FAILED

6. **Test lustre export option**:
   - Set application's `scratch_copy_strategy = "lustre-export"`
   - Ensure system has `GCP_CONFIGURATION.lustre_export_enabled = True`
   - Run analysis
   - Verify lustre-export command used instead of rsync
   - Check for retry behavior if export fails

7. **Test application-level strategy**:
   - Create two apps: one with `scratch_copy_strategy="rsync"`, another with `"lustre-export"`
   - Run both analyses simultaneously
   - Verify first uses rsync, second uses lustre-export
   - Check logs show correct copy method

### Unit Tests

**Add to `tests/test_scratch_storage.py`**:
- `test_scratch_disabled_backward_compat`
- `test_scratch_enabled_dual_dirs`
- `test_scratch_symlinks_created`
- `test_get_execution_directory`
- `test_scratch_copy_command_rsync`
- `test_scratch_copy_command_lustre`
- `test_should_use_lustre_export_for_scratch`
- `test_application_level_strategy_override`
- `test_cleanup_command`
- `test_path_conversion_for_get_command`
- `test_retry_logic`

---

## Key Design Change: Application-Level Strategy

### Why Application-Level?

Originally, `SCRATCH_USE_LUSTRE_EXPORT` was a system-wide setting that applied to ALL analyses. This was too inflexible because:
- Different pipelines have different performance characteristics
- Some apps might benefit from lustre-export, others work fine with rsync
- No way to test both strategies simultaneously
- Forces all apps to use same copy method

### New Pattern (Following `gcp_lustre_export`)

Just like how applications opt-in to lustre export for their outputs via `gcp_lustre_export=True`, they now choose their scratch copy strategy:

```python
class FastGenomicsPipeline(AbstractApplication):
    # Large outputs, use lustre-export for speed
    scratch_copy_strategy = "lustre-export"

class QuickQCPipeline(AbstractApplication):
    # Small outputs, rsync is fine (default)
    scratch_copy_strategy = "rsync"  # or omit, it's the default
```

### Decision Logic

1. **Application declares preference** via `scratch_copy_strategy`
2. **System validates availability** via `GCP_CONFIGURATION.lustre_export_enabled`
3. **Falls back gracefully** to rsync if lustre-export unavailable

This matches the existing `_should_export_to_gcs()` pattern exactly.

### Removed System Setting

- ❌ **Removed**: `SCRATCH_USE_LUSTRE_EXPORT` from `settings.py`
- ✅ **Replaced with**: Application-level `scratch_copy_strategy` attribute

---

## Implementation Order

### Stage 1: Foundation
1. Remove `SCRATCH_USE_LUSTRE_EXPORT` from [settings.py](isabl_cli/settings.py)
2. Add helpers to [data.py](isabl_cli/data.py)
3. Add `scratch_copy_strategy` attribute to [app.py](isabl_cli/app.py)
4. Update `_patch_analysis` in [app.py](isabl_cli/app.py)
5. **Test**: Directory creation and symlinks

### Stage 2: Execution
6. Add helper methods to [app.py](isabl_cli/app.py) including `_should_use_lustre_export_for_scratch()`
7. Update `write_command_script` in [app.py](isabl_cli/app.py)
8. Update `run_analyses` path conversion
9. **Test**: Full analysis lifecycle with both rsync and lustre-export strategies

### Stage 3: Reliability
10. Add retry logic to [gcp_lustre.py](isabl_cli/gcp_lustre.py)
11. Update `_set_analysis_permissions` in [api.py](isabl_cli/api.py)
12. **Test**: Error handling and retries

### Stage 4: Verification
13. Add unit tests
14. Add integration tests
15. Manual testing checklist

---

## Deployment Plan

1. **Deploy with scratch disabled** (default `None`)
2. **Verify backward compatibility** (run existing test suite)
3. **Enable in staging environment**
4. **Monitor metrics**: copy duration, failure rate, disk usage
5. **Enable in production**

---

## Key Design Decisions

### Why symlinks vs copying logs?
- **Real-time access** during execution (no sync delay)
- **Zero I/O overhead** (no periodic copying)
- **Automatic cleanup** (symlinks removed with parent dir)

### Why store scratch path in memory vs database?
- **Derivable from config** (no schema migration needed)
- **Ephemeral by design** (scratch is temporary)
- **Backward compatible** (no database changes)

### Why pass scratch via `analysis["storage_url"]`?
- **Zero app changes** (apps already use this field)
- **Clean interface** (no new fields to understand)
- **Database consistency** (storage_url in DB = final location)

### Why retry for lustre export?
- **Transient API failures** (GCP best practice)
- **User requirement** (exponential backoff requested)
- **Reliability** (prevents analysis failure on temporary issues)

---

## Open Questions (for clarification if needed)

1. **Rsync permissions**: Use `--chmod=ugo=rwX` (Isabl convention) or `--perms` (preserve exact)?
2. **Lustre export timeout**: Current 3 hours sufficient for large exports?
3. **Scratch structure**: Mirror final path or use flat structure?
4. **Log rotation**: Should there be automatic rotation for long-running analyses?

---

## Risk Mitigation

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Copy failure after successful pipeline | Low | Medium | Retry logic, scratch preserved, clear errors |
| Symlink breakage | Very Low | Low | Synchronous creation, test coverage |
| Disk space exhaustion | Medium | Medium | Immediate cleanup, sizing docs, monitoring |
| Backward compatibility break | Very Low | High | Default disabled, extensive testing, phased rollout |

---

## Summary

This implementation adds scratch storage support to isabl_cli with:
- **Minimal invasiveness**: 5 files, mostly additive changes
- **Backward compatibility**: Disabled by default, existing behavior preserved
- **Transparency**: Applications receive scratch path automatically
- **Reliability**: Retry logic, error handling, cleanup
- **Application-level flexibility**: Each app chooses rsync vs lustre-export strategy
- **Follows existing patterns**: Mirrors `gcp_lustre_export` design
- **User-friendly**: Real-time log access via symlinks
- **Lustre import removed**: Only export functionality retained

### Key Difference from Initial Implementation

**Original Design**: System-wide `SCRATCH_USE_LUSTRE_EXPORT` setting (all apps use same strategy)

**Updated Design**: Application-level `scratch_copy_strategy` attribute (each app chooses independently)

This matches the existing pattern where applications opt-in to features like `gcp_lustre_export`, giving fine-grained control over which pipelines use which copy method.

The design follows existing Lustre import/export patterns and integrates seamlessly with the current architecture.
