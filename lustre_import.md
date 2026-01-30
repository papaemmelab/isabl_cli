# Plan: Shared Input Imports from GCP to Lustre

## Summary

Refactor the input import system so that GCP inputs are copied to a **shared location** on Lustre (`/shared_inputs/{hash}/`) instead of per-analysis. Uses **file-based reference counting** with `flock` locking to prevent duplicate imports and safe cleanup. Import strategy (`rsync` vs `lustre-import`) is configured **per-application**. Import command is **prepended** to `head_job.sh`; cleanup is in **both success and failure paths** (configurable via `lustre_cleanup_on_failure` setting).

---

## Files to Modify

| File | Change |
|------|--------|
| [settings.py](isabl_cli/settings.py) | Add `shared_inputs_path`, `lustre_cleanup_on_failure` to `GCP_CONFIGURATION`; register 2 new commands |
| [lustre_inputs.py](isabl_cli/lustre_inputs.py) | Refactor to shared paths (`/shared_inputs/{md5}/`), add `get_cleanup_command()`, update `get_import_command()` |
| [gcp_lustre.py](isabl_cli/gcp_lustre.py) | Add `run_shared_import()`, `run_shared_cleanup()`, `_rsync_import()` with flock locking |
| [commands.py](isabl_cli/commands.py) | Add `lustre-shared-import` and `lustre-shared-cleanup` CLI commands |
| [app.py](isabl_cli/app.py) | Add `input_import_strategy` attribute; thread `LustreInputs` to `write_command_script`; prepend import, append cleanup in both success/failure paths |
| [tests/test_gcp_lustre.py](tests/test_gcp_lustre.py) | Tests for shared import, cleanup, ref counting, locking |

---

## Lustre Filesystem Layout

```
/scratch/shared_inputs/
  {md5_of_gcs_dir}/              # data directory
    file1.bam
    .import_complete             # marker: import done
  {md5_of_gcs_dir}.refs/         # ref tracking (sibling dir)
    123                          # analysis PK 123 uses this
    456                          # analysis PK 456 uses this
  {md5_of_gcs_dir}.lock          # flock file
```

---

## Implementation Steps

### Step 1: `settings.py`

Add to `GCP_CONFIGURATION` defaults:
```python
"shared_inputs_path": "/shared_inputs",    # base dir on lustre for shared inputs
"lustre_cleanup_on_failure": True,         # release refs even when pipeline fails
```

Add to `SYSTEM_COMMANDS` list:
```python
"isabl_cli.commands.lustre_shared_import",
"isabl_cli.commands.lustre_shared_cleanup",
```

### Step 2: `lustre_inputs.py` — Refactor to Shared Paths

**Change `input_dir`** from `/{analysis_pk}/inputs` to use `shared_inputs_path` from config.

**Change path resolution:** `add()` and `get()` return `/scratch/shared_inputs/{md5_of_gcs_parent_dir}/{filename}` instead of `/scratch/{pk}/inputs/{subdir}/{filename}`.

**Hash function:** Full MD5 hex digest of the GCS parent directory **URI only** (instant, no file I/O). Example:
```python
# Input GCS URI: gs://bucket/genomes/hg38/
# MD5 = hashlib.md5("gs://bucket/genomes/hg38/".encode()).hexdigest()
# Result: "a1b2c3d4e5f6g7h8i9j0k1l2m3n4o5p6"
# Lustre path: /scratch/shared_inputs/a1b2c3d4e5f6g7h8i9j0k1l2m3n4o5p6/

# Two files from same directory share the hash:
lustre.add("gs://bucket/genomes/hg38/chr1.fa")    # → /scratch/shared_inputs/a1b2.../chr1.fa
lustre.add("gs://bucket/genomes/hg38/chr2.fa")    # → /scratch/shared_inputs/a1b2.../chr2.fa

# Different directory gets different hash:
lustre.add("gs://bucket/genomes/hg37/chr1.fa")    # → /scratch/shared_inputs/z9y8x7w6.../chr1.fa
```

**New method `get_cleanup_command()`:** Returns `isabl lustre-shared-cleanup --analysis-pk {pk} --specs '{json}'`.

**Update `get_import_command(strategy)`:** Returns `isabl lustre-shared-import --analysis-pk {pk} --strategy {strategy} --specs '{json}'`. Takes strategy parameter.

### Step 3: `gcp_lustre.py` — Shared Import/Cleanup with Locking

**`run_shared_import(analysis_pk, import_specs, strategy="lustre-import")`:**
For each `(gcs_dir, lustre_dir)` in specs:
1. Create parent dirs, refs dir
2. Acquire exclusive `flock` on `{hash}.lock`
3. Check `.import_complete` marker — if exists, skip import
4. If not exists: import via strategy (`_rsync_import` or `_lustre_import_data`)
5. Write `.import_complete` marker
6. Create ref file `{hash}.refs/{analysis_pk}`
7. Release lock

**`run_shared_cleanup(analysis_pk, import_specs)`:**
For each `(gcs_dir, lustre_dir)`:
1. Acquire exclusive `flock` on `{hash}.lock`
2. Remove ref file `{hash}.refs/{analysis_pk}`
3. If refs dir empty → `shutil.rmtree` data dir, refs dir, lock file
4. Release lock

**`_rsync_import(gcs_path, local_dest)`:** Runs `gsutil -m rsync -r {gcs_path} {local_dest}`.

**`_lustre_import_data(gcp_config, gcs_path, lustre_path)`:** Wraps existing `initiate_import` + `wait_for_imports` for a single directory.

### Step 4: `commands.py` — New CLI Commands

**`lustre-shared-import`:**
- Options: `--analysis-pk`, `--strategy` (choice: lustre-import/rsync), `--specs` (JSON)
- Calls `gcp_lustre.run_shared_import()`

**`lustre-shared-cleanup`:**
- Options: `--analysis-pk`, `--specs` (JSON)
- Calls `gcp_lustre.run_shared_cleanup()`
- Catches exceptions and warns (non-fatal) so cleanup failures don't crash the job

### Step 5: `app.py` — Integration

**New class attribute** on `AbstractApplication`:
```python
input_import_strategy = "lustre-import"  # or "rsync"
```

**In `run_analyses()`** (after `get_command()` call, ~line 1082):
- Retrieve `self._lustre_inputs` (set by app developer in `get_command()`)
- Pass to `write_command_script(analysis, command, lustre_inputs=lustre_inputs)`
- Reset `self._lustre_inputs = None`

**In `write_command_script()`:**
- Accept optional `lustre_inputs` parameter
- Generate import command (if lustre_inputs has entries)
- Generate input cleanup command
- Refactor command chain building to parts-based approach:

```
SUCCESS PATH:
  import -> started -> pipeline -> export -> copy -> finished -> input_cleanup -> scratch_cleanup

FAILURE PATH:
  failed_status && input_cleanup (if lustre_cleanup_on_failure=True)
```

The failure block changes from just `{ failed }` to:
```bash
{ failed && input_cleanup }   # when cleanup_on_failure=True
{ failed }                     # when cleanup_on_failure=False
```

### Step 6: Tests

- `test_shared_import_first_analysis`: verify data imported + ref created + marker written
- `test_shared_import_skip_existing`: verify second analysis skips import, only adds ref
- `test_shared_cleanup_with_remaining_refs`: remove one ref, data persists
- `test_shared_cleanup_last_ref`: remove last ref, data + refs + lock deleted
- `test_lustre_inputs_shared_paths`: verify path resolution uses `/shared_inputs/{md5}/`
- `test_lustre_inputs_commands`: verify import and cleanup command generation
- `test_write_command_script_with_imports`: verify head_job.sh contains import prepended, cleanup in both paths

---

## Generated `head_job.sh` Structure

```bash
{

    umask g+wrx && date && cd /scratch/analyses/00/01/123 && export TMP=... &&
    ( isabl lustre-shared-import --analysis-pk 123 --strategy lustre-import --specs '[...]' ) &&
    isabl patch-status --key 123 --status STARTED &&
    my_pipeline --bam /scratch/shared_inputs/a1b2.../file.bam &&
    ( isabl lustre-export ... ) &&
    ( rsync ... ) &&
    isabl patch-status --key 123 --status SUCCEEDED &&
    ( isabl lustre-shared-cleanup --analysis-pk 123 --specs '[...]' ) &&
    ( rm -rf /scratch/analyses/00/01/123 ) &&
    date

} || {

    isabl patch-status --key 123 --status FAILED &&
    ( isabl lustre-shared-cleanup --analysis-pk 123 --specs '[...]' )

}
```

---

## Locking Details

- **Granularity:** Per shared-input directory (one lock per GCS source dir)
- **Mechanism:** `fcntl.flock(fd, LOCK_EX)` — POSIX, supported on Lustre
- **Scope:** Held only during check+import or check+delete (not during pipeline)
- **Race handling:** Two jobs for same data → first imports, second skips. Both add refs atomically.

---

## App Developer Usage (unchanged API)

```python
class MyApp(AbstractApplication):
    input_import_strategy = "rsync"  # optional, default is "lustre-import"

    def get_command(self, analysis, inputs, settings):
        from isabl_cli.lustre_inputs import LustreInputs
        lustre = LustreInputs(analysis, settings)
        bam_path = lustre.add(inputs["bam_file"])
        self._lustre_inputs = lustre
        return f"my_pipeline --bam {bam_path}"
```

---

## Verification

1. **Unit tests:** Run `pytest tests/test_gcp_lustre.py -v`
2. **Integration check:** Create a mock application with `LustreInputs`, call `write_command_script`, verify the generated script contains import/cleanup commands in correct order
3. **Concurrent test:** Use threading to simulate two analyses importing the same data — verify only one import occurs and both refs exist
4. **Cleanup test:** Remove refs one by one, verify data deleted only after last ref removed
