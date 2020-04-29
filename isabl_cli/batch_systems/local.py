import subprocess
import os

from click import progressbar

from isabl_cli import api
from isabl_cli import utils


def submit_local(app, command_tuples):
    """Submit analyses locally and serially."""
    # make sure analyses are in submitted status to avoid
    # merging project level analyses in every success
    ret = []
    api.patch_analyses_status([i for i, _ in command_tuples], "SUBMITTED")
    label = f"Running {len(command_tuples)} analyses..."

    with progressbar(command_tuples, label=label) as bar:
        for i, j in bar:
            log = app.get_command_log_path(i)
            err = app.get_command_err_path(i)
            api.patch_analysis_status(i, "STARTED")
            oldmask = os.umask(0o22)
            status = app._get_after_completion_status(i)

            with open(log, "w") as stdout, open(err, "w") as stderr:
                try:
                    subprocess.check_call(
                        args=j,
                        cwd=i["storage_url"],
                        stdout=stdout,
                        stderr=stderr,
                        shell=True,
                    )

                except subprocess.CalledProcessError:
                    status = "FAILED"

            os.umask(oldmask)
            i = api.patch_analysis_status(i, status)
            ret.append((i, i["status"]))

    return ret
