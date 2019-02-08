import logging
from pathlib import Path
from subprocess import Popen, PIPE
from __init__ import __version__, __author__, __email__


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    logging.info(f"Version: {__version__}")
    logging.info(f"Author: {__author__}")
    logging.info(f"Email: {__email__}")
    quit()


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


def run_subprocess(cmd: str, get_stdout: bool = False) -> str:
    if get_stdout:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out = out.decode().strip()
        err = err.decode().strip()
        if out != "":
            return out
        elif err != "":
            return err
        else:
            return ""
    else:
        p = Popen(cmd, shell=True)
        p.wait()
