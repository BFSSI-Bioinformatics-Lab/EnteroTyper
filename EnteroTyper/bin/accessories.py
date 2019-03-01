import shutil
import logging
from subprocess import Popen, PIPE


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


def dependency_check(dependency: str) -> bool:
    """
    Checks if a given program is present in the user's $PATH
    :param dependency: String of program name
    :return: True if program is in $PATH, False if not
    """
    check = shutil.which(dependency)
    if check is not None:
        return True
    else:
        return False


def check_all_dependencies(dependencies: [str]) -> bool:
    # Dependency check
    logging.info("Conducting dependency check...")
    dependency_dict = dict()
    dependencies_met = True
    for dependency in dependencies:
        dependency_dict[dependency] = dependency_check(dependency)
    if False in dependency_dict.values():
        dependencies_met = False
        logging.error("ERROR: Cannot locate some dependencies in $PATH...")
        for key, value in dependency_dict.items():
            if not value:
                logging.error(f"Dependency missing: {key}")
        return dependencies_met
    else:
        for key, value in dependency_dict.items():
            logging.debug(f"Dependency {key}: {value}")
        logging.info("Dependencies OK")
        return dependencies_met
