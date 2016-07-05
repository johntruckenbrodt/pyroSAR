
import subprocess as sp


# custom error in case a command fails
class InstallationError(Exception):
    def __init__(self, command, errormessage):
        Exception.__init__(self, "\n\nexecuted command:\n{}\n\n{}".format(" ".join(command), errormessage))


# function to execute shell commands
def execute(cmd, cwd=None, env=None):
    proc = sp.Popen(cmd, stdin=None, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd, env=env)
    out, err = proc.communicate()
    if proc.returncode and err:
        if cmd[0] == "make":
            rollback(cwd)
        raise InstallationError(cmd, err)


# function to clean up temporary files created during compilation
def rollback(directory):
    try:
        sp.check_call(["make", "clean"], cwd=directory)
        sp.check_call(["make", "distclean"], cwd=directory)
    except sp.CalledProcessError:
        pass
