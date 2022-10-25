import subprocess


def run_command(command, verbose=False):
    """Run a shell command"""
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True
    )
    logout = ""
    while True:
        output = process.stdout.readline().strip()
        output = output.decode("utf-8")
        if output == "" and process.poll() is not None:
            break
        if output:
            logout += str(output.strip()) + "\n"
            if verbose:
                print((str(output.strip())))
    rc = process.poll()
    return logout