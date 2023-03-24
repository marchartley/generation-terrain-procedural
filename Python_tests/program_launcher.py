from typing import Optional

import psutil
import time
import datetime
import os

def checkIfProcessRunning(processNames) -> Optional[psutil.Process]:
    """
    Check if there is any running process that contains the given name processName.
    """
    processNames = [s.lower() for s in processNames]
    # Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            for processName in processNames:
                if processName in proc.name().lower():
                    if proc.status() == 'zombie':
                        proc.kill()
                    else:
                        return proc
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return None


def adios_process(cmd: str):
    if cmd:
        print(f"Launching {cmd} again!")
        os.system(" ".join(cmd))
    return


if __name__ == "__main__":
    main = psutil.Process(os.getpid())
    # main.nice(psutil.IDLE_PRIORITY_CLASS)
    launchedMyself = False
    process: Optional[psutil.Process] = None

    checkingName = "interface"
    try:
        while True:
            process = checkIfProcessRunning(["interface"])
            print(process)
            if process:
                process_cmd = process.as_dict()["cmdline"]
                psutil.wait_procs([process], callback=lambda p: adios_process(process_cmd))
            else:
                os.system("/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/build/interface")
                launchedMyself = True
                time.sleep(5)
    except KeyboardInterrupt:
        if launchedMyself and process is not None:
            process.kill()
        raise

