#!/usr/bin/env python3
"""Build the Qt application with conservative CPU and memory limits (Linux)."""

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess

MIB = 1024 * 1024
ROOT = Path(__file__).resolve().parents[1]


def memory_available():
    """Available physical memory, further bounded by visible cgroup v2 limits."""
    budgets = []
    try:
        fields = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        budgets.append(int(fields["MemAvailable"].split()[0]) * 1024)
    except (OSError, ValueError, KeyError):
        pass
    # Check ancestor limits too: a child can say 'max' while its parent is capped.
    cgroup_root = Path("/sys/fs/cgroup")
    paths = {cgroup_root}
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            if line.startswith("0::"):
                current = cgroup_root / line[3:].lstrip("/")
                if ".." not in current.parts:
                    paths.add(current)
                    paths.update(p for p in current.parents if cgroup_root in (p, *p.parents))
    except OSError:
        pass
    for path in paths:
        try:
            limit = (path / "memory.max").read_text().strip()
            if limit != "max":
                used = int((path / "memory.current").read_text())
                budgets.append(max(0, int(limit) - used))
        except (OSError, ValueError):
            pass
    return min(budgets) if budgets else None


def positive_int(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", choices=("debug", "release"), default="debug")
    parser.add_argument("--build-dir", type=Path, help="default: build/safe-<config>")
    parser.add_argument("--qmake", default=os.environ.get("QMAKE"), help="Qt 6 qmake executable")
    parser.add_argument("--jobs", type=positive_int, default=4, help="upper limit; memory may reduce it (default: 4)")
    parser.add_argument("--memory-per-job-mib", type=positive_int, default=2560,
                        help="compiler memory allowance (default: 2560 MiB)")
    parser.add_argument("--reserve-mib", type=positive_int, default=2048,
                        help="leave this much available RAM unused (default: 2048 MiB)")
    parser.add_argument("--qmake-arg", action="append", default=[], help="extra qmake argument, e.g. CONFIG+=no_pch")
    parser.add_argument("--dry-run", action="store_true", help="print commands without creating files or building")
    args = parser.parse_args()
    qmake = args.qmake or shutil.which("qmake6") or shutil.which("qmake")
    if not qmake:
        parser.error("Qt 6 qmake not found; pass --qmake /path/to/Qt/6.x/gcc_64/bin/qmake")
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 1
    available = memory_available()
    if available is None:
        jobs = 1
        print("Available RAM unknown; using one compiler job.", flush=True)
    else:
        budget = available - args.reserve_mib * MIB
        jobs = min(args.jobs, cpus, max(0, budget // (args.memory_per_job_mib * MIB)))
        print(f"Available RAM: {available // MIB} MiB; reserve: {args.reserve_mib} MiB; "
              f"allowance per job: {args.memory_per_job_mib} MiB.", flush=True)
        if jobs < 1:
            parser.error("insufficient available RAM for one compiler plus reserve; close other applications and retry")
    build_dir = (args.build_dir or ROOT / "build" / f"safe-{args.config}").resolve()
    configure = [qmake, str(ROOT / "generation-terrain-procedural.pro"),
                 f"CONFIG+={args.config}", *args.qmake_arg]
    build = ["make", f"-j{jobs}"]
    print(f"Build directory: {build_dir}\n{shlex.join(configure)}\n{shlex.join(build)}", flush=True)
    if args.dry_run:
        return
    version = subprocess.check_output([qmake, "-query", "QT_VERSION"], text=True).strip()
    if not version.startswith("6."):
        parser.error(f"Qt 6 required; selected qmake reports Qt {version}")
    build_dir.mkdir(parents=True, exist_ok=True)
    # Parent make invocations and environment flags must not undo the job limit.
    env = os.environ.copy()
    for key in ("MAKEFLAGS", "MFLAGS", "GNUMAKEFLAGS"):
        env.pop(key, None)
    subprocess.run(configure, cwd=build_dir, env=env, check=True)
    subprocess.run(build, cwd=build_dir, env=env, check=True)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        raise SystemExit(error.returncode)
