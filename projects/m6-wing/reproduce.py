"""Build, test and reproduce this project independently from the repository root."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build"
CASE = json.loads((ROOT / "profiles/benchmark.json").read_text(encoding="utf-8"))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def build_environment():
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    scratch = BUILD / "tmp"
    scratch.mkdir(parents=True, exist_ok=True)
    for name in ("TEMP", "TMP", "TMPDIR"):
        env[name] = str(scratch)
    if os.name == "nt" and not env.get("INCLUDE"):
        installer = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) / "Microsoft Visual Studio/Installer/vswhere.exe"
        if installer.is_file():
            where = subprocess.run([str(installer), "-latest", "-products", "*", "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64", "-property", "installationPath"], capture_output=True, text=True, check=True)
            install = where.stdout.strip()
            batch = Path(install) / "Common7/Tools/VsDevCmd.bat"
            if install and batch.is_file():
                # Capture only into this process; never write or print the environment.
                command = f'call "{batch}" -no_logo -arch=x64 -host_arch=x64 >nul && set'
                result = subprocess.run(["cmd.exe", "/d", "/s", "/c", command], capture_output=True, text=True, check=True)
                for line in result.stdout.splitlines():
                    if "=" in line and not line.startswith("="):
                        key, value = line.split("=", 1)
                        env[key] = value
    mpi = shutil.which("mpiexec", path=env.get("PATH"))
    if os.name == "nt" and mpi:
        candidates = [Path(mpi).resolve().parents[2] / "sdk", Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) / "Microsoft SDKs/MPI"]
        for sdk in candidates:
            if (sdk / "Include/mpi.h").is_file():
                env.setdefault("MSMPI_INC", str(sdk / "Include"))
                env.setdefault("MSMPI_LIB64", str(sdk / "Lib/x64"))
                break
    return env


def execute(command, label, env):
    BUILD.mkdir(parents=True, exist_ok=True)
    with (BUILD / (label + ".log")).open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"{label} failed ({result.returncode}); see build/{label}.log")


def executable():
    name = "cylinder" if CASE["ranks"] == 1 else "m6"
    suffix = ".exe" if os.name == "nt" else ""
    for path in (BUILD / (name + suffix), BUILD / "Release" / (name + suffix)):
        if path.is_file():
            return path
    raise FileNotFoundError("No built solver. Run: python reproduce.py build")


def preserve_probe_output():
    """The archived probe deliberately refuses old output; retain it before a new test."""
    previous = BUILD / "skip-postflow-assembly-probe"
    if not previous.exists():
        return
    history = BUILD / "test-history"
    for path in (BUILD, previous, history):
        if path.is_symlink() or getattr(path, "is_junction", lambda: False)():
            raise ValueError("Test output must not use filesystem links")
        if not path.resolve().is_relative_to(ROOT):
            raise ValueError("Test output escapes this project")
    if not previous.is_dir():
        raise ValueError("Probe output path is not a directory; preserve and inspect it")
    history.mkdir(exist_ok=True)
    destination = history / ("skip-postflow-assembly-" + str(time.time_ns()))
    if not destination.resolve().is_relative_to(ROOT):
        raise ValueError("Test history escapes this project")
    previous.rename(destination)
    log = BUILD / "ctest.log"
    if log.is_file() and not log.is_symlink():
        shutil.copyfile(log, destination / "previous-ctest.log")


def verify_results():
    metrics = json.loads((ROOT / "results/metrics.json").read_text(encoding="utf-8"))
    provenance = json.loads((ROOT / "results/PROVENANCE.json").read_text(encoding="utf-8"))
    for record in provenance["files"]:
        path = (ROOT / record["path"]).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError("Invalid result path")
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError("Archived result changed: " + record["path"])
    if CASE["ranks"] == 1:
        summaries = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results/raw").glob("formal-*/summary.json"))]
        if len(summaries) != 9 or any(s["iterations"] != 66 or not s["converged"] for s in summaries):
            raise ValueError("Expected nine converged 66-step samples")
        median = statistics.median(s["iteration_seconds"] for s in summaries)
        if not math.isclose(median, metrics["iteration_seconds_median"], abs_tol=1e-12, rel_tol=0):
            raise ValueError("Median timing does not match archived summaries")
        value = {"samples": 9, "iterations": 66, "median_iteration_seconds": median}
    else:
        with (ROOT / "results/raw/history.csv").open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != 400 or [int(r["iteration"]) for r in rows] != list(range(1, 401)):
            raise ValueError("Expected complete original-initial 400-step history")
        first = next(r for r in rows if int(r["iteration"]) >= 100 and float(r["fy_window_ratio"]) <= .01)
        confirmation = int(first["iteration"])
        if confirmation != 306 or confirmation - 99 != metrics["main_window_start"]:
            raise ValueError("Stability milestone differs")
        for step, key in ((207, "start_iteration_seconds"), (306, "confirmation_iteration_seconds")):
            if not math.isclose(float(rows[step-1]["seconds"]), metrics[key], abs_tol=1e-8, rel_tol=0):
                raise ValueError("Timing milestone differs: " + str(step))
        summary = json.loads((ROOT / "results/raw/summary.json").read_text(encoding="utf-8"))
        if not math.isclose(summary["iteration_seconds"], metrics["full400_iteration_seconds"], abs_tol=1e-8, rel_tol=0):
            raise ValueError("Full integration timing differs from final summary")
        value = {"rows": 400, "window_start": 207, "window_confirmation": 306, "ratio": float(first["fy_window_ratio"])}
    print(json.dumps({"status": "passed", "scope": "archived_result_integrity_and_accounting", **value}, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("doctor", "build", "test", "verify-results", "check-run", "import-data", "run", "plot"))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--native", action="store_true", help="GCC/Clang host-specific optimization; rebuild for every target")
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--input", type=Path, help="Default: data/external (M6) or its cylinder.bin")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.action == "verify-results":
        verify_results()
        return
    if args.action == "check-run":
        if args.output is None:
            parser.error("check-run requires --output PATH_TO_NEW_RUN")
        subprocess.run([sys.executable, "-B", "-X", "utf8", str(ROOT / "scripts/validate_run.py"), str(args.output.resolve())], check=True)
        return
    if args.action == "import-data":
        if args.archive is None:
            parser.error("import-data requires --archive PATH_TO_AUTHORIZED_ARCHIVE")
        subprocess.run([sys.executable, "-B", "-X", "utf8", str(ROOT / "scripts/import_data.py"), str(args.archive.resolve())], check=True)
        return
    if args.action == "plot":
        subprocess.run([sys.executable, "-B", "-X", "utf8", str(ROOT / "scripts/plot_results.py")], check=True)
        return
    env = build_environment()
    if args.action == "doctor":
        print(json.dumps({"platform": platform.platform(), "python": platform.python_version(), "cmake": shutil.which("cmake", path=env.get("PATH")), "compiler": shutil.which("cl" if os.name == "nt" else "c++", path=env.get("PATH")), "mpi_launcher": shutil.which("mpiexec", path=env.get("PATH")), "mpi_required": CASE["ranks"] > 1, "external_data_present": (ROOT / "data/external").exists()}, ensure_ascii=False, indent=2))
    elif args.action == "build":
        if args.jobs < 1:
            parser.error("--jobs must be positive")
        # Explicitly reset this cached variable when returning from --native.
        command = ["cmake", "-S", str(ROOT), "-B", str(BUILD), "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_CXX_FLAGS=", "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"]
        if not (BUILD / "CMakeCache.txt").exists() and shutil.which("ninja", path=env.get("PATH")):
            command += ["-G", "Ninja"]
        config = ROOT / "profiles/build.json"
        if config.is_file():
            command += json.loads(config.read_text(encoding="utf-8"))["definitions"]
        if args.native:
            if os.name == "nt":
                parser.error("--native is for GCC/Clang; use the portable MSVC build")
            command += ["-DCMAKE_CXX_FLAGS=-march=native -mtune=generic -ffp-contract=off"]
        execute(command, "configure", env)
        execute(["cmake", "--build", str(BUILD), "--config", "Release", "--parallel", str(args.jobs)], "compile", env)
        print(json.dumps({"status": "built", "binary": str(executable()), "native": args.native}, ensure_ascii=False))
    elif args.action == "test":
        preserve_probe_output()
        execute(["ctest", "--test-dir", str(BUILD), "-C", "Release", "--output-on-failure"], "ctest", env)
        print(json.dumps({"status": "passed", "scope": "numerical_core_tests", "log": "build/ctest.log"}))
    elif args.action == "run":
        data_dir = ROOT / "data/external"
        input_path = (args.input or (data_dir / "cylinder.bin" if CASE["ranks"] == 1 else data_dir)).resolve()
        manifest = json.loads((ROOT / "data/manifest.json").read_text(encoding="utf-8"))
        for row in manifest["files"]:
            path = input_path if CASE["ranks"] == 1 else input_path / row["path"]
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                    digest.update(block)
            if path.stat().st_size != row["bytes"] or digest.hexdigest() != row["sha256"]:
                raise ValueError("Original initial input mismatch: " + row["path"])
        output = (args.output or ROOT / "runs" / time.strftime("%Y%m%d-%H%M%S")).resolve()
        if output.exists() and any(output.iterdir()):
            raise FileExistsError("Output is not empty; choose a fresh directory")
        command = [str(executable()), "--input", str(input_path), "--output", str(output), *CASE["options"]]
        if CASE["ranks"] > 1:
            command = ["mpiexec", "-n", str(CASE["ranks"]), *command]
        if args.dry_run:
            print(json.dumps({"command": command, "creates_output": False}, ensure_ascii=False, indent=2))
            return
        output.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        with (output / "stdout.log").open("w", encoding="utf-8") as stdout, (output / "stderr.log").open("w", encoding="utf-8") as stderr:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
        report = {"returncode": result.returncode, "command": command, "process_seconds": time.perf_counter()-started, "timing_scope": "current_local_run_not_archived_formal_measurement"}
        write_json(output / "run-record.json", report)
        if result.returncode:
            raise RuntimeError("Solver failed; see output logs")
        print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, FileNotFoundError, FileExistsError) as error:
        raise SystemExit(str(error)) from error
