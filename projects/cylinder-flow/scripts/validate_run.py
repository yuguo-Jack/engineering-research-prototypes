"""Compare new run files with this project's archived numerical anchors; never compare speed."""
from pathlib import Path
import argparse
import csv
import json
import math

ROOT = Path(__file__).resolve().parents[1]
RELATIVE_TOLERANCE = 1e-5
ABSOLUTE_TOLERANCE = 1e-7


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    number = float(value)
    require(math.isfinite(number), "Non-finite numerical value")
    return number


def close(actual, expected, name):
    actual, expected = finite(actual), finite(expected)
    require(math.isclose(actual, expected, rel_tol=RELATIVE_TOLERANCE, abs_tol=ABSOLUTE_TOLERANCE),
            f"Numerical regression at {name}: {actual} vs archived {expected}")


def validate_directory(project, output):
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    case = json.loads((project / "profiles/benchmark.json").read_text(encoding="utf-8"))
    metrics = json.loads((project / "results/metrics.json").read_text(encoding="utf-8"))
    with (output / "history.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    iterations = 66 if case["ranks"] == 1 else 400
    require(summary["iterations"] == iterations and len(rows) == iterations, "Incomplete or different iteration count")
    require([int(r["iteration"]) for r in rows] == list(range(1, iterations + 1)), "History must start at the original initial field")
    require(summary["order"] == 2 and summary["cfl"] == 50, "Different spatial order or CFL")
    record = output / "run-record.json"
    if record.is_file():
        require(json.loads(record.read_text(encoding="utf-8"))["returncode"] == 0, "Recorded solver process failed")
    details = {}
    if case["ranks"] == 1:
        require(summary["converged"] is True, "Solver did not converge")
        require(0 <= finite(summary["residual_scaled_max"]) <= 1e-9, "Residual above the configured criterion")
        close(summary["cd"], metrics["Cd"], "final Cd")
        for row in rows:
            for key in ("res_mass", "res_xmom", "res_ymom", "res_energy", "scaled_max"):
                finite(row[key])
        details = {"iterations": iterations, "Cd": summary["cd"], "residual_scaled_max": summary["residual_scaled_max"]}
    else:
        require(summary["ranks"] == 32 and summary["cells"] == 1485672, "Different partition count or mesh")
        require(summary["restart_iteration_offset"] == 0 and summary["updates_this_run"] == 400, "Restarted or partial run")
        force_history = []
        first_confirmation = None
        for row in rows:
            for key in ("fx", "fy", "fz", "scaled_flow", "res_rho", "res_mx", "res_my", "res_mz", "res_energy", "res_k", "res_omega"):
                finite(row[key])
            force_history.append(float(row["fy"]))
            if len(force_history) >= 100:
                window = force_history[-100:]
                ratio = (max(window) - min(window)) / max(abs(sum(window) / 100), 1e-30)
                require(math.isclose(ratio, finite(row["fy_window_ratio"]), rel_tol=1e-8, abs_tol=1e-12), "Saved window ratio does not match force history")
                if first_confirmation is None and ratio <= .01:
                    first_confirmation = int(row["iteration"])
        require(first_confirmation == metrics["window_confirmation"], "First 100-step force window moved; review this platform/run")
        require(summary["first_100_step_fy_window_confirmation"] == first_confirmation and
                summary["first_100_step_fy_window_start"] == first_confirmation - 99, "Summary window differs from history")
        for anchor in metrics["accuracy_rows"]:
            row = rows[int(anchor["iteration"]) - 1]
            for key in ("fx", "fy", "fz"):
                close(row[key], anchor[key + "_actual_N"], key + " at iteration " + anchor["iteration"])
        for key in ("fx", "fy", "fz"):
            close(summary[key], rows[-1][key], "summary " + key)
        details = {"iterations": iterations, "window_start": first_confirmation - 99, "window_confirmation": first_confirmation}
    return {"status": "passed", "scope": "run_file_numerical_regression_against_archived_prototype", "relative_tolerance": RELATIVE_TOLERANCE,
            "absolute_tolerance": ABSOLUTE_TOLERANCE, "checks_speed": False, "checks_complete_pressure_field": False,
            "official_acceptance_claimed": False, **details}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Directory containing the new summary.json and history.csv")
    args = parser.parse_args()
    try:
        result = validate_directory(ROOT, args.output.resolve())
    except (ValueError, KeyError, FileNotFoundError) as error:
        raise SystemExit(str(error)) from error
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
