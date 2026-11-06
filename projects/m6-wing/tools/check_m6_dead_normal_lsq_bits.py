"""Compare actual Solver::gradients live outputs in two compile-time variants."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def compare(left, right, name):
    if left != right:
        first = next((i for i, (a, b) in enumerate(zip(left, right)) if a != b), min(len(left), len(right)))
        raise AssertionError(f"{name}: bit mismatch at byte {first}; lengths={len(left)},{len(right)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("off", type=Path)
    parser.add_argument("on", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--all-scratch-equal", action="store_true",
                        help="Also compare dead scalar scratch, for original-source versus default-OFF regression")
    args = parser.parse_args()
    off, on, destination = args.off, args.on, args.destination
    destination.mkdir(parents=True, exist_ok=True)
    records = []
    # The caller supplies the build directory; no C-drive temporary files.
    with tempfile.TemporaryDirectory(prefix="gradient-pair-", dir=destination) as temporary:
        work = Path(temporary)
        for label, program in (("off", off), ("on", on)):
            result = subprocess.run([str(program.resolve()), str(work / label)], capture_output=True, text=True)
            (destination / f"{label}.stdout.log").write_text(result.stdout, encoding="utf-8")
            (destination / f"{label}.stderr.log").write_text(result.stderr, encoding="utf-8")
            if result.returncode:
                raise RuntimeError(f"{label} actual gradient probe failed: {result.stdout}\n{result.stderr}")
        left_manifest = (work / "off/manifest.json").read_bytes()
        compare(left_manifest, (work / "on/manifest.json").read_bytes(), "case identity")
        manifest = json.loads(left_manifest)
        assert manifest["nonlinear_updates"] == 0 and manifest["ranks"] == 1
        assert len(manifest["cases"]) == 22
        for case in manifest["cases"]:
            name = case["name"]
            live_off = (work / f"off/{name}.live.bin").read_bytes()
            live_on = (work / f"on/{name}.live.bin").read_bytes()
            compare(live_off, live_on, name + " live arrays")
            scratch_off = (work / f"off/{name}.scratch.bin").read_bytes()
            scratch_on = (work / f"on/{name}.scratch.bin").read_bytes()
            same_scratch = scratch_off == scratch_on
            if args.all_scratch_equal or not case["optimized_eligible"]:
                compare(scratch_off, scratch_on, name + " fallback scratch")
            records.append(dict(name=name, optimized_eligible=case["optimized_eligible"],
                                live_bytes=len(live_off), live_sha256=hashlib.sha256(live_off).hexdigest(),
                                dead_scratch_equal=same_scratch))
        # This proves the ON target exercised the optimized branch, rather than
        # accidentally compiling both probes with the original seven-variable loop.
        if not args.all_scratch_equal:
            assert any(x["optimized_eligible"] and not x["dead_scratch_equal"] for x in records)
        # A deliberate live-data error must be rejected by this comparison code.
        negative_control_rejected = False
        try:
            compare(b"actual-live-state", b"actual-live-statf", "deliberate bit error")
        except AssertionError:
            negative_control_rejected = True
        assert negative_control_rejected
    report = dict(status="actual_small_mesh_gradient_live_arrays_bit_identical",
                  nonlinear_updates=0, cases=records, cases_passed=len(records),
                  live_bytes_compared=sum(x["live_bytes"] for x in records),
                  scratch_difference_cases=sum(not x["dead_scratch_equal"] for x in records),
                  original_source_default_off_regression=args.all_scratch_equal,
                  negative_control_rejected=negative_control_rejected,
                  limitations="One-rank synthetic 27-cell mesh; no halo peers. Real original/own209 32-rank gates remain required.")
    (destination / "result.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
