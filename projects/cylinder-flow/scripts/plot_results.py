"""Replot owned numerical histories without changing archived figures."""
from pathlib import Path
import csv
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
profile = json.loads((root / "profiles/benchmark.json").read_text(encoding="utf-8"))
if profile["ranks"] == 1:
    folders = sorted((root / "results/raw").glob("formal-*"))
    summaries = [json.loads((p / "summary.json").read_text(encoding="utf-8")) for p in folders]
    with (folders[0] / "history.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    axes[0].semilogy([int(r["iteration"]) for r in rows], [float(r["scaled_max"]) for r in rows], color="#0b746b")
    axes[0].set(xlabel="Iteration", ylabel="Scaled residual", title="Cylinder convergence")
    axes[1].bar(range(1, 10), [s["iteration_seconds"] for s in summaries], color="#0b746b")
    axes[1].set(xlabel="Formal sample", ylabel="Integration time (s)", title="Nine archived single-core runs")
else:
    with (root / "results/raw/history.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    steps = [int(r["iteration"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    axes[0].plot(steps, [float(r["fy"]) for r in rows], color="#0b746b")
    axes[0].axvspan(207, 306, alpha=.15, color="#d77735", label="First confirmed window")
    axes[0].legend()
    axes[0].set(xlabel="Iteration", ylabel="Fy (N)", title="M6 archived force history")
    axes[1].semilogy(steps, [float(r["scaled_flow"]) for r in rows], color="#0b746b")
    axes[1].set(xlabel="Iteration", ylabel="Scaled flow residual", title="M6 convergence")
for ax in axes:
    ax.grid(alpha=.2)
target = root / "data/generated/summary.png"
target.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(target, dpi=180)
plt.close(fig)
print(target)
