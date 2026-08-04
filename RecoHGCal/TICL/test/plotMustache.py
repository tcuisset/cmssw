#!/usr/bin/env python3
"""Plot Mustache and dynamic-dPhi boundaries from MustacheESProductDumper."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def metadata(path):
    result = {}
    with path.open() as source:
        for line in source:
            if not line.startswith("#"):
                break
            key, value = line[1:].strip().split("=", 1)
            result[key] = float(value)
    return result


parser = argparse.ArgumentParser()
parser.add_argument("scan", type=Path, nargs="?", default=Path("mustache_scan.csv"))
parser.add_argument("-o", "--output", type=Path, default=Path("mustache_eta_phi.png"))
parser.add_argument("--show", action="store_true")
args = parser.parse_args()

meta = metadata(args.scan)
data = np.genfromtxt(args.scan, delimiter=",", names=True, comments="#", skip_header=len(meta))
etas = np.unique(data["eta"])
phis = np.unique(data["phi"])
shape = (len(etas), len(phis))

fig, ax = plt.subplots(figsize=(8, 6))
styles = (
    ("in_mustache", "tab:blue", "-", "Mustache"),
    ("in_dynamic_dphi", "tab:orange", "--", "dynamic $\\Delta\\phi$"),
    ("in_combined", "black", "-", "combined"),
)
for field, color, linestyle, _ in styles:
    accepted = data[field].reshape(shape)
    ax.contour(etas, phis, accepted.T, levels=[0.5], colors=[color], linestyles=[linestyle], linewidths=2)

ax.plot(meta["seed_eta"], meta["seed_phi"], marker="*", color="crimson", markersize=12)
if meta["cluster_energy"] > 0:
    energy_label = f"candidate E = {meta['cluster_energy']:g} GeV"
else:
    energy_label = f"candidate $E_T$ = {meta['cluster_et']:g} GeV"
ax.set(
    xlabel=r"candidate $\eta$",
    ylabel=r"candidate $\phi$ [rad]",
    title=rf"Mustache acceptance: seed $\eta={meta['seed_eta']:g}$, $\phi={meta['seed_phi']:g}$; {energy_label}",
)
ax.grid(alpha=0.25)
handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=2, label=label)
           for _, color, linestyle, label in styles]
handles.append(Line2D([0], [0], marker="*", color="crimson", linestyle="", markersize=10, label="seed"))
ax.legend(handles=handles)
fig.tight_layout()
fig.savefig(args.output, dpi=160)
print(f"Wrote {args.output}")
if args.show:
    plt.show()
