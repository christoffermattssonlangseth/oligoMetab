#!/usr/bin/env python3
"""Check which panel genes resolve in each dataset without loading the matrices (reads var/obs with h5py only).

    python scripts/check_panel_symbols.py path1.h5ad path2.h5ad ...

Prints, per file: species guess, n_obs, obs columns, the panel genes that are missing, and the alias that was used
where the primary symbol was absent. Used once when the panel was designed; kept for re-checks when datasets change.
"""
from __future__ import annotations

import os
import sys

import h5py
import pandas as pd
from anndata.experimental import read_elem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import oligometab as om  # noqa: E402


class _Fake:
    """Minimal stand-in for AnnData so oligometab.resolve() can be reused on var/obs read with h5py."""

    def __init__(self, var_names, obs):
        self.var_names = pd.Index(var_names)
        self.obs = obs


def main(paths):
    summary = []
    for p in paths:
        with h5py.File(p, "r") as f:
            var = read_elem(f["var"])
            obs = read_elem(f["obs"])
        var_names = var.index.astype(str)
        species = str(obs["species"].iloc[0]) if "species" in obs else ("human" if var_names.str.isupper().mean() > 0.8 else "mouse")
        if "species" not in obs:
            obs = obs.assign(species=species)
        fake = _Fake(var_names, obs)
        r = om.resolve(fake, om.ALL_GENES + om.CONTEXT)
        miss = [g for g in om.ALL_GENES if r[g] is None and om.sym(species, g) is not None]
        alias = {g: v for g, v in r.items() if v is not None and v != om.sym(species, g)}
        print(f"\n=== {os.path.basename(p)}  ({species}, {obs.shape[0]} obs, {len(var_names)} vars)")
        print("obs columns:", list(obs.columns)[:40])
        for c in ["cell_type_coarse", "cellType", "cell_type", "group", "condition", "age_group", "model", "sample", "sample_name", "run"]:
            if c in obs and obs[c].nunique() < 60:
                print(f"  {c}: {sorted(map(str, obs[c].unique()))[:60]}")
        print(f"panel genes present: {len(om.ALL_GENES) - len(miss)} / {len(om.ALL_GENES)}")
        print("missing:", miss)
        print("resolved via alias:", alias)
        summary.append({"file": os.path.basename(p), "species": species, "n_obs": obs.shape[0], "n_present": len(om.ALL_GENES) - len(miss),
                        "missing": ",".join(miss)})
    print("\n" + pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:])
