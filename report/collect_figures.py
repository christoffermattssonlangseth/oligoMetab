#!/usr/bin/env python3
"""Pull selected figures out of the executed notebooks into report/figures/ as PNG files.

Each entry maps an output file name to (notebook path, substring that identifies the code cell, index of the image
within that cell's outputs; negative indices count from the end). Run from the repo root after the notebooks have
been executed:

    python report/collect_figures.py
"""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

SPATIAL = ROOT / "notebooks/analysis/analysis_spatial_metabolism.ipynb"
PUBLIC = ROOT / "notebooks/analysis/analysis_public_datasets_metabolism.ipynb"

FIGURES = {
    # in-house spatial notebook
    "eae_celltype_scores.png": (SPATIAL, 'heat(s_ct, "Xenium EAE: mean pathway score by cell type', -1),
    "eae_vs_control_pathways.png": (SPATIAL, 'heat(d, "Xenium EAE: pathway-score difference EAE − CONTROL per cell type', 0),
    "eae_lesion_distance_lines.png": (SPATIAL, 'g = sns.relplot(data=ld, x="bin", y="score", hue="population"', 0),
    "eae_lesion_distance_cholesterol.png": (SPATIAL, 'heat(m.sub(far, axis=0), f"Xenium EAE, {pw}', 6),
    "eae_da_vs_homeostatic.png": (SPATIAL, 'fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))', 0),
    "eae_neighbourhood_ratio.png": (SPATIAL, 'sns.boxplot(data=long, x="target", y="ratio"', 0),
    "eae_lesion_matched.png": (SPATIAL, 'sns.pointplot(data=pb_plot, x="bin", y="median_ratio", hue="target"', 0),
    "eae_spatial_peak.png": (SPATIAL, 'spatial_panel(ad_eae, "sample_name", s, ["C4b", "Slc16a1", "Ldha", "Hk2", "Hcar2", "Cpt1a", "Plin2"]', 0),
    "visium_age_scores_wm.png": (SPATIAL, 'g = sns.catplot(data=tidy, x="age_group", y="score", col="pathway", col_wrap=6, kind="strip"', -1),
    "falcao_cluster_scores.png": (SPATIAL, 'heat(s_cl, "Falcão: mean pathway score by cluster', -1),
    # public datasets notebook
    "public_celltype_z.png": (PUBLIC, 'heat(z.T, "Pathway score by cell type, z-scored within dataset', 0),
    "public_oligo_pathway_diff.png": (PUBLIC, 'heat(d, f"{ct}: pathway-score difference (disease / aged / demyelinated vs reference)"', 0),
    "public_microglia_pathway_diff.png": (PUBLIC, 'heat(d, f"{ct}: pathway-score difference (disease / aged / demyelinated vs reference)"', 2),
    "public_astrocyte_pathway_diff.png": (PUBLIC, 'heat(d, f"{ct}: pathway-score difference (disease / aged / demyelinated vs reference)"', 3),
    "public_oligo_genes_energy.png": (PUBLIC, 'heat(d, f"{ct}: pseudobulk log2 fold change, {label} genes"', 0),
    "public_oligo_genes_lipid.png": (PUBLIC, 'heat(d, f"{ct}: pseudobulk log2 fold change, {label} genes"', 1),
    "public_c4b_pathway_enrichment.png": (PUBLIC, 'heat(E, "Pathway-level enrichment among C4b-correlated genes in oligodendrocytes', 0),
    "public_c4b_high_vs_neg.png": (PUBLIC, 'heat(D, "Pathway score in C4b-high minus C4b-negative oligodendrocytes', 0),
}


def images_in_cell(cell: dict) -> list:
    out = []
    for o in cell.get("outputs", []):
        data = o.get("data", {})
        if "image/png" in data:
            out.append(base64.b64decode(data["image/png"]))
    return out


def main() -> int:
    missing = []
    for name, (nb_path, needle, idx) in FIGURES.items():
        if not nb_path.exists():
            missing.append((name, "notebook missing")); continue
        nb = json.loads(nb_path.read_text(encoding="utf-8"))
        found = False
        for cell in nb["cells"]:
            if cell["cell_type"] != "code" or needle not in "".join(cell["source"]):
                continue
            imgs = images_in_cell(cell)
            if not imgs:
                continue
            try:
                (OUT / name).write_bytes(imgs[idx])
                found = True
                break
            except IndexError:
                pass
        if not found:
            missing.append((name, "no image"))
        else:
            print(f"wrote {name}")
    for name, why in missing:
        print(f"skipped {name}: {why}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
