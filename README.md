# OligoMetab

Notebook-based analysis of the **metabolic machinery of oligodendrocytes** (ketone-body synthesis and utilisation, lactate / monocarboxylate shuttles, glycolysis and pyruvate fate, TCA cycle and oxidative phosphorylation, fatty-acid β-oxidation, lipid and cholesterol synthesis and transport, nutrient-sensing regulators) across the datasets assembled in [OligoC4b](https://github.com/christoffermattssonlangseth/OligoC4b): in-house Xenium (mouse AD time course, mouse EAE), Visium (aging mouse brain), the Falcão et al. 2018 EAE scRNA-seq, and sixteen public datasets (mouse AD, aging and demyelination, spatial AD, human AD and MS). The questions and statistics are the ones OligoC4b asked about the complement system, asked of metabolism, with the C4b⁺ disease-associated oligodendrocyte state of OligoC4b as the anchor.

## Layout

```
notebooks/
  src/        percent-format sources (edit these; convert with scripts/nb_from_py.py)
  analysis/   executed notebooks with outputs
scripts/
  oligometab.py             gene panel (183 genes, 21 pathway groups) and all shared helpers
  check_panel_symbols.py    which panel genes resolve in a given .h5ad (no matrix loading)
  nb_from_py.py             percent-format .py -> .ipynb
  run_notebooks.sh          execute the notebooks in place (kernel `sc`)
  check_notebooks.py        notebook hygiene (local + CI)
docs/
  gene_panel.md             the panel, pathway groups, score sets, what is measurable where
  notebooks.md              plain-language guide to both notebooks and the statistics used
  metabolism_findings.md    what the analyses showed
  public_datasets.md        the datasets (re-used from OligoC4b) and what changes for metabolism
results/                    summary tables written by the notebooks (CSV)
```

## Data

No data are built here. The notebooks read the objects produced by OligoC4b:

| Variable | Object |
| --- | --- |
| `OLIGOMETAB_XENIUM_AD_H5AD` | `Xenium_AD_mouse.h5ad` (347-gene panel; only Apoe, Apod, Acsbg1 of the metabolic panel are on it) |
| `OLIGOMETAB_XENIUM_EAE_H5AD` | Xenium EAE 5K object (raw counts, lesion-distance annotation; 100 of 183 panel genes) |
| `OLIGOMETAB_VISIUM_AGING_H5AD` | `visum_aging_brain.h5ad` |
| `OLIGOMETAB_FALCAO_H5AD` | `falcao_et_al_2018.h5ad` |
| `OLIGOMETAB_PUBLIC_PROCESSED_DIR` | folder with the sixteen harmonised public `.h5ad` files (falls back to `OLIGOC4B_PUBLIC_PROCESSED_DIR`) |

Copy `.env.example` to `.env` and fill in the paths; missing files make the corresponding notebook section skip itself.

## Quick start

```bash
conda env create -f environment.yml && conda activate oligometab
python -m ipykernel install --user --name sc --display-name "sc"      # the notebooks are bound to a kernel named `sc`
cp .env.example .env                                                  # fill in paths
python scripts/nb_from_py.py notebooks/src/analysis_spatial_metabolism.py notebooks/analysis/analysis_spatial_metabolism.ipynb sc
sh scripts/run_notebooks.sh                                           # executes both notebooks in place, logs in logs/
python scripts/check_notebooks.py                                     # before committing
```

## Notebooks

| Notebook | Focus |
| --- | --- |
| `analysis/analysis_spatial_metabolism.ipynb` | Panel content, cell-type expression, disease / age / lesion-distance effects, relation to the C4b⁺ state and spatial neighbourhoods in Xenium AD, Xenium EAE, Visium aging and Falcão |
| `analysis/analysis_public_datasets_metabolism.ipynb` | The same questions across the sixteen public datasets, including human MS and AD by lesion type / Braak stage |

Plain-language summaries and the statistics used are in [`docs/notebooks.md`](docs/notebooks.md); the gene panel in [`docs/gene_panel.md`](docs/gene_panel.md).

## Findings in brief

See [`docs/metabolism_findings.md`](docs/metabolism_findings.md).
