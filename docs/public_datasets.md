# Public datasets

OligoMetab re-uses, unchanged, the sixteen public datasets that OligoC4b downloaded, processed and annotated (`OligoC4b/scripts/oligoc4b_public.py`, `OligoC4b/notebooks/build/build_public_datasets.ipynb`). Nothing is rebuilt here: the notebooks read the harmonised `.h5ad` files from `OLIGOMETAB_PUBLIC_PROCESSED_DIR` (falls back to `OLIGOC4B_PUBLIC_PROCESSED_DIR`). Accessions, loaders, format notes and the cell-type annotation method are documented in OligoC4b (`docs/public_datasets.md`, `docs/cell_type_annotation.md`); the short version:

| Dataset | GEO | Species, modality | Comparison |
| --- | --- | --- | --- |
| Park et al. 2023, AD hippocampus | GSE224398 | mouse scRNA-seq | App^NL-G-F^ vs control, 1/3/6 mo |
| Aging snRNA-seq, hippocampus + caudate putamen | GSE212576 | mouse snRNA-seq | old vs young |
| Ximerakis et al. 2019, whole brain | GSE129788 | mouse scRNA-seq | old vs young |
| Kaya et al. 2022, aged white vs grey matter | GSE202579 | mouse scRNA-seq | WM vs GM at 24 mo |
| Zhou et al. 2020, 5XFAD | GSE140511 | mouse snRNA-seq | 5XFAD vs non-Tg, 7 and 15 mo |
| Chen et al. 2020, Spatial Transcriptomics | GSE152506 | mouse ST | App^NL-G-F^ vs WT, 3–18 mo |
| LPC + cuprizone corpus callosum | GSE293850 | mouse snRNA-seq | LPC / cuprizone vs controls |
| Serpina3n cKO + cuprizone | GSE319903 | mouse snRNA-seq | cuprizone vs normal diet |
| Jäkel et al. 2019, MS white matter | GSE118257 | human snRNA-seq | MS lesion types vs control |
| Absinta et al. 2021, chronic active MS | GSE180759 | human snRNA-seq | lesion edge / core / periplaque vs control |
| Schirmer et al. 2019, MS lesions | UCSC Cell Browser `ms` | human snRNA-seq | lesion types vs control |
| Lerma-Martin et al. 2024, subcortical MS | GSE279180 / GSE279181 | human snRNA-seq + Visium | chronic active / inactive vs control |
| Senescent-like glia in MS, 2025 | GSE277435 | human Visium | MS subtypes (no controls) |
| Leng et al. 2021, SFG + entorhinal cortex | GSE147528 | human snRNA-seq | Braak 6 vs Braak 0 |
| Sadick et al. 2022, PFC astrocyte/oligodendrocyte-enriched | GSE167494 | human snRNA-seq | AD vs non-symptomatic |

Every object has raw counts in `layers["counts"]`, log-normalised `X` and the harmonised `obs` columns `dataset`, `species`, `modality`, `sample`, `group`, `group_ref`, `cell_type_coarse`, `cell_type_original`, `condition_original`.

## What changes for metabolism compared with complement

- All panel genes are quantifiable in human data (the C4A/C4B multi-mapping problem does not apply), so the human MS and AD datasets contribute real disease contrasts, not just detection tables.
- The C4b-anchored analyses (co-expression, genome-wide rank, C4b-high vs C4b-negative) are still mouse-only, for the same reason as in OligoC4b.
- The caveats on the marker-based coarse annotation, on ambient microglial RNA in demyelinating-lesion nuclei (LPC, cuprizone), on the Kaya comparison being WM vs GM in aged brain, and on underpowered pseudobulk tests (Park, human sets with 1–2 control donors) carry over unchanged.
