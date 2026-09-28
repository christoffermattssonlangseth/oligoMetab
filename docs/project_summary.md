# OligoMetab: project summary

*One page that says what was asked, what was done, what came out, and where everything is. Written 2026-09-28. The full numbers live in `docs/metabolism_findings.md`, the notebooks and `results/`.*

## The ask

Emulate the OligoC4b analysis (complement / *C4b*⁺ oligodendrocytes across in-house spatial data and sixteen public datasets) for **metabolism**: ketone bodies, glycolysis and lactate, TCA / OXPHOS, fatty acids, lipids and their regulators. The analysis is about metabolism **across all cell types and conditions**; the *C4b*⁺ disease-associated oligodendrocyte state is one axis among several, not the anchor.

## What was built

| Piece | Where |
| --- | --- |
| Gene panel: 183 mouse genes in 21 pathway groups, human symbol mapping, aliases; 12 pathway score sets | `scripts/oligometab.py`, documented in `docs/gene_panel.md` |
| Shared helpers: detection, pseudobulk tests, pathway scores, *C4b* co-expression, genome-wide Spearman ranks (chunked), pathway rank enrichment, spatial neighbourhood tests | `scripts/oligometab.py` |
| In-house notebook: Xenium AD, Xenium EAE (5K, 107 samples), Visium aging, Falcão sorted cells | `notebooks/analysis/analysis_spatial_metabolism.ipynb` (source `notebooks/src/`) |
| Public notebook: the sixteen OligoC4b datasets, every coarse cell type tested | `notebooks/analysis/analysis_public_datasets_metabolism.ipynb` |
| Summary tables (50 CSVs) | `results/` |
| Findings, notebook guide, dataset notes | `docs/metabolism_findings.md`, `docs/notebooks.md`, `docs/public_datasets.md` |
| Report with figures and PDF | `report/report.md`, `report/OligoMetab_metabolism_report.pdf` |
| Tooling: percent-format → ipynb converter, run script, output dumper, symbol checker, hygiene check + CI | `scripts/` |

Both notebooks were executed on the group's analysis machine (kernel `sc`, scanpy 1.9.8) and end with a written interpretation cell. Data are the OligoC4b objects; nothing was rebuilt.

## What is measurable where

- Xenium mouse AD (347 genes): 3 of 183 panel genes (*Apoe*, *Apod*, *Acsbg1*). Metabolism is not measurable there.
- Xenium mouse EAE (5K): 100 of 183; no ketogenesis / ketolysis enzymes, 2 of 15 OXPHOS, 3 of 12 cholesterol-synthesis genes.
- Visium, Falcão and all public data: 177–183 of 183. Unlike *C4A*/C4B, the metabolic genes are quantifiable in human data.

## Findings

1. **Division of labour is stable across platforms and species.** Neurons and astrocytes carry the highest glycolysis / TCA / OXPHOS scores; oligodendrocytes the highest lipid-, cholesterol- and myelin-lipid synthesis (newly formed oligodendrocytes the highest cholesterol synthesis); microglia and macrophages the highest lipid transport / storage and pentose phosphate and the lowest cholesterol synthesis; OPC / COP / NFOL the highest ketone utilisation of the lineage. *Hcar2* (ketone / niacin receptor) is microglial; *Slc16a1* (MCT1) is in ~40 % of oligodendrocytes in situ.
2. **EAE is a tissue-wide metabolic shift.** In 92 EAE vs 15 control samples, cholesterol synthesis falls and lipid transport / storage rises in essentially every cell type, TCA falls in most, glycolysis rises only in microglia and DA oligodendrocytes. Toward lesions, glycolysis and TCA fall in neurons, oligodendrocytes and vascular cells, lipid storage rises almost everywhere, and cholesterol synthesis falls most in OPCs (−0.33), newly formed oligodendrocytes (−0.29) and DA oligodendrocytes (−0.16).
3. **The disease-associated oligodendrocyte state** (paired within 107 samples; confirmed in sorted Falcão cells): less cholesterol and myelin-lipid synthesis (*Hmgcs1*, *Fdps*, *Cyp51*, *Dhcr7*, *Sqle* down), more lipid uptake and storage (*Plin4* 5×, *Plin2* 3×, *Abca1* 2.5×, *Lpl*, *Cd36*), a *Hk2* / *Pdk4* / *Ldha* glycolytic-switch signature, *Bdh1* induction, less MCT1. The glycolysis / TCA *scores* point opposite ways in Xenium and Falcão because of gene-set composition; the gene tables are the safer read-out.
4. **Ketone bodies are not part of the oligodendrocyte program.** The enzymes belong to astrocytes and OPCs, *Hmgcs2* is barely expressed, *Oxct1* / *Acat1* are flat or lower; only *Bdh1* rises in EAE-associated oligodendrocytes, and it falls in 5 of 14 public contrasts.
5. **The *C4b* axis.** The one metabolic feature that follows *C4b* in every mouse dataset is lipid handling: *Plin4* (rank 2 of 3,673 panel genes in EAE DA oligodendrocytes), *Apod* (ranks 8–42 in five public datasets), *Abca1*; *C4b*-high oligodendrocytes have a higher lipid-storage score in all seven public mouse datasets. Energy metabolism is unchanged or slightly anti-correlated; cholesterol synthesis is anti-correlated in sorted cells.
6. **Spatial partners.** *Hcar2*⁺ cells (1.7×), *Hcar2*⁺ myeloid cells (1.8×), *Slc16a3*⁺ (MCT4) astrocytes (2.3×), *Hk2*⁺, *Pdk1*⁺ and *Plin2*⁺ cells are enriched within 30 µm of *C4b*-high oligodendrocytes in EAE, surviving lesion-distance matching and restriction to non-lesion tissue (as C5aR1⁺ cells did in OligoC4b).
7. **Public data.** Pathway shifts are small (mostly within ±0.05) but the same genes recur in oligodendrocytes (*Lpl* and *Hk2* up in 10 of 14 contrasts; *Apoe*, *Plin2*, *Fabp5* up; *Slc16a1*, *Gal3st1*, *Acsl6*, *Ppara*, *Scd1* down). Microglia give the most reproducible signature of all (*Lpl* up in 12 of 13 contrasts, median log2FC 2.8; *Apoe*, *Plin2*, *Pparg*, *Cd36*, *Hcar2*). Astrocytes lose ketone utilisation in every MS dataset. Human MS oligodendrocytes have, if anything, higher cholesterol synthesis and more *Hif1a*, *Txnip* and *Ddit4*.
8. **Aging (Visium)** is a mild version of the EAE pattern in white matter: less lipid / cholesterol synthesis (*Fa2h*, *Ugt8a*, *Cyp51*, *Hmgcs1* down 10–19 %), more *Apoe* and *Plin2*, less *Ddit4* and *Slc2a1*.

## Caveats

Xenium AD: one section per genotype × age and three genes. EAE 5K panel: ketone, most OXPHOS and cholesterol-synthesis genes absent. In situ transcripts can come from a neighbouring cell (*Hcar2*, *Cd36*, *Slc16a3* in "oligodendrocytes"); sorted cells are the control, and they are the replicates there. Pathway scores depend on set composition. Pseudobulk tests need ≥3 samples per group; several public contrasts fall short. Nuclei detect metabolic transcripts far less than whole cells. Microglial ambient RNA inflates *Apoe* / *Lpl* / *Hcar2* in lesion nuclei. Human *C4B* is not quantifiable, so the *C4b* axis is mouse-only.

## Suggested next steps

Lipid-droplet (PLIN2 / PLIN4) and cholesterol-synthesis readouts in *C4b*⁺ oligodendrocytes; MCT1 protein and lactate supply near lesions; human *C4A*/C4B re-quantification so the *C4b* axis can be tested in MS; flux measurements on sorted DA vs homeostatic oligodendrocytes to settle the glycolysis / TCA direction.

## How to re-run

```bash
cp .env.example .env            # paths to the OligoC4b objects
python scripts/nb_from_py.py notebooks/src/analysis_spatial_metabolism.py notebooks/analysis/analysis_spatial_metabolism.ipynb sc
python scripts/nb_from_py.py notebooks/src/analysis_public_datasets_metabolism.py notebooks/analysis/analysis_public_datasets_metabolism.ipynb sc
sh scripts/run_notebooks.sh     # sequential; ~10 min each on the analysis machine
python report/collect_figures.py && python report/build_report.py
```

Process notes: the two notebooks must run one at a time on the shared 48 GB machine; SciPy's Mann–Whitney defaults to an exact method for small gene sets, which cost hours until forced to the asymptotic method (`oligometab._mwu_p`); `nbconvert --inplace` writes only at the end, so progress is followed through `results/*.csv` timestamps.
