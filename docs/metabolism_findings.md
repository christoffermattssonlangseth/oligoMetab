# Metabolism across cell types, disease and the C4b⁺ oligodendrocyte state: what the analyses showed

Question (September 2026): *do for metabolism (ketone bodies, glycolysis, and the rest) what OligoC4b did for complement*: which cell types run which pathway, how that changes with disease, age, demyelination and lesion distance in every annotated cell type, and how it relates to the C4b⁺ disease-associated oligodendrocyte state.

Notebooks: `notebooks/analysis/analysis_spatial_metabolism.ipynb` (Xenium AD, Xenium EAE, Visium aging, Falcão scRNA-seq) and `notebooks/analysis/analysis_public_datasets_metabolism.ipynb` (sixteen public datasets). Each ends with a written interpretation; this page is the condensed version. Gene panel: [`gene_panel.md`](gene_panel.md) (183 genes, 21 pathway groups). Pathway *scores* are `scanpy.tl.score_genes` values (set mean minus a size-matched control set, log scale) compared by difference; single genes are compared by log2 fold change of pseudobulk means.

## What is measurable where

| Platform | Metabolic panel genes present |
| --- | --- |
| Xenium mouse AD (347 genes) | 3 of 183: Apoe, Apod, Acsbg1. Metabolism cannot be read from this dataset. |
| Xenium mouse EAE (5K panel) | 100 of 183. Missing: ketogenesis (Hmgcs2, Hmgcl), Oxct1 / Acat1, 13 of 15 OXPHOS subunits, 9 of 12 cholesterol-synthesis genes, Apoe. Ketone scores are therefore not computable there. |
| Visium aging, Falcão, all public datasets | 177–183 of 183 (whole transcriptome; losses are lowly expressed genes filtered at build time). Unlike C4A/C4B, every metabolic gene is quantifiable in the human data. |

## 1. Division of labour between cell types (in-house)

- **Xenium EAE, 27 cell types.** Astrocytes and neurons carry the highest glycolysis and TCA scores; homeostatic oligodendrocytes the highest lipid-synthesis score; newly formed oligodendrocytes the highest cholesterol-synthesis score; microglia, macrophages and above all foamy myeloid cells the highest lipid-transport / storage score together with the lowest cholesterol synthesis. Slc16a1 (MCT1) is detected in 39 % of oligodendrocytes and 84 % of DA astrocytes; Hcar2 (ketone / niacin receptor) is myeloid (8 % of microglia, 13 % of activated myeloid cells, ≤1.5 % of oligodendrocytes); Bdh1 is astrocytic (51 %) more than oligodendroglial (18 %).
- **Falcão sorted cells.** Along the lineage, OPC/COP/NFOL have the highest ketone-utilisation and TCA scores, mature oligodendrocytes the highest cholesterol- and myelin-lipid synthesis, microglia the highest glycolysis and pentose-phosphate scores and the lowest cholesterol synthesis.

## 2. EAE: a tissue-wide metabolic shift (Xenium EAE, 92 EAE vs 15 control samples, pseudobulk per cell type)

- **Cholesterol synthesis falls and lipid transport / storage rises in essentially every cell type** with enough control samples: oligodendrocytes (−0.17 / +0.14), DA oligodendrocytes (−0.08 / +0.23), OPCs, newly formed oligodendrocytes, astrocytes (−0.22 / +0.13), DA astrocytes (−0.30), microglia (−0.14 / +0.24), activated myeloid cells, neurons, endothelial, vascular, fibroblasts, stromal, ependymal and Schwann cells (−0.33). The **TCA score falls in most** of them (−0.06 to −0.12). Glycolysis rises only in microglia (+0.11) and DA oligodendrocytes (+0.03) and falls in OPCs, neurons and endothelial cells.
- Gene level: Hcar2, Cd36, Hk2, Plin2, Plin4, Abca1, Lpl, Slc16a3 and Txnip rise across many cell types (2–7-fold for Hcar2 and Cd36); Ppargc1a, Acss2, Acsl6, Mlxipl and Ldlr fall in several. Because Xenium is in situ, the myeloid genes among the "oligodendrocyte" hits (Hcar2, Cd36, Slc16a3) can be spill-over from adjacent processes; the sorted-cell check is in section 4.
- **Lesion distance.** Toward lesions, glycolysis and TCA scores fall in neurons, homeostatic oligodendrocytes, endothelial and vascular cells, lipid transport / storage rises in nearly every cell type (largest in endothelial cells, microglia, OPCs and oligodendrocytes), and cholesterol synthesis falls most steeply in OPCs (−0.33), newly formed oligodendrocytes (−0.29) and DA oligodendrocytes (−0.16), while rising in homeostatic oligodendrocytes near lesions (+0.17). Myeloid cells are glycolytic at every distance and become more lactate-producing toward lesions.

## 3. The disease-associated oligodendrocyte state (Xenium EAE, paired within 107 samples)

DA oligodendrocytes versus homeostatic oligodendrocytes of the same sample: glycolysis −0.10, TCA −0.12, lipid synthesis −0.15, cholesterol synthesis −0.16, lipid transport / storage +0.15 (positive in 99 % of samples), β-oxidation unchanged. Genes up: Hk2 (3.5×), Mlxipl, Hcar2, Cd36, Lpl, Slc16a3, Slc2a4, Txnip, Fabp7, Plin2, Abca1, Ldha, Pdk1. Genes down: Slc2a3, Eno2, Ppargc1a, Got2, Mpc2, Sirt3, Prkaa2, Dnm1l, Pcx, Fh1, Mfn2, Elovl1, Ogdh, Pfkp. C4b-high versus C4b-negative lineage cells show the same pattern.

## 4. Sorted cells confirm the core intrinsically (Falcão 2018, EAE-MOL vs control-MOL clusters)

- Cholesterol synthesis −1.0 (Hmgcs1, Fdps, Cyp51, Dhcr7, Dhcr24, Sqle each 20–30 % lower), myelin-lipid synthesis −0.30, lipid synthesis −0.18 (Scd1, Elovl), OXPHOS −0.17, ketone utilisation −0.28 (Acat2, Acss2 lower).
- Lipid transport / storage +0.95 (Plin4 5×, Plin2 3×, Abca1 2.5×), β-oxidation +0.20 (Cpt1a 2.5×), TCA +0.25, glycolysis +0.14 (Ldha 2.2×, Pfkp), Pdk4 4×, Bdh1 2.3×, Slc16a1 (MCT1) −0.24 log2.
- So the shared, intrinsic core of the disease-associated state is **less cholesterol / myelin-lipid synthesis, more lipid uptake and storage, a Pdk4 / Ldha / Hk2-type glycolytic-switch signature, Bdh1 induction and less MCT1**. Where the two datasets disagree (glycolysis / TCA / β-oxidation *scores* down in Xenium DA oligodendrocytes, up in sorted EAE-MOL) the reason is set composition: the 15-gene Xenium glycolysis set is weighted toward Eno2, Pfkp and Slc2a3, which fall, the whole-transcriptome set toward Aldoa, Tpi1, Gapdh and Ldha, which rise. The gene-level tables are the safer read-out.
- Ketone bodies: Bdh1 rises in DA states, Oxct1 / Acat1 are flat or lower, Hmgcs2 is barely expressed in oligodendrocytes; Hcar2 stays a myeloid receptor (2 % of EAE-MOL cells, a handful).

## 5. The C4b axis

- Inside DA oligodendrocytes (Xenium EAE) the metabolic genes that track C4b are lipid-droplet / lipid-handling genes: Plin4 is rank 2 of 3,673 panel genes, Apod 42, Cers2 72, Fasn 86. No metabolic pathway is enriched as a set, and across the whole lineage of EAE animals glycolysis, TCA and mitochondrial-biogenesis genes are significantly *anti*-correlated with C4b.
- In sorted MOL the same holds with more power: lipid transport / storage is the top pathway (Apod ρ = 0.46, Plin4 0.40, Abca1 0.32), cholesterol synthesis is significantly anti-correlated, glycolysis, TCA and β-oxidation mildly positive (Aldoa 0.34, Pdk4 0.30, Cpt1a 0.22, Bdh1 0.22).
- In Visium the spot-level C4b correlates are white-matter lipid genes (Apoe, Apod, Scd2, Fa2h, Ugt8a), i.e. composition. Within white-matter-rich spots, C4b-high spots have higher lipid-synthesis / storage / myelin-lipid and pentose-phosphate scores and lower glycolysis / TCA / OXPHOS in all six sections.

## 6. Spatial partners of C4b-high oligodendrocytes (Xenium EAE, 30 µm, 87 samples)

Hcar2⁺ cells (median ratio 1.7), Hcar2⁺ myeloid cells (1.8), Slc16a3⁺ (MCT4, lactate-exporting) astrocytes (2.3), Hk2⁺ (1.8), Pdk1⁺ (1.3) and Plin2⁺ (1.4) cells are enriched around C4b-high compared with C4b-negative oligodendrocytes (paired Wilcoxon p < 10⁻⁵); Slc2a1⁺ cells are not (0.97). The enrichment survives matching the source oligodendrocytes for lesion distance (1.2–1.75, 76 samples) and restriction to non-lesion tissue (1.3–2.2). As with C5aR1⁺ myeloid cells in OligoC4b, C4b-high oligodendrocytes mark glycolytic, lactate-exporting, lipid-loaded micro-niches also outside lesions. The non-oligodendrocyte neighbourhood of C4b-high cells scores lower on glycolysis / TCA and higher on lipid storage, which mostly reflects its composition (myeloid rather than neuronal / astrocytic).

## 7. Aging (Visium, 6 sections)

Whole-section glycolysis, TCA and OXPHOS scores rise slightly with age (r ≈ 0.84–0.86, p ≈ 0.03–0.04; slopes of a few hundredths over 15 months), driven by housekeeping genes (Pkm, Aldoa, Gapdh, Hk1, Ldhb, Ndufv1, Uqcrc1, Sdhb). In white-matter-rich spots the lipid-synthesis score falls (p = 0.03) and Fa2h, Ugt8a, Cyp51 and Hmgcs1 decline by 10–19 %, while Apoe (+0.12 log2), Plin2 (+0.48 log2) and Bdh1 rise and Ddit4 (−0.63 log2), Slc2a1 (−0.15 log2) and Sirt1 fall. Aging white matter shows a milder version of the EAE pattern: less lipid / cholesterol synthesis, more lipid storage.

## 8. Public datasets (sixteen, about one million cells / nuclei / spots)

- **Same division of labour** as in situ: neurons and astrocytes highest on glycolysis / TCA / OXPHOS scores, oligodendrocytes on myelin-lipid, cholesterol and lipid synthesis, microglia and immune cells on lipid transport / storage and pentose phosphate; ketone utilisation highest in astrocytes and OPCs. Hcar2 is microglial (4–9 % in droplet data); Slc16a1 (MCT1) is in 29–57 % of whole oligodendrocytes but far fewer nuclei (modality effect).
- **Pathway-level disease effects are small** (mostly within ±0.05 in oligodendrocytes, versus ±0.15–0.25 in situ in EAE). The clear ones: Schirmer MS oligodendrocytes more glycolytic (+0.16) and oxidative (+0.11) with more lipid storage (+0.08); aged oligodendrocytes (Ximerakis) lower OXPHOS (−0.06), cholesterol synthesis (−0.07) and myelin-lipid synthesis (−0.14); ketone synthesis / utilisation slightly lower in MS oligodendrocytes in Absinta, Jäkel and Lerma-Martin; Braak-6 oligodendrocytes (Leng) less lipid synthesis (−0.09); cuprizone / LPC oligodendrocytes more lipid storage (+0.13 / +0.27, descriptive).
- **Gene-level effects are consistent.** In oligodendrocytes over 14 contrasts: Lpl up in 10, Hk2 up in 10, Apoe / Plin2 / Fabp5 up in 6–7; Slc16a1 (MCT1) down in 5, Gal3st1 down in 6, Acsl6, Ppara, Lss, Srebf1, Cpt1c and Scd1 down in 4–5. Bdh1 falls in 5 of 14 (unlike its induction in EAE oligodendrocytes), so ketone handling is not unidirectional. **Microglia carry the most reproducible signature: Lpl up in 12 of 13 contrasts (median log2FC 2.8), Apoe and Plin2 in 8, Pparg and Cd36 in 6–8, Hcar2 in 6 of 10**, the lipid-laden disease-associated microglial program (lipid-storage score up in LPC, cuprizone, Jäkel, aged white matter, Schirmer, Lerma-Martin, 5XFAD and aging). Astrocytes gain Cd36 (10 of 14), Fabp7 (9), Abca1 and Lpl, lose Fa2h, Got1, Plin4, Mvk and Fdps, and lose ketone utilisation in every MS dataset and in Braak-6 AD.
- **Human MS by lesion type.** Lesion-associated microglia (all four snRNA-seq sets) are more glycolytic, more lipid-storing, less cholesterol-synthesising and less ketone-utilising than control white-matter microglia. MS oligodendrocytes have modestly *higher* cholesterol-synthesis scores than controls (Absinta lesion edges, Jäkel active lesions, Lerma-Martin), the opposite of mouse EAE, and more Hif1a (14–18 % vs 8 %), Txnip (37–59 % vs 32 %) and Ddit4 (9–10 % vs 2 %) in Lerma-Martin.
- **Spatial.** Chen 2020 AD sections have slightly higher lipid-synthesis and myelin-lipid scores than WT; across spots C4b correlates with Apoe (rank 9 of 21,586), Apod (39), Lpl (191), Fasn, Cers2, Fabp7, and the lipid pathways are the enriched ones. Lerma-Martin MS Visium lesion spots have much lower lipid / myelin-lipid scores than control white matter, i.e. loss of oligodendrocytes.
- **The C4b axis.** Metabolic correlations with C4b inside public mouse oligodendrocytes are weak (|mean rho| < 0.1 per pathway). Apod recurs (rank 8 in 5XFAD, 11 in aged WM, 18 after LPC / cuprizone, 42 in both aging sets), then Abca1, Acaca, Fdft1, Fa2h; lipid transport / storage is the enriched pathway in LPC / cuprizone, Serpina3n-cKO and 5XFAD; OXPHOS is mildly anti-correlated in aging and 5XFAD. **C4b-high oligodendrocytes have a higher lipid-transport / storage score in all seven mouse datasets** (cell-level p < 0.001; paired across samples p < 0.05 in LPC / cuprizone and 5XFAD), lower ketone utilisation in aged white matter and after LPC / cuprizone, lower glycolysis in 5XFAD.

## 9. Bottom line

1. Metabolic identity of the cell types is stable across platforms and species: neurons / astrocytes oxidative-glycolytic, oligodendrocytes lipid- and cholesterol-synthesising, microglia lipid-storing and glycolytic, OPCs the most ketone-utilising of the lineage.
2. Neuroinflammation (EAE) shifts the whole tissue, not one cell type: less cholesterol synthesis, less TCA, more lipid uptake and storage in nearly every cell type; near lesions the cholesterol-synthesis loss is steepest in OPCs and newly formed oligodendrocytes.
3. The disease-associated oligodendrocyte state, in situ and in sorted cells, is **less cholesterol / myelin-lipid synthesis, more lipid uptake and storage (Plin2 / Plin4, Abca1, Lpl, Cd36), a Hk2 / Pdk4 / Ldha glycolytic-switch signature, Bdh1 induction, less MCT1**. In human MS oligodendrocytes cholesterol synthesis is, if anything, up, and hypoxia / nutrient-stress genes (Hif1a, Txnip, Ddit4) are up.
4. The metabolic feature that travels with C4b everywhere is lipid handling (Apod, Plin4, Abca1), not energy metabolism and not ketone bodies. Ketone-body enzymes and the receptor Hcar2 belong to astrocytes / OPCs and microglia respectively; Hcar2⁺ and Slc16a3⁺ (MCT4) cells are enriched around C4b-high oligodendrocytes in EAE, as C5aR1⁺ cells were in OligoC4b.
5. Aging white matter shows a mild version of the same shift (less lipid / cholesterol synthesis, more Apoe / Plin2, less Ddit4 / Slc2a1).

## Caveats

1. Xenium AD has one section per genotype × age and only three metabolic genes.
2. The EAE 5K panel misses ketogenesis, ketolysis, most OXPHOS and most cholesterol-synthesis genes; those scores are absent or rest on 2–3 genes there.
3. In situ transcripts inside a segment can come from a neighbouring cell (Hcar2, Cd36, Slc16a3 in "oligodendrocytes"); the sorted Falcão cells are the control for that, and they are the replicates there (few animals), so Falcão p-values are descriptive.
4. Pathway scores are relative to size-matched control genes and depend on set composition (section 4); compare by difference, and check the gene tables.
5. The Visium age trend rests on 6 sections; pseudobulk tests need ≥3 samples per group and several public contrasts do not have that.
6. Coarse cell-type labels of the public datasets are marker-based for eight of them; only oligodendrocyte and microglia labels were validated (see OligoC4b `cell_type_annotation.md`).
7. Nuclei detect cytoplasmic / mitochondrial metabolic transcripts far less than whole cells; do not compare detection rates across modalities. Microglial ambient RNA in demyelinating-lesion nuclei (LPC, cuprizone) inflates Apoe, Lpl and Hcar2 in "oligodendrocytes" there.
8. Human C4B is not quantifiable (OligoC4b), so the C4b axis is mouse-only; everything else in the human data is direct.
