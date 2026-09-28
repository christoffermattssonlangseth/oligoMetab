# The metabolic gene panel

183 mouse genes in 21 pathway groups, defined in `scripts/oligometab.py` (`PANEL`). Human datasets use the upper-cased symbol unless `MOUSE_TO_HUMAN` says otherwise (e.g. `Gpi1` → `GPI`, `G6pdx` → `G6PD`, `Fh1` → `FH`, `Cyp51` → `CYP51A1`, `Atp5a1` → `ATP5F1A`, `Scd1` → `SCD`; `Scd2` has no one-to-one human ortholog and is skipped). `ALIASES` covers annotation-version differences (`Atp5a1`/`Atp5f1a`, `Hcar2`/`Gpr109a`).

Pathway **scores** (`SCORE_SETS`, `pathway_scores()`) are `scanpy.tl.score_genes`: the mean log-expression of the set minus that of a size-matched random control set drawn from expression-matched bins. They are on the log scale, can be negative, and are compared by difference. Single genes are compared by log2 fold change of pseudobulk means.

| Pathway group | n | Genes |
| --- | --- | --- |
| Ketone body synthesis | 4 | Hmgcs2, Hmgcl, Bdh1, Bdh2 |
| Ketone body utilisation | 5 | Oxct1, Acat1, Acat2, Acss1, Acss2 |
| Ketone / SCFA receptors | 3 | Hcar2, Ffar3, Ffar2 |
| Monocarboxylate transporters | 5 | Slc16a1, Slc16a7, Slc16a3, Slc16a6, Bsg |
| Glucose transporters | 4 | Slc2a1, Slc2a3, Slc2a4, Slc2a8 |
| Glycolysis | 20 | Hk1, Hk2, Hk3, Gpi1, Pfkl, Pfkm, Pfkp, Pfkfb2, Pfkfb3, Aldoa, Aldoc, Tpi1, Gapdh, Pgk1, Pgam1, Eno1, Eno2, Pkm, Ldha, Ldhb |
| Pyruvate fate | 11 | Pdha1, Pdhb, Pdk1, Pdk2, Pdk3, Pdk4, Pdp1, Mpc1, Mpc2, Pcx, Pck2 |
| Pentose phosphate pathway | 4 | G6pdx, Pgd, Tkt, Taldo1 |
| Glycogen | 4 | Gys1, Pygb, Gbe1, Agl |
| TCA cycle | 13 | Cs, Aco2, Idh2, Idh3a, Ogdh, Dlst, Sucla2, Suclg1, Sdha, Sdhb, Fh1, Mdh1, Mdh2 |
| Anaplerosis / amino acids | 8 | Got1, Got2, Glud1, Gls, Glul, Idh1, Slc25a1, Slc25a11 |
| Oxidative phosphorylation | 15 | Ndufa4, Ndufs1, Ndufv1, Ndufb8, Uqcrc1, Uqcrc2, Uqcrfs1, Cycs, Cox4i1, Cox5a, Cox6c, Cox7c, Atp5a1, Atp5b, Atp5o |
| Mitochondrial biogenesis / dynamics | 8 | Ppargc1a, Tfam, Nrf1, Opa1, Mfn2, Dnm1l, Pink1, Prkn |
| Fatty acid uptake / activation | 8 | Cd36, Slc27a1, Fabp5, Fabp7, Acsl1, Acsl3, Acsl6, Acsbg1 |
| Fatty acid beta-oxidation | 16 | Cpt1a, Cpt1c, Cpt2, Slc25a20, Acadm, Acadl, Acadvl, Acads, Hadha, Hadhb, Echs1, Acaa2, Decr1, Eci1, Acox1, Ehhadh |
| Lipid synthesis | 10 | Acly, Acaca, Fasn, Elovl1, Elovl5, Elovl6, Scd1, Scd2, Srebf1, Mlxipl |
| Cholesterol synthesis | 12 | Hmgcs1, Hmgcr, Mvk, Fdps, Fdft1, Sqle, Lss, Cyp51, Dhcr24, Dhcr7, Srebf2, Insig1 |
| Cholesterol / lipid transport | 12 | Ldlr, Abca1, Abca2, Abcg1, Apoe, Apod, Lpl, Nr1h2, Nr1h3, Plin2, Plin3, Plin4 |
| Myelin lipid synthesis | 6 | Ugt8a, Gal3st1, Cers2, Sptlc1, Pigt, Fa2h |
| Nutrient sensing / regulators | 15 | Hif1a, Epas1, Prkaa1, Prkaa2, Mtor, Rptor, Tsc2, Ddit4, Txnip, Sirt1, Sirt3, Ppara, Ppard, Pparg, Nfe2l2 |

## Focus genes (shown in every summary figure)

Hmgcs2, Bdh1, Oxct1, Acat1, Hcar2, Slc16a1, Slc16a7, Slc16a3, Slc2a1, Slc2a3, Hk1, Hk2, Pfkp, Aldoc, Pkm, Ldha, Ldhb, Pdk1, Cs, Sdha, Cox4i1, Atp5a1, Ppargc1a, Cpt1a, Acadm, Hadha, Fasn, Hmgcr, Srebf2, Plin2, Apoe, Hif1a, Txnip, Ddit4

## Score sets

| Score | Genes |
| --- | --- |
| Ketone synthesis | Hmgcs2, Hmgcl, Bdh1, Bdh2 |
| Ketone utilisation | Oxct1, Acat1, Acat2, Acss1, Acss2 |
| Glycolysis | Slc2a1, Slc2a3, Slc2a4, Slc2a8, Hk1, Hk2, Hk3, Gpi1, Pfkl, Pfkm, Pfkp, Pfkfb2, Pfkfb3, Aldoa, Aldoc, Tpi1, Gapdh, Pgk1, Pgam1, Eno1, Eno2, Pkm, Ldha, Ldhb |
| Pyruvate to lactate | Ldha, Ldhb, Pdk1, Pdk2, Pdk3, Pdk4, Slc16a1, Slc16a3, Slc16a7 |
| Pentose phosphate | G6pdx, Pgd, Tkt, Taldo1 |
| TCA cycle | Cs, Aco2, Idh2, Idh3a, Ogdh, Dlst, Sucla2, Suclg1, Sdha, Sdhb, Fh1, Mdh1, Mdh2 |
| OXPHOS | Ndufa4, Ndufs1, Ndufv1, Ndufb8, Uqcrc1, Uqcrc2, Uqcrfs1, Cycs, Cox4i1, Cox5a, Cox6c, Cox7c, Atp5a1, Atp5b, Atp5o |
| Beta-oxidation | Cpt1a, Cpt1c, Cpt2, Slc25a20, Acadm, Acadl, Acadvl, Acads, Hadha, Hadhb, Echs1, Acaa2, Decr1, Eci1, Acox1, Ehhadh |
| Lipid synthesis | Acly, Acaca, Fasn, Elovl1, Elovl5, Elovl6, Scd1, Scd2, Srebf1, Mlxipl |
| Cholesterol synthesis | Hmgcs1, Hmgcr, Mvk, Fdps, Fdft1, Sqle, Lss, Cyp51, Dhcr24, Dhcr7, Srebf2, Insig1 |
| Lipid transport / storage | Ldlr, Abca1, Abca2, Abcg1, Apoe, Apod, Lpl, Nr1h2, Nr1h3, Plin2, Plin3, Plin4 |
| Myelin lipids | Ugt8a, Gal3st1, Cers2, Sptlc1, Pigt, Fa2h |

## What is measurable where

| Platform | Panel genes present |
| --- | --- |
| Xenium mouse AD (347-gene panel) | 3 / 183: Apoe, Apod, Acsbg1 |
| Xenium mouse EAE (5K panel) | 100 / 183; missing ketogenesis (Hmgcs2, Hmgcl), Oxct1 / Acat1, most OXPHOS subunits, Apoe, most of the cholesterol pathway |
| Visium aging, Falcão, all public datasets | 177–183 / 183 (losses are lowly expressed genes filtered at build time, e.g. Prkn, Hcar2, Hmgcs2 in some human sets) |

Re-check with `python scripts/check_panel_symbols.py <files.h5ad>`.
