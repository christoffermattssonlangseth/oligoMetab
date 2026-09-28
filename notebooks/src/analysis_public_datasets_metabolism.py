# %% [markdown]
# # Metabolism of oligodendrocytes across sixteen public datasets: ketone bodies, glycolysis and lactate, mitochondrial oxidation, fatty acids and lipids
#
# Companion to `analysis_spatial_metabolism.ipynb` (in-house Xenium / Visium / Falcão data). This notebook asks the questions that
# `OligoC4b/analysis_public_datasets_complement.ipynb` asked about the complement system of the **metabolic machinery** instead, in the
# sixteen public datasets harmonised by OligoC4b: mouse AD (Park 2023, Zhou 2020 5XFAD), mouse aging (aging snRNA-seq HIP+CP, Ximerakis 2019,
# Kaya 2022 aged WM vs GM), toxic demyelination (LPC + cuprizone, Serpina3n-cKO cuprizone), mouse spatial AD (Chen 2020), human AD
# (Leng 2021, Sadick 2022), human MS single-nucleus (Jäkel 2019, Absinta 2021, Schirmer 2019, Lerma-Martin 2024) and human MS Visium
# (Lerma-Martin 2024, senescent-glia 2025). All are whole-transcriptome, so every panel gene is measurable, and, unlike C4A/C4B, the metabolic
# genes are quantifiable in the human data.
#
# **Gene panel** (`scripts/oligometab.py`, mouse symbols, 183 genes in 21 pathway groups): ketone-body synthesis (Hmgcs2, Hmgcl, Bdh1) and
# utilisation (Oxct1, Acat1), ketone / SCFA receptors (Hcar2, Ffar3), monocarboxylate transporters (Slc16a1 = MCT1, Slc16a7 = MCT2,
# Slc16a3 = MCT4), glucose transporters and glycolysis, pyruvate fate (Pdk1–4, Pdha1, Mpc1/2), pentose phosphate pathway, glycogen, TCA cycle,
# anaplerosis, oxidative phosphorylation, mitochondrial biogenesis, fatty-acid uptake and β-oxidation, lipid and cholesterol synthesis,
# cholesterol / lipid transport and storage (Apoe, Plin2), myelin lipid synthesis, and nutrient-sensing regulators (Hif1a, Prkaa1/2, Mtor, Txnip, Ddit4).
#
# Questions, per dataset:
#
# 1. Which coarse cell types express each pathway (pathway scores and fraction of cells with ≥1 UMI for the focus genes)?
# 2. Does expression change with disease / age / demyelination in every coarse cell type (oligodendrocytes, OPCs, microglia, astrocytes, neurons, endothelial, vascular, immune, ependymal), pseudobulk per sample, gene-wise and pathway-wise?
# 3. In the human datasets, how do the pathways behave by the authors' lesion / Braak labels?
# 4. In the spatial datasets, how do the pathways change with genotype, age or lesion type across spots?
# 5. Finally, as one axis among these: inside mouse oligodendrocytes, which metabolic genes and pathways co-vary with C4b, the marker of the disease-associated state characterised in OligoC4b, and do C4b-high cells differ from C4b-negative ones in pathway scores?
#
# Mouse symbols are used throughout; human genes are mapped through `oligometab.MOUSE_TO_HUMAN` / `ALIASES` (e.g. `Gpi1` → `GPI`, `Atp5a1` → `ATP5F1A` / `ATP5A1`).
# Datasets are processed one at a time and only summary tables are kept in memory. Pathway scores are `scanpy.tl.score_genes` (mean of the set minus a
# size-matched random control set, on the log scale) and are compared by difference; single genes are compared by log2 fold change of pseudobulk means.

# %%
import os, sys, gc, traceback, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "../../scripts")
try:
    from dotenv import load_dotenv
    load_dotenv("../../.env")
except Exception:
    pass
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, seaborn as sns, matplotlib.pyplot as plt
from scipy import stats
import oligometab as om
from oligometab import PANEL, ALL_GENES, FOCUS, SCORE_SETS, CONTEXT, PATHWAY_OF, short, heat

sc.settings.verbosity = 0
sc.set_figure_params(dpi=80, frameon=False)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 400)

PATHS = om.processed_paths()
only = os.getenv("OLIGOMETAB_ONLY_DATASETS")          # optional comma-separated subset, for smoke tests
if only:
    PATHS = {k: v for k, v in PATHS.items() if k in only.split(",")}
print(f"{len(PATHS)} datasets in {om.PROCESSED_DIR}")
for k, v in PATHS.items():
    print(f"  {k:42s} {v}")

# the focus genes in two blocks for readable heatmaps
FOCUS_ENERGY = ["Hmgcs2", "Bdh1", "Oxct1", "Acat1", "Hcar2", "Slc16a1", "Slc16a7", "Slc16a3", "Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkp", "Aldoc", "Pkm",
                "Ldha", "Ldhb", "Pdk1", "Cs", "Sdha", "Cox4i1", "Atp5a1", "Ppargc1a"]
FOCUS_LIPID = ["Cpt1a", "Acadm", "Hadha", "Fasn", "Hmgcr", "Srebf2", "Plin2", "Apoe", "Hif1a", "Txnip", "Ddit4"]
GENES = FOCUS + ["C4b", "Serpina3n"]
CTS_MAIN = ["Oligodendrocyte", "OPC", "Microglia", "Astrocyte", "Immune (lymphoid/myeloid)", "Neuron"]
CTS_TEST = [c for c in om.COARSE_TYPES if c != "Other"]          # every coarse cell type is tested; small ones drop out via the cell / sample filters
CTS_GENE_HEAT = ["Oligodendrocyte", "OPC", "Microglia", "Astrocyte", "Neuron"]

panel_table = pd.DataFrame([(k, len(v), ", ".join(v)) for k, v in PANEL.items()], columns=["pathway", "n_genes", "genes"]).set_index("pathway")
display(panel_table)

# %% [markdown]
# ## 1. Per-dataset pass
#
# For each dataset: load, record which panel genes resolve, compute pathway scores for every cell, the detection and pathway-score tables by cell type,
# the pseudobulk group comparisons (gene-wise and pathway-wise) for oligodendrocytes / OPCs / microglia / astrocytes, the C4b co-expression, genome-wide
# ranks and pathway-level rank enrichment inside oligodendrocytes, the C4b-high vs C4b-negative pathway-score comparison, the by-condition tables for
# human data and the spot-level analyses for spatial data. Per-dataset figures are shown here; cross-dataset summaries follow in sections 2–7 (cell types, disease contrasts in every cell type, human conditions, spatial, then the C4b axis).

# %%
DET, DE, DE_SCORE, SCORE_CT = [], [], [], []
CO, RANKS, ENRICH, TOP, C4B_SCORE, COND, COND_SCORE, SPATIAL, INFO, AVAIL, ERRORS = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
np.random.seed(0)


def pv(df, index, columns, values):
    """NaN-preserving pivot (pivot_table drops rows whose values are all NaN, e.g. p-values of underpowered contrasts)."""
    idx = [index] if isinstance(index, str) else list(index)
    return df.drop_duplicates(idx + [columns]).set_index(idx + [columns])[values].unstack(columns)


def star(d, p):
    """Annotation matrix: rounded value plus '*' where p < 0.05 (p aligned to d)."""
    p = p.reindex(index=d.index, columns=d.columns)
    return (d.round(2).astype(str).replace("nan", "") + np.where(p.values < 0.05, "*", "")).values


def c4b_high_vs_neg_scores(ol, P, sample_col="sample", high_q=0.75, min_cells=20):
    """Per pathway: mean score in C4b-high (top quartile of C4b+) vs C4b-negative oligodendrocytes, per-sample paired Wilcoxon."""
    c4b = om.expr(ol, "C4b")
    pos = c4b > 0
    if pos.sum() < 50:
        return pd.DataFrame()
    thr = np.quantile(c4b[pos], high_q)
    high, neg = c4b >= thr, c4b == 0
    rows = []
    smp = ol.obs[sample_col].astype(str).values
    for pw in P.columns:
        v = P[pw].values
        per_sample = []
        for s in np.unique(smp):
            m = smp == s
            if (m & high).sum() >= min_cells and (m & neg).sum() >= min_cells:
                per_sample.append((v[m & high].mean(), v[m & neg].mean()))
        ps = np.array(per_sample)
        p = stats.wilcoxon(ps[:, 0], ps[:, 1]).pvalue if len(ps) >= 3 and not np.allclose(ps[:, 0], ps[:, 1]) else np.nan
        d_cell = v[high].mean() - v[neg].mean()
        d_cell_p = stats.mannwhitneyu(v[high], v[neg]).pvalue
        rows.append({"pathway": pw, "mean_C4b_high": v[high].mean(), "mean_C4b_neg": v[neg].mean(), "diff_cells": d_cell, "MWU_p_cells": d_cell_p,
                     "n_samples_paired": len(ps), "diff_paired_median": np.median(ps[:, 0] - ps[:, 1]) if len(ps) else np.nan, "wilcoxon_p_samples": p})
    return pd.DataFrame(rows).set_index("pathway")


for name, path in PATHS.items():
    try:
        a = sc.read_h5ad(path)
        a.obs["cell_type_coarse"] = a.obs["cell_type_coarse"].astype(str)
        species, modality = str(a.obs["species"].iloc[0]), str(a.obs["modality"].iloc[0])
        ref, alts = om.alt_groups(a)
        INFO[name] = {"n_cells": a.n_obs, "n_genes": a.n_vars, "species": species, "modality": modality, "n_samples": a.obs["sample"].nunique(),
                      "groups": a.obs["group"].value_counts().to_dict()}
        print(f"\n{'=' * 120}\n{name}: {a.n_obs} cells x {a.n_vars} genes | {species} {modality} | samples={INFO[name]['n_samples']} | groups={INFO[name]['groups']}")
        r_all = om.resolve(a, ALL_GENES)
        AVAIL[name] = {g: r_all[g] is not None for g in ALL_GENES}
        miss = [g for g in ALL_GENES if r_all[g] is None]
        print(f"panel genes resolved: {len(ALL_GENES) - len(miss)} / {len(ALL_GENES)}; missing: {miss}")

        # -- pathway scores for every cell (log scale, set mean minus matched control)
        P = om.pathway_scores(a, SCORE_SETS)
        print("pathway scores computed:", list(P.columns))

        if modality != "spatial":
            # -- 1) detection and pathway scores by cell type
            t = om.detection_table(a, GENES, "cell_type_coarse", min_cells=30)
            t.insert(0, "dataset", name)
            DET.append(t.reset_index().rename(columns={"group": "cell_type"}))
            display(t.loc[t.index.isin(CTS_MAIN), ["n_cells"] + [g for g in FOCUS_ENERGY if g in t.columns]].round(3))
            s_ct = P.groupby(a.obs["cell_type_coarse"].values).mean()
            n_ct = a.obs["cell_type_coarse"].value_counts()
            s_ct = s_ct.loc[n_ct[n_ct >= 30].index.intersection(s_ct.index)]
            s_ct.insert(0, "dataset", name)
            SCORE_CT.append(s_ct.reset_index().rename(columns={"index": "cell_type"}))
            keep_ct = [c for c in CTS_TEST if c in s_ct.index]
            heat(s_ct.loc[keep_ct].drop(columns="dataset"), f"{short(name)}: mean pathway score by cell type", "pathway score (log scale)",
                 fmt=".2f", cmap="RdBu_r", center=0, row_labels=keep_ct)

            # -- 2) pseudobulk group comparison: genes and pathway scores
            for ct in CTS_TEST:
                m_ct = (a.obs["cell_type_coarse"] == ct).values
                if m_ct.sum() < 50 or not alts:
                    continue
                sub = a[m_ct]
                n = sub.obs.groupby("sample", observed=True).size()
                ok_samples = n[n >= 20].index
                pb = om.pseudobulk(sub, ALL_GENES)
                pb = pb[pb["sample"].isin(ok_samples)]
                pbs = om.pseudobulk(sub, [], values=P.loc[sub.obs_names])
                pbs = pbs[pbs["sample"].isin(ok_samples)]
                for alt in alts:
                    for g in [c for c in pb.columns if c in ALL_GENES]:
                        ma, mr, p, na, nr = om.group_test(pb, g, "group", alt, ref)
                        DE.append({"dataset": name, "cell_type": ct, "gene": g, "pathway": PATHWAY_OF[g], "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr,
                                   "log2FC": om.log2fc(ma, mr), "MWU_p": p, "n_alt": na, "n_ref": nr})
                    for pw in P.columns:
                        ma, mr, p, na, nr = om.group_test(pbs, pw, "group", alt, ref)
                        DE_SCORE.append({"dataset": name, "cell_type": ct, "pathway": pw, "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr,
                                         "diff": ma - mr, "MWU_p": p, "n_alt": na, "n_ref": nr})
            d = pd.DataFrame([x for x in DE_SCORE if x["dataset"] == name])
            if len(d):
                piv = pv(d, ["cell_type", "alt"], "pathway", "diff")
                pp = pv(d, ["cell_type", "alt"], "pathway", "MWU_p")
                heat(piv, f"{short(name)}: pathway-score difference ({' / '.join(alts)} vs {ref}), pseudobulk per sample", "score difference (* MWU p < 0.05)",
                     annot=star(piv, pp), vmin=-0.3, vmax=0.3, row_labels=[f"{c} | {al} vs {ref}" for c, al in piv.index])

            # -- 3) C4b program inside oligodendrocytes (mouse; human C4B is not quantifiable)
            m_ol = (a.obs["cell_type_coarse"] == "Oligodendrocyte").values
            ol = a[m_ol].copy()
            r = om.resolve(ol, ["C4b"])
            n_c4b = int((om.expr(ol, "C4b") > 0).sum()) if r["C4b"] is not None else 0
            if ol.n_obs >= 100 and r["C4b"] is not None and n_c4b >= 50 and species == "mouse":
                co = om.coexpression_with_anchor(ol, "C4b", [g for g in ALL_GENES + ["Serpina3n"]])
                rho = om.genome_wide_rho(ol, r["C4b"], max_cells=40000)
                rk = om.rank_panel(rho, ol, ALL_GENES + ["Serpina3n"])
                en = om.pathway_rank_enrichment(rho, ol, PANEL)
                CO[name], RANKS[name], ENRICH[name], TOP[name] = co, rk, en, list(rho.head(25).index)
                print(f"oligodendrocytes: {ol.n_obs} ({n_c4b / ol.n_obs:.1%} C4b+); {len(rho)} genes ranked by Spearman correlation with C4b")
                print("top-25 C4b-correlated genes:", ", ".join(TOP[name]))
                print("pathway-level rank enrichment (mean rho of pathway genes; MWU vs all other genes):"); display(en.round(4))
                best = rk.sort_values("rank").head(15)
                print("15 best-ranked metabolic genes:"); display(best.round(4))
                cs = c4b_high_vs_neg_scores(ol, P.loc[ol.obs_names])
                C4B_SCORE[name] = cs
                print("pathway scores in C4b-high vs C4b-negative oligodendrocytes:"); display(cs.round(4))
            else:
                print(f"oligodendrocytes: {ol.n_obs}; C4b analyses skipped ({'human, C4B not quantifiable' if species == 'human' else f'{n_c4b} C4b+ cells'})")
            del ol

            # -- 4) human / multi-level conditions: by the authors' labels
            key = None
            if "condition_original" in a.obs and a.obs["condition_original"].nunique() > 1:
                key = "condition_original"
            elif "braak" in a.obs:
                key = "braak"
            elif len(alts) > 1:
                key = "group"
            if key is not None:
                for ct in ["Oligodendrocyte", "Microglia", "Astrocyte"]:
                    m_ct = (a.obs["cell_type_coarse"] == ct).values
                    if m_ct.sum() < 50:
                        continue
                    sub = a[m_ct]
                    t = om.detection_table(sub, GENES, key, min_cells=30)
                    COND[(name, ct)] = t
                    s = P.loc[sub.obs_names].groupby(sub.obs[key].astype(str).values).mean()
                    nn = sub.obs[key].astype(str).value_counts()
                    s = s.loc[nn[nn >= 30].index.intersection(s.index)]
                    COND_SCORE[(name, ct)] = s
                    if ct == "Oligodendrocyte":
                        print(f"-- {ct}: mean pathway score by {key}"); display(s.round(3))

            # -- dotplot of the focus panel by cell type
            r = om.resolve(a, FOCUS + ["Plp1", "Hexb", "Aqp4"])
            hv = [r[g] for g in FOCUS + ["Plp1", "Hexb", "Aqp4"] if r[g] is not None]
            keep = a.obs["cell_type_coarse"].isin(CTS_MAIN)
            sc.pl.dotplot(a[keep], var_names=hv, groupby="cell_type_coarse", standard_scale="var", color_map="Reds",
                          figsize=(0.32 * len(hv) + 2, 3), show=False, title=short(name))
            plt.show()

        else:
            # -- spatial: means by condition, pathway scores by condition, spot-level C4b co-expression, maps
            r = om.resolve(a, ["C4b"])
            key = "group_age" if "age" in a.obs else ("condition_original" if "condition_original" in a.obs else "group")
            if key == "group_age":
                a.obs["group_age"] = a.obs["group"].astype(str) + "_" + a.obs["age"].astype(str)
            print("groups:", a.obs["group"].value_counts().to_dict(), "| key:", key, a.obs[key].value_counts().to_dict())
            mt = om.mean_table(a, GENES, key, min_cells=30)
            SPATIAL[(name, "means")] = mt
            s = P.groupby(a.obs[key].astype(str).values).mean()
            SPATIAL[(name, "scores")] = s
            print(f"mean pathway score by {key}:"); display(s.round(3))
            # pseudobulk per sample for the two-level group
            pbs = om.pseudobulk(a, [], values=P)
            for alt in alts:
                for pw in P.columns:
                    ma, mr, p, na, nr = om.group_test(pbs, pw, "group", alt, ref)
                    DE_SCORE.append({"dataset": name, "cell_type": "spot", "pathway": pw, "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr,
                                     "diff": ma - mr, "MWU_p": p, "n_alt": na, "n_ref": nr})
            pb = om.pseudobulk(a, ALL_GENES)
            for alt in alts:
                for g in [c for c in pb.columns if c in ALL_GENES]:
                    ma, mr, p, na, nr = om.group_test(pb, g, "group", alt, ref)
                    DE.append({"dataset": name, "cell_type": "spot", "gene": g, "pathway": PATHWAY_OF[g], "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr,
                               "log2FC": om.log2fc(ma, mr), "MWU_p": p, "n_alt": na, "n_ref": nr})
            if r["C4b"] is not None and (om.expr(a, "C4b") > 0).sum() >= 50 and species == "mouse":
                co = om.coexpression_with_anchor(a, "C4b", ALL_GENES + ["Serpina3n"])
                rho = om.genome_wide_rho(a, r["C4b"], max_cells=40000)
                rk = om.rank_panel(rho, a, ALL_GENES + ["Serpina3n"])
                en = om.pathway_rank_enrichment(rho, a, PANEL)
                SPATIAL[(name, "coexpr")] = co.join(rk[["rank"]])
                SPATIAL[(name, "enrich")] = en
                TOP[name] = list(rho.head(30).index)
                print("pathway-level enrichment among C4b-correlated genes across spots:"); display(en.round(4))
                print("15 best-ranked metabolic genes across spots:"); display(rk.sort_values("rank").head(15).round(4))
                print("top-30 C4b-correlated genes across spots:", ", ".join(TOP[name]))
            else:
                print("C4b not quantifiable here; spot-level C4b co-expression skipped")
            if "spatial" in a.obsm:
                pick = a.obs.groupby("group", observed=True)["sample"].agg(lambda s: s.value_counts().index[0])
                rr = om.resolve(a, ["C4b", "Slc16a1", "Ldha", "Hk2", "Cox4i1", "Apoe", "Plin2"])
                hv = [rr[g] for g in ["C4b", "Slc16a1", "Ldha", "Hk2", "Cox4i1", "Apoe", "Plin2"] if rr[g] is not None]
                for grp, smp in pick.items():
                    sub = a[a.obs["sample"] == smp]
                    span = np.ptp(sub.obsm["spatial"], axis=0).max()
                    sc.pl.spatial(sub, color=hv, spot_size=span / 80, cmap="magma", vmax="p99", ncols=len(hv), show=False)
                    plt.suptitle(f"{short(name)}: {smp} ({grp})", y=1.02); plt.show()
                    for pw in ["Glycolysis", "OXPHOS", "Ketone utilisation"]:
                        if pw in P.columns:
                            sub.obs[pw] = P.loc[sub.obs_names, pw].values
                    sc.pl.spatial(sub, color=[pw for pw in ["Glycolysis", "OXPHOS", "Ketone utilisation"] if pw in sub.obs], spot_size=span / 80,
                                  cmap="RdBu_r", vcenter=0, ncols=3, show=False)
                    plt.suptitle(f"{short(name)}: {smp} ({grp}) pathway scores", y=1.02); plt.show()
    except Exception as exc:
        ERRORS[name] = traceback.format_exc()
        print(f"!! {name} failed: {type(exc).__name__}: {exc}")
    finally:
        for v in ["a", "P", "sub", "ol"]:
            if v in globals():
                del globals()[v]
        gc.collect()

DET = pd.concat(DET, ignore_index=True) if DET else pd.DataFrame()
DE, DE_SCORE = pd.DataFrame(DE), pd.DataFrame(DE_SCORE)
SCORE_CT = pd.concat(SCORE_CT, ignore_index=True) if SCORE_CT else pd.DataFrame()
print("\ndone:", len(INFO), "datasets;", len(ERRORS), "failed", list(ERRORS))
for k, v in ERRORS.items():
    print(f"\n--- {k}\n{v}")

# %% [markdown]
# ## 2. Panel availability and which cell types express the pathways

# %%
avail = pd.DataFrame(AVAIL).reindex(ALL_GENES)
avail.index.name = "gene"
n_missing = (~avail).sum()
print("genes not resolved per dataset:"); display(n_missing.rename("n_missing").to_frame().T)
miss_long = {ds: [g for g in ALL_GENES if not avail.loc[g, ds]] for ds in avail.columns}
display(pd.Series({short(k): ", ".join(v) for k, v in miss_long.items() if v}, name="missing genes").to_frame())

# %%
# pathway score by cell type, averaged over the mouse single-cell / single-nucleus datasets (each dataset z-scored across cell types first)
if len(SCORE_CT):
    mouse_ds = [k for k, v in INFO.items() if v["species"] == "mouse" and v["modality"] != "spatial"]
    blocks = []
    for ds in mouse_ds:
        s = SCORE_CT[SCORE_CT.dataset == ds].set_index("cell_type").drop(columns="dataset")
        s = s.loc[[c for c in CTS_TEST if c in s.index]]
        blocks.append((s - s.mean()) / s.std(ddof=0))
    z = pd.concat(blocks).groupby(level=0).mean().loc[[c for c in CTS_TEST if c in pd.concat(blocks).index]]
    heat(z.T, "Pathway score by cell type, z-scored within dataset and averaged over the mouse datasets", "z-score across cell types", fmt=".1f",
         cmap="RdBu_r", center=0, row_labels=list(z.columns), col_labels=list(z.index), annot_size=8)
    # raw scores in oligodendrocytes and microglia per dataset
    for ct in [c for c in CTS_TEST if c in SCORE_CT.cell_type.unique()]:
        s = SCORE_CT[SCORE_CT.cell_type == ct].set_index("dataset").drop(columns="cell_type")
        heat(s, f"{ct}: mean pathway score per dataset (log scale; 0 = same as size-matched control genes)", "pathway score", fmt=".2f", cmap="RdBu_r", center=0)

# %%
# detection of the focus genes in oligodendrocytes and microglia
if len(DET):
    for ct in ["Oligodendrocyte", "Microglia"]:
        h = DET[DET.cell_type == ct].set_index("dataset")
        for label, block in [("energy", FOCUS_ENERGY), ("lipid / regulators", FOCUS_LIPID)]:
            cols = [g for g in block if g in h.columns]
            heat(h[cols], f"{ct}: fraction of cells with ≥1 UMI ({label} genes)", "fraction of cells detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=1)
    display(DET[DET.cell_type.isin(["Oligodendrocyte", "Microglia"])].set_index(["dataset", "cell_type"])[["n_cells"] + [g for g in FOCUS if g in DET.columns]].round(3))

# %% [markdown]
# ## 3. Disease / age / demyelination effects (pseudobulk per sample)
#
# Groups are compared with a Mann–Whitney test when both have ≥3 samples (`*`); otherwise the value is descriptive. Datasets with more than one contrast
# appear once per contrast. Pathway scores are compared by difference; single genes by log2 fold change of the pseudobulk means.

# %%
if len(DE_SCORE):
    DE_SCORE["contrast"] = DE_SCORE["dataset"].map(short) + "  [" + DE_SCORE["alt"] + " vs " + DE_SCORE["ref"] + "]"
    for ct in [c for c in CTS_TEST + ["spot"] if c in DE_SCORE.cell_type.unique()]:
        sub_ = DE_SCORE[DE_SCORE.cell_type == ct]
        if not len(sub_):
            continue
        d = pv(sub_, "contrast", "pathway", "diff")
        d = d[[c for c in SCORE_SETS if c in d.columns]]
        p = pv(sub_, "contrast", "pathway", "MWU_p")
        heat(d, f"{ct}: pathway-score difference (disease / aged / demyelinated vs reference)", "score difference (* Mann–Whitney p < 0.05)",
             annot=star(d, p), vmin=-0.3, vmax=0.3, row_labels=list(d.index))

# %%
if len(DE):
    DE["contrast"] = DE["dataset"].map(short) + "  [" + DE["alt"] + " vs " + DE["ref"] + "]"
    for ct in [c for c in CTS_GENE_HEAT if c in DE.cell_type.unique()]:
        for label, block in [("energy", FOCUS_ENERGY), ("lipid / regulators", FOCUS_LIPID)]:
            sub_ = DE[DE.cell_type == ct]
            if not len(sub_):
                continue
            d = pv(sub_, "contrast", "gene", "log2FC")
            d = d[[g for g in block if g in d.columns]]
            p = pv(sub_, "contrast", "gene", "MWU_p")
            heat(d, f"{ct}: pseudobulk log2 fold change, {label} genes", "log2 fold change  (* Mann–Whitney p < 0.05)", annot=star(d, p), vmin=-2, vmax=2,
                 row_labels=list(d.index))

# %%
# consistency across contrasts, per cell type: how many contrasts move each gene up / down (|log2FC| > 0.5), and how many of those are significant
if len(DE):
    CONS = {}
    for ct in [c for c in CTS_TEST if c in DE.cell_type.unique()]:
        d_ = DE[DE.cell_type == ct].copy()
        d_["sig_up"] = (d_.MWU_p < 0.05) & (d_.log2FC > 0); d_["sig_down"] = (d_.MWU_p < 0.05) & (d_.log2FC < 0)
        cons = d_.groupby(["pathway", "gene"]).agg(n_contrasts=("log2FC", "size"), n_up=("log2FC", lambda x: int((x > 0.5).sum())),
                                                   n_down=("log2FC", lambda x: int((x < -0.5).sum())), median_log2FC=("log2FC", "median"),
                                                   n_sig_up=("sig_up", "sum"), n_sig_down=("sig_down", "sum"))
        cons["net"] = cons["n_up"] - cons["n_down"]; CONS[ct] = cons
        print(f"\n=== {ct}: most consistently UP across {d_.contrast.nunique()} contrasts"); display(cons.sort_values(["net", "median_log2FC"], ascending=False).head(12).round(3))
        print(f"=== {ct}: most consistently DOWN"); display(cons.sort_values(["net", "median_log2FC"]).head(12).round(3))
    sig = DE[(DE.MWU_p < 0.05) & (DE.cell_type != "spot")].sort_values(["cell_type", "pathway", "gene", "dataset"])
    print(f"\n{len(sig)} significant gene-level contrasts across all cell types (full table in results/public_gene_pseudobulk_contrasts.csv); oligodendrocytes:")
    display(sig[sig.cell_type == "Oligodendrocyte"][["dataset", "gene", "pathway", "alt", "ref", "mean_alt", "mean_ref", "log2FC", "MWU_p", "n_alt", "n_ref"]].round(4).reset_index(drop=True))

# %% [markdown]
# ## 4. Human datasets by the authors' condition labels
#
# Pathway scores and focus-gene detection in oligodendrocytes and microglia by lesion type (MS) or Braak stage / diagnosis (AD). The metabolic genes are
# quantifiable in human nuclei (unlike C4A/C4B), so these are direct read-outs.

# %%
for (name, ct), s in COND_SCORE.items():
    if INFO[name]["species"] != "human" or ct not in ("Oligodendrocyte", "Microglia"):
        continue
    heat(s[[c for c in SCORE_SETS if c in s.columns]], f"{short(name)} — {ct}: mean pathway score by condition", "pathway score", fmt=".2f", cmap="RdBu_r", center=0, row_labels=list(s.index))
for (name, ct), t in COND.items():
    if INFO[name]["species"] != "human" or ct != "Oligodendrocyte":
        continue
    print(f"=== {short(name)} — {ct}: fraction detected by condition"); display(t[["n_cells"] + [g for g in FOCUS if g in t.columns]].round(3))

# %% [markdown]
# ## 5. Spatial datasets: metabolism across spots

# %%
for (name, kind), t in SPATIAL.items():
    print(f"=== {short(name)} — {kind}")
    display(t.round(4) if kind != "coexpr" else t.sort_values("rank")[["pathway", "spearman_rho", "rank", "frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "fisher_p"]].round(4).head(40))

# %% [markdown]
# ## 6. One axis among many: the C4b⁺ oligodendrocyte program and metabolism (mouse datasets)
#
# Spearman correlation of every gene with C4b inside oligodendrocytes, then (a) the correlation of the focus genes, (b) pathway-level rank enrichment
# (mean rho of the pathway's genes; Mann–Whitney against all other ranked genes), (c) the best-ranked metabolic genes per dataset, and (d) pathway
# scores in C4b-high vs C4b-negative oligodendrocytes with a per-sample paired test.

# %%
if CO:
    R = pd.DataFrame({k: v["spearman_rho"] for k, v in CO.items()}).T
    for label, block in [("energy", FOCUS_ENERGY), ("lipid / regulators", FOCUS_LIPID + ["Serpina3n"])]:
        cols = [g for g in block if g in R.columns]
        heat(R[cols], f"Spearman correlation with C4b inside oligodendrocytes ({label} genes; Serpina3n = positive control)", "Spearman rho with C4b", fmt=".2f", vmin=-0.2, vmax=0.2)
    E = pd.DataFrame({k: v["mean_rho"] for k, v in ENRICH.items()}).T
    Ep = pd.DataFrame({k: v["MWU_p"] for k, v in ENRICH.items()}).T.reindex(columns=E.columns)
    E = E[[c for c in PANEL if c in E.columns]]; Ep = Ep[E.columns]
    annot = E.round(2).astype(str).replace("nan", "") + np.where(Ep < 0.05, "*", "")
    heat(E, "Pathway-level enrichment among C4b-correlated genes in oligodendrocytes (mean rho of pathway genes; * MWU p < 0.05 vs all other genes)",
         "mean Spearman rho with C4b", annot=annot.values, vmin=-0.1, vmax=0.1)
    RK = pd.DataFrame({k: v["rank"] for k, v in RANKS.items()}).T
    n_ranked = pd.Series({k: v["n_ranked"].iloc[0] for k, v in RANKS.items()})
    print("metabolic genes ranked in the top 300 C4b-correlated genes of any dataset (rank; 1 = most C4b-correlated):")
    top = RK.loc[:, (RK <= 300).any(axis=0)].T
    top["pathway"] = [PATHWAY_OF.get(g, "context") for g in top.index]
    display(top.sort_values(list(RK.index)[0]).astype({c: "Int64" for c in RK.index}))
    print("genes ranked (per dataset):"); display(n_ranked.to_frame("n_ranked").T)
    print("\ntop-25 C4b-correlated genes per dataset (any gene, for context):")
    for k, v in TOP.items():
        if INFO[k]["modality"] != "spatial":
            print(f"  {short(k)}: {', '.join(v)}")

# %%
if C4B_SCORE:
    D = pd.DataFrame({k: v["diff_cells"] for k, v in C4B_SCORE.items() if len(v)}).T
    Dp = pd.DataFrame({k: v["wilcoxon_p_samples"] for k, v in C4B_SCORE.items() if len(v)}).T.reindex(columns=D.columns)
    Dn = pd.DataFrame({k: v["n_samples_paired"] for k, v in C4B_SCORE.items() if len(v)}).T.reindex(columns=D.columns)
    D = D[[c for c in SCORE_SETS if c in D.columns]]; Dp = Dp[D.columns]
    annot = D.round(2).astype(str).replace("nan", "") + np.where(Dp < 0.05, "*", "")
    heat(D, "Pathway score in C4b-high minus C4b-negative oligodendrocytes (cell-level difference; * paired Wilcoxon across samples p < 0.05)",
         "score difference", annot=annot.values, vmin=-0.3, vmax=0.3)
    display(pd.concat(C4B_SCORE, names=["dataset", "pathway"]).round(4))

# %% [markdown]
# ## 7. Cross-dataset summary

# %%
display(pd.DataFrame(INFO).T)
if len(DE_SCORE):
    summ = pv(DE_SCORE[DE_SCORE.cell_type == "Oligodendrocyte"], "contrast", "pathway", "diff")
    summ = summ[[c for c in SCORE_SETS if c in summ.columns]]
    print("Oligodendrocyte pathway-score differences per contrast (also shown as a heatmap in section 3):"); display(summ.round(3))
    sig = DE_SCORE[(DE_SCORE.MWU_p < 0.05) & (DE_SCORE.cell_type != "spot")].sort_values(["cell_type", "pathway", "dataset"])
    print(f"{len(sig)} significant pathway-level contrasts:"); display(sig[["dataset", "cell_type", "pathway", "alt", "ref", "mean_alt", "mean_ref", "diff", "MWU_p", "n_alt", "n_ref"]].round(4).reset_index(drop=True))

# %%
# tables for the report / docs
os.makedirs("../../results", exist_ok=True)
if len(DE): DE.to_csv("../../results/public_gene_pseudobulk_contrasts.csv", index=False)
if len(DE_SCORE): DE_SCORE.to_csv("../../results/public_pathway_pseudobulk_contrasts.csv", index=False)
if len(DET): DET.to_csv("../../results/public_detection_by_celltype.csv", index=False)
if len(SCORE_CT): SCORE_CT.to_csv("../../results/public_pathway_scores_by_celltype.csv", index=False)
if CO: pd.concat(CO, names=["dataset", "gene"]).to_csv("../../results/public_c4b_coexpression_oligodendrocytes.csv")
if RANKS: pd.concat(RANKS, names=["dataset", "gene"]).to_csv("../../results/public_c4b_rank_oligodendrocytes.csv")
if ENRICH: pd.concat(ENRICH, names=["dataset", "pathway"]).to_csv("../../results/public_c4b_pathway_enrichment_oligodendrocytes.csv")
if C4B_SCORE: pd.concat({k: v for k, v in C4B_SCORE.items() if len(v)}, names=["dataset", "pathway"]).to_csv("../../results/public_c4b_high_vs_neg_pathway_scores.csv")
if COND_SCORE: pd.concat({f"{k[0]}|{k[1]}": v for k, v in COND_SCORE.items()}, names=["dataset|cell_type", "condition"]).to_csv("../../results/public_pathway_scores_by_condition.csv")
print("results written to ../../results/")

# %% [markdown]
# ### Interpretation
#
# **Scope.** Sixteen public datasets (about one million cells, nuclei and spots), the harmonised objects from OligoC4b, 177–183 of the 183 panel genes resolved in each. Every coarse cell type was tested; pathway scores are compared by difference and genes by pseudobulk log2 fold change, with Mann–Whitney across samples where both groups have ≥3 samples.
#
# **Division of labour.** Across the seven mouse single-cell / single-nucleus datasets the same profile emerges as in the in-house Xenium EAE data: neurons and astrocytes carry the highest glycolysis, TCA and OXPHOS scores; oligodendrocytes the highest myelin-lipid, cholesterol- and lipid-synthesis scores (with OPCs next); microglia and immune cells the highest lipid-transport / storage and pentose-phosphate scores and the lowest cholesterol synthesis; ketone-body utilisation is highest in astrocytes and OPCs and low in microglia. Hcar2 is detected almost only in microglia (4–9 % in the droplet datasets) and Slc16a1 (MCT1) in 29–57 % of oligodendrocytes in the whole-cell datasets but far less in nuclei, a modality effect that also lowers every detection rate in the human snRNA-seq sets.
#
# **Disease, age and demyelination effects are small at the pathway level and consistent at the gene level.** Pathway-score differences in oligodendrocytes are mostly within ±0.05 (compared with ±0.15–0.25 for the in-situ EAE contrasts), with a few clear ones: MS oligodendrocytes in Schirmer 2019 are more glycolytic (+0.16) and oxidative (+0.11) with more lipid storage (+0.08); aged oligodendrocytes in Ximerakis 2019 have lower OXPHOS (−0.06), cholesterol synthesis (−0.07) and myelin-lipid synthesis (−0.14) and slightly more lipid storage (+0.03); ketone-body synthesis and utilisation scores are slightly but significantly lower in MS oligodendrocytes in three of four human MS datasets (Absinta, Jäkel, Lerma-Martin); Braak-6 oligodendrocytes in Leng 2021 have less lipid synthesis (−0.09); after cuprizone and LPC (2 samples each, descriptive) lipid storage rises (+0.13 / +0.27) and lipid / myelin-lipid synthesis fall. At the gene level the same handful of genes moves in the same direction in most of the 14 oligodendrocyte contrasts: **Lpl up in 10, Hk2 up in 10, Apoe, Plin2 and Fabp5 up in 6–7; Slc16a1 (MCT1) down in 5, Gal3st1 down in 6, Acsl6, Ppara, Lss, Srebf1, Cpt1c, Scd1 down in 4–5.** Bdh1 falls in 5 of 14 public oligodendrocyte contrasts, in contrast to its induction in EAE-associated oligodendrocytes (in-house notebook), so ketone-body handling does not move in one direction across conditions. Microglia show the most reproducible signature of all: **Lpl up in 12 of 13 contrasts (median log2FC 2.8), Apoe and Plin2 up in 8, Pparg and Cd36 in 6–8, Hcar2 in 6 of 10**, i.e. the lipid-laden, PPARγ / LXR-driven disease-associated microglial program, visible as a lipid-transport / storage score increase in LPC, cuprizone, Jäkel MS, aged white matter, Schirmer MS, Lerma-Martin MS, 5XFAD and aging. Astrocytes gain Cd36 (10 of 14), Fabp7 (9), Abca1 and Lpl, and lose Fa2h, Got1, Plin4 and the cholesterol-synthesis genes Mvk and Fdps; their ketone-utilisation score falls in every MS dataset and in Braak-6 AD. Neurons, endothelial and vascular cells move little; endothelial cells lose glycogen-handling genes (Gys1, Pygb) in most contrasts.
#
# **Human MS by lesion type.** Microglia at chronic active lesion edges and in periplaque white matter (Absinta) and in all lesion types (Jäkel, Schirmer, Lerma-Martin) have higher glycolysis and lipid-storage and lower cholesterol-synthesis and ketone-utilisation scores than control white-matter microglia. Oligodendrocytes in MS lesions show modestly higher cholesterol-synthesis scores than controls (Absinta lesion edges, Jäkel active lesions, Lerma-Martin), the opposite of the mouse EAE oligodendrocytes, and higher detection of the hypoxia / nutrient-stress genes Hif1a (14–18 % vs 8 % in Lerma-Martin), Txnip (37–59 % vs 32 %) and Ddit4 (9–10 % vs 2 %). Jäkel control nuclei are deeper than the lesion nuclei, so its detection tables are confounded by depth; the Schirmer chronic-inactive oligodendrocyte glycolysis / OXPHOS increase rests on ~1,000 nuclei from few donors.
#
# **Spatial datasets.** In the Chen 2020 AD Spatial Transcriptomics data the AD sections have slightly higher lipid-synthesis and myelin-lipid scores than WT (p < 0.05 over 10 vs 10 sections) and across spots C4b correlates with Apoe (rank 9 of 21,586 genes), Apod (39), Lpl (191), Fasn, Cers2 and Fabp7, with lipid transport / storage, myelin-lipid synthesis and lipid synthesis the top enriched pathways; glycolysis, cholesterol synthesis and OXPHOS are also mildly enriched, TCA and the ketone pathways are not. In the Lerma-Martin MS Visium data lesion spots have much lower lipid-synthesis (−0.31), myelin-lipid (−0.22) and lipid-storage scores than control white matter, which is loss of oligodendrocytes in demyelinated tissue rather than a per-cell change.
#
# **The C4b axis in public mouse oligodendrocytes.** Correlations with C4b inside oligodendrocytes are weak for metabolic genes everywhere (|mean rho| < 0.1 per pathway). The one metabolic gene that recurs among the C4b correlates is Apod (rank 8 in 5XFAD, 11 in aged white matter, 18 after LPC / cuprizone, 42 in both aging datasets), followed by Abca1 (67–124 in LPC and aging), Acaca, Fdft1 and Fa2h (aging), and lipid transport / storage is the enriched pathway in LPC / cuprizone, Serpina3n-cKO cuprizone and 5XFAD; Apoe's top ranks in the demyelination datasets are microglial ambient RNA in lesion nuclei (see OligoC4b). OXPHOS genes are mildly anti-correlated with C4b in the aging and 5XFAD datasets. C4b-high versus C4b-negative oligodendrocytes have a **higher lipid-transport / storage score in all seven mouse datasets** (cell-level p < 0.001 everywhere; paired across samples p < 0.05 in LPC / cuprizone and 5XFAD), a lower ketone-utilisation score in aged white matter and after LPC / cuprizone, and lower glycolysis in 5XFAD. Together with the in-house data this makes lipid uptake and storage (Apod, Plin2 / Plin4, Abca1, Lpl) the one metabolic feature that travels with the C4b⁺ state across aging, amyloid, toxic demyelination and EAE, while energy metabolism (glycolysis, TCA, OXPHOS) is unchanged or slightly lower and ketone-body handling is not part of the program.
#
# **Caveats.** Nuclei detect metabolic transcripts (many of them cytoplasmic and mitochondrial) far less than whole cells, so mouse droplet datasets (Park, Ximerakis, Kaya) and human snRNA-seq should not be compared on detection rates; pseudobulk tests with fewer than three samples per group (Park, LPC / cuprizone, Serpina3n-cKO, several human contrasts) are descriptive; the Kaya comparison is white vs grey matter in aged brain; marker-based cell-type labels were validated only for oligodendrocytes and microglia; microglial ambient RNA in demyelinating-lesion nuclei inflates Apoe, Lpl and Hcar2 in "oligodendrocytes" from LPC / cuprizone; and pathway scores depend on set composition, so gene-level tables should be checked before drawing conclusions from a score.
