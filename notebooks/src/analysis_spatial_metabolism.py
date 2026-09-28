# %% [markdown]
# # Metabolism of oligodendrocytes in the in-house spatial datasets: ketone bodies, glycolysis / lactate, mitochondrial oxidation, fatty acids and lipids
#
# This notebook does for the metabolic machinery what `OligoC4b/analysis_spatial_complement_C5_C1q_Cfb.ipynb` did for the complement system. It asks,
# dataset by dataset:
#
# 1. **Is the gene measurable at all?** Xenium is targeted, so the first result is which metabolic genes are on each panel.
# 2. **Which cell types express which pathway, and does that change with disease / age, in every cell type?** Pathway scores (`scanpy.tl.score_genes`, set mean
#    minus a size-matched control set, log scale) and single genes for all annotated cell types; pseudobulk per sample wherever replicates exist; along the
#    lesion-distance gradient in EAE for each cell type.
# 3. **As one axis among these, how does metabolism relate to the C4b⁺ oligodendrocyte state?** Cell-level co-expression with C4b inside oligodendrocytes, the rank of every
#    metabolic gene among all panel genes correlated with C4b, a pathway-level rank-enrichment test, and pathway scores in disease-associated versus
#    homeostatic oligodendrocytes.
# 4. **Spatially**, whether cells expressing the ketone receptor Hcar2, the astrocytic lactate exporter Slc16a3 (MCT4) or the hypoxia-inducible hexokinase
#    Hk2 sit closer to C4b-high than to C4b-negative oligodendrocytes, with the lesion-distance-matched control from OligoC4b, and whether the
#    neighbourhood of C4b-high oligodendrocytes is itself more glycolytic.
#
# Datasets: Xenium mouse AD (TgCRND8 time course, 347-gene panel), Xenium mouse EAE (5K panel, 107 samples), Visium aging mouse brain (whole
# transcriptome), and the Falcão et al. 2018 EAE scRNA-seq (Smart-seq2, sorted oligodendrocyte-lineage cells) as the intrinsic-expression check.
# The gene panel and all helpers live in `scripts/oligometab.py` (183 genes in 21 pathway groups; see `docs/gene_panel.md`).

# %%
import os, sys, gc, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "../../scripts")
try:
    from dotenv import load_dotenv
    load_dotenv("../../.env")
except Exception:
    pass
import numpy as np, pandas as pd, scipy.sparse as sp, h5py, scanpy as sc, seaborn as sns, matplotlib.pyplot as plt
from scipy import stats
import oligometab as om
from oligometab import PANEL, ALL_GENES, FOCUS, SCORE_SETS, CONTEXT, PATHWAY_OF, heat

sc.settings.verbosity = 0
sc.set_figure_params(dpi=80, frameon=False)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 300)

PATHS = {
    "Xenium AD":     om.resolve_path("OLIGOMETAB_XENIUM_AD_H5AD",    "../../data/Xenium_AD_mouse.h5ad"),
    "Xenium EAE":    om.resolve_path("OLIGOMETAB_XENIUM_EAE_H5AD",   "../../data/RREAE_5k_raw_only_integration_processed.h5ad"),
    "Visium aging":  om.resolve_path("OLIGOMETAB_VISIUM_AGING_H5AD", "../../data/visum_aging_brain.h5ad"),
    "Falcão EAE sc": om.resolve_path("OLIGOMETAB_FALCAO_H5AD",       "../../data/falcao_et_al_2018.h5ad"),
}
for k, v in PATHS.items():
    print(f"{k:13s} -> {v if v else 'NOT FOUND (section will be skipped)'}")

FOCUS_ENERGY = ["Hmgcs2", "Bdh1", "Oxct1", "Acat1", "Hcar2", "Slc16a1", "Slc16a7", "Slc16a3", "Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkp", "Aldoc", "Pkm",
                "Ldha", "Ldhb", "Pdk1", "Cs", "Sdha", "Cox4i1", "Atp5a1", "Ppargc1a"]
FOCUS_LIPID = ["Cpt1a", "Acadm", "Hadha", "Fasn", "Hmgcr", "Srebf2", "Plin2", "Apoe", "Hif1a", "Txnip", "Ddit4"]
availability, results = {}, {}
os.makedirs("../../results", exist_ok=True)


def availability_by_pathway(adata):
    r = om.resolve(adata, ALL_GENES)
    rows = []
    for pw, genes in PANEL.items():
        pres = [g for g in genes if r[g] is not None]
        rows.append({"pathway": pw, "n_panel": len(genes), "n_present": len(pres), "present": ", ".join(pres)})
    return pd.DataFrame(rows).set_index("pathway")


def dotplot(adata, genes, groupby, title=None, **kw):
    r = om.resolve(adata, genes)
    hv = [r[g] for g in genes if r[g] is not None]
    if not hv:
        print("no genes to plot"); return
    n_groups = adata.obs[groupby].nunique()
    sc.pl.dotplot(adata, var_names=hv, groupby=groupby, standard_scale="var", color_map="Reds", figsize=(1.0 + 0.32 * len(hv), 0.8 + 0.26 * n_groups),
                  show=False, title=title, **kw)
    plt.show()


def score_by(adata, P, key, min_cells=30):
    g = adata.obs[key].astype(str).values
    s = P.groupby(g).mean()
    n = pd.Series(g).value_counts()
    return s.loc[n[n >= min_cells].index.intersection(s.index)]


def paired_pathway_test(adata, P, sample_col, mask_a, mask_b, label_a, label_b, min_cells=20):
    """Per sample: mean pathway score in cells of mask_a vs mask_b; paired Wilcoxon across samples."""
    smp = adata.obs[sample_col].astype(str).values
    rows = []
    for pw in P.columns:
        v = P[pw].values
        pairs = []
        for s in np.unique(smp):
            m = smp == s
            if (m & mask_a).sum() >= min_cells and (m & mask_b).sum() >= min_cells:
                pairs.append((v[m & mask_a].mean(), v[m & mask_b].mean()))
        ps = np.array(pairs)
        p = stats.wilcoxon(ps[:, 0], ps[:, 1]).pvalue if len(ps) >= 3 else np.nan
        rows.append({"pathway": pw, f"mean_{label_a}": v[mask_a].mean(), f"mean_{label_b}": v[mask_b].mean(), "diff_cells": v[mask_a].mean() - v[mask_b].mean(),
                     "n_samples_paired": len(ps), "median_paired_diff": np.median(ps[:, 0] - ps[:, 1]) if len(ps) else np.nan,
                     "frac_samples_positive": float((ps[:, 0] > ps[:, 1]).mean()) if len(ps) else np.nan, "wilcoxon_p": p})
    return pd.DataFrame(rows).set_index("pathway")


def paired_gene_test(adata, genes, sample_col, mask_a, mask_b, label_a, label_b, min_cells=20):
    """Per sample pseudobulk means of each gene in two cell populations; paired Wilcoxon across samples; log2FC of the mean of means."""
    df = om.expr_df(adata, genes)
    smp = adata.obs[sample_col].astype(str).values
    A, B = [], []
    for s in np.unique(smp):
        m = smp == s
        if (m & mask_a).sum() >= min_cells and (m & mask_b).sum() >= min_cells:
            A.append(df[m & mask_a].mean()); B.append(df[m & mask_b].mean())
    A, B = pd.DataFrame(A), pd.DataFrame(B)
    rows = []
    for g in df.columns:
        p = stats.wilcoxon(A[g], B[g]).pvalue if len(A) >= 3 and not np.allclose(A[g], B[g]) else np.nan
        rows.append({"gene": g, "pathway": PATHWAY_OF.get(g, "context"), f"mean_{label_a}": A[g].mean(), f"mean_{label_b}": B[g].mean(),
                     "log2FC": om.log2fc(A[g].mean(), B[g].mean()), "n_samples_paired": len(A), "wilcoxon_p": p})
    return pd.DataFrame(rows).set_index("gene")


def bin_order(values):
    """Sort lesion-distance bins like '0–10µm', '10–25µm', ... '>500µm' by their leading number."""
    import re
    def key(b):
        m = re.search(r"(\d+)", str(b))
        return (1, 0) if str(b).startswith(">") else ((0, int(m.group(1))) if m else (2, 0))
    return sorted([v for v in values if str(v) not in ("nan", "None")], key=key)


def spatial_panel(adata, sample_col, sample, genes, spot_size=15, vmax="p99", title=None, extra_obs=None):
    r = om.resolve(adata, genes)
    hv = [r[g] for g in genes if r[g] is not None]
    sub = adata[adata.obs[sample_col].astype(str) == str(sample)].copy()
    if extra_obs is not None:
        for c in extra_obs.columns:
            sub.obs[c] = extra_obs.loc[sub.obs_names, c].values
    sc.pl.spatial(sub, color=hv, spot_size=spot_size, cmap="magma", vmax=vmax, ncols=len(hv), show=False)
    if title:
        plt.suptitle(title, y=1.02)
    plt.show()
    if extra_obs is not None:
        sc.pl.spatial(sub, color=list(extra_obs.columns), spot_size=spot_size, cmap="RdBu_r", vcenter=0, ncols=len(extra_obs.columns), show=False)
        plt.suptitle(f"{title} — pathway scores", y=1.02); plt.show()

# %% [markdown]
# ---
# ## 1. Xenium mouse AD (TgCRND8 time course, 347-gene panel)
#
# The 247-gene brain panel plus the 99-gene add-on was designed around cell identity and neuroinflammation, not metabolism. Panel content is the result here.
# The matrix is `normalize_total(target_sum=100)` + `log1p` from `build_Xenium_AD_mouse.ipynb`.

# %%
ds = "Xenium AD"
ad_ad = None
if PATHS[ds]:
    ad_ad = sc.read_h5ad(PATHS[ds])
    ad_ad.obs.index = ad_ad.obs.index.astype(str); ad_ad.obs_names_make_unique()
    ad_ad.obs["model_month"] = ad_ad.obs["model"].astype(str) + "-" + ad_ad.obs["age_months"].astype(str) + "m"
    order = (ad_ad.obs[["model_month", "model", "age_months"]].drop_duplicates().sort_values(["model", "age_months"], ascending=[False, True])["model_month"].tolist())
    ad_ad.obs["model_month"] = pd.Categorical(ad_ad.obs["model_month"], categories=order)
    availability[ds] = {g: (v is not None) for g, v in om.resolve(ad_ad, ALL_GENES).items()}
    av = availability_by_pathway(ad_ad)
    print(ad_ad.shape); display(av[["n_panel", "n_present", "present"]])
    genes_ad = om.present(ad_ad, ALL_GENES)
    print("metabolic genes on the panel:", genes_ad)
else:
    print("skipped")

# %% [markdown]
# ### 1.1 The three measurable genes (Apoe, Apod, Acsbg1) by cell type and over the disease course
#
# Apod (lipid-binding, a known correlate of the C4b⁺ program), Apoe (cholesterol transport) and Acsbg1 (very-long-chain acyl-CoA synthetase, oligodendrocyte-
# enriched) are the only metabolic panel genes present. One section per model × age, so trends are descriptive.

# %%
if ad_ad is not None and genes_ad:
    dotplot(ad_ad, genes_ad + ["C4b", "Serpina3n", "Plp1"], "cellType", title="Xenium AD: metabolic panel genes by cell type")
    summ = om.detection_table(ad_ad, genes_ad + ["C4b"], "cellType")
    results[ds] = {"by_celltype": summ}
    print("fraction of cells with ≥1 transcript:"); display(summ.round(3))
    rows = []
    for ct in ["Oligodendrocytes", "Microglia", "Astrocytes"]:
        sub = ad_ad[ad_ad.obs["cellType"] == ct]
        df = om.expr_df(sub, genes_ad + ["C4b"]).join(sub.obs[["model", "age_months"]])
        m = df.groupby(["model", "age_months"], observed=True)[genes_ad + ["C4b"]].mean().reset_index(); m["cellType"] = ct; rows.append(m)
    course = pd.concat(rows); results[ds]["course"] = course
    tidy = course.melt(id_vars=["model", "age_months", "cellType"], var_name="gene", value_name="mean_expr")
    g = sns.relplot(data=tidy, x="age_months", y="mean_expr", hue="model", style="model", col="gene", row="cellType", kind="line", marker="o",
                    facet_kws={"sharey": False}, height=2.2, aspect=1.1)
    g.set_titles("{row_name} | {col_name}"); plt.show()
    # co-expression with C4b inside oligodendrocytes and rank among the 347 panel genes
    ol = ad_ad[ad_ad.obs["cellType"] == "Oligodendrocytes"].copy()
    co = om.coexpression_with_anchor(ol, "C4b", genes_ad + ["Serpina3n"])
    rho = om.genome_wide_rho(ol, "C4b", min_frac=0.005, max_cells=60000)
    rk = om.rank_panel(rho, ol, genes_ad + ["Serpina3n"])
    results[ds]["coexpr_oligo"] = co.join(rk[["rank", "n_ranked"]])
    print("co-expression with C4b in oligodendrocytes (rank among panel genes; 1 = most C4b-correlated):"); display(results[ds]["coexpr_oligo"].round(4))
    print("top 15 C4b-correlated panel genes:", ", ".join(rho.head(15).index))
    del ol

# %% [markdown]
# ---
# ## 2. Xenium mouse EAE (5K pan-tissue panel; RR-EAE and chronic EAE, 107 samples)
#
# The 5K panel carries about 100 of the 183 panel genes, including Slc16a1 (MCT1), Slc16a3 (MCT4), Hk1/Hk2, Pfkp, Pkm, Ldha/Ldhb, Pdk1, Bdh1, Hcar2, Cpt1a,
# Fasn, Hmgcr, Srebf1/2, Plin2 and the regulators. Missing are ketogenesis (Hmgcs2, Hmgcl), Oxct1/Acat1, most OXPHOS subunits, Apoe and most of the
# cholesterol pathway. The matrix is raw counts and is normalised (`target_sum=1e4`) and log1p-transformed here as in OligoC4b; the sparse matrix is read with
# `h5py` to skip the neighbour graphs.

# %%
ds = "Xenium EAE"
ad_eae = None
if PATHS[ds]:
    from anndata import AnnData
    from anndata.experimental import read_elem
    with h5py.File(PATHS[ds], "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        shape = tuple(f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X/data"][:].astype(np.float32), f["X/indices"][:], f["X/indptr"][:]), shape=shape)
        spatial = f["obsm/spatial"][:]
    obs.index = obs.index.astype(str); var.index = var.index.astype(str)
    ad_eae = AnnData(X=X, obs=obs, var=var); ad_eae.obsm["spatial"] = spatial; ad_eae.obs_names_make_unique(); del X
    ad_eae.obs["cell_type"] = ad_eae.obs["cell_type"].astype(str).replace({"DA oligodendrocytes": "DA Oligodendrocytes", "Astrocyte": "Astrocytes"})
    print(ad_eae.shape, "| raw counts? max =", ad_eae.X.max())
    sc.pp.normalize_total(ad_eae, target_sum=1e4); sc.pp.log1p(ad_eae)
    availability[ds] = {g: (v is not None) for g, v in om.resolve(ad_eae, ALL_GENES).items()}
    av = availability_by_pathway(ad_eae); display(av[["n_panel", "n_present", "present"]])
    genes_eae = om.present(ad_eae, ALL_GENES)
    print(ad_eae.obs["condition"].value_counts().to_dict(), ad_eae.obs["model"].value_counts().to_dict(), "| samples:", ad_eae.obs["sample_name"].nunique())
    P_eae = om.pathway_scores(ad_eae, SCORE_SETS, min_genes=3)
    print("pathway scores available on this panel:", list(P_eae.columns))
    print("not scorable (<3 genes on panel):", [k for k in SCORE_SETS if k not in P_eae.columns])
else:
    print("skipped")

# %% [markdown]
# ### 2.1 Which cell types express which pathway?

# %%
if ad_eae is not None:
    focus_eae = [g for g in FOCUS if g in genes_eae]
    dotplot(ad_eae, focus_eae + ["C4b", "Serpina3n", "Mbp"], "cell_type", title="Xenium EAE: metabolic focus genes by cell type")
    summ = om.detection_table(ad_eae, focus_eae + ["C4b"], "cell_type")
    results[ds] = {"by_celltype_detection": summ}
    print("fraction of cells with ≥1 transcript:"); display(summ.round(3))
    s_ct = score_by(ad_eae, P_eae, "cell_type"); results[ds]["by_celltype_scores"] = s_ct
    heat(s_ct, "Xenium EAE: mean pathway score by cell type (log scale; 0 = same as matched control genes)", "pathway score", fmt=".2f", row_labels=list(s_ct.index))

# %% [markdown]
# ### 2.2 EAE vs control, per cell type (pseudobulk per sample)
#
# Each sample is collapsed to one mean per cell type (≥20 cells), then EAE and CONTROL samples are compared with a Mann–Whitney test: pathway scores by
# difference, genes by log2 fold change.

# %%
if ad_eae is not None:
    vc = ad_eae.obs["cell_type"].value_counts()
    cts = [c for c in vc.index if vc[c] >= 1500 and c != "unclear"]          # every annotated cell type with enough cells
    print(len(cts), "cell types tested:", cts)
    m_cts = ad_eae.obs["cell_type"].isin(cts).values
    sub = ad_eae[m_cts]
    n = sub.obs.groupby(["sample_name", "cell_type"], observed=True).size().rename("n_cells").reset_index()
    n[["sample_name", "cell_type"]] = n[["sample_name", "cell_type"]].astype(str)
    pb = om.pseudobulk(sub, genes_eae, "sample_name", ["cell_type", "condition", "model"]).merge(n, on=["sample_name", "cell_type"]).query("n_cells >= 20")
    pbs = om.pseudobulk(sub, [], "sample_name", ["cell_type", "condition", "model"], values=P_eae.loc[sub.obs_names]).merge(n, on=["sample_name", "cell_type"]).query("n_cells >= 20")
    rows, rows_s = [], []
    for ct in cts:
        p_ct, s_ct = pb[pb.cell_type == ct], pbs[pbs.cell_type == ct]
        for g in genes_eae:
            ma, mr, p, na, nr = om.group_test(p_ct, g, "condition", "EAE", "CONTROL")
            rows.append({"cell_type": ct, "gene": g, "pathway": PATHWAY_OF[g], "mean_EAE": ma, "mean_CONTROL": mr, "log2FC": om.log2fc(ma, mr), "MWU_p": p, "n_EAE": na, "n_CONTROL": nr})
        for pw in P_eae.columns:
            ma, mr, p, na, nr = om.group_test(s_ct, pw, "condition", "EAE", "CONTROL")
            rows_s.append({"cell_type": ct, "pathway": pw, "mean_EAE": ma, "mean_CONTROL": mr, "diff": ma - mr, "MWU_p": p, "n_EAE": na, "n_CONTROL": nr})
    de, de_s = pd.DataFrame(rows), pd.DataFrame(rows_s)
    results[ds]["EAE_vs_control_genes"], results[ds]["EAE_vs_control_scores"] = de, de_s
    d = de_s.pivot(index="cell_type", columns="pathway", values="diff").loc[cts, list(P_eae.columns)]
    d = d[d.notna().any(axis=1)]                      # cell types without control samples have nothing to show
    p = de_s.pivot(index="cell_type", columns="pathway", values="MWU_p").loc[d.index, d.columns]
    annot = d.round(2).astype(str) + np.where(p < 0.05, "*", "")
    heat(d, "Xenium EAE: pathway-score difference EAE − CONTROL per cell type (pseudobulk per sample; * MWU p < 0.05)", "score difference", annot=annot.values, vmin=-0.3, vmax=0.3, row_labels=list(d.index))
    print("cell types without enough CONTROL samples (descriptive only, not shown):", [c for c in cts if c not in d.index])
    for label, block in [("energy", FOCUS_ENERGY), ("lipid / regulators", FOCUS_LIPID)]:
        cols = [g for g in block if g in genes_eae]
        d = de.pivot(index="cell_type", columns="gene", values="log2FC").loc[cts, cols]; d = d[d.notna().any(axis=1)]
        p = de.pivot(index="cell_type", columns="gene", values="MWU_p").loc[d.index, cols]
        annot = d.round(1).astype(str) + np.where(p < 0.05, "*", "")
        heat(d, f"Xenium EAE: pseudobulk log2FC EAE vs CONTROL, {label} genes (* MWU p < 0.05)", "log2 fold change", annot=annot.values, vmin=-2, vmax=2, row_labels=list(d.index))
    sig = de[(de.MWU_p < 0.05) & (de.log2FC.abs() > 0.5)].sort_values(["cell_type", "log2FC"])
    print(f"{len(sig)} gene-level changes with p < 0.05 and |log2FC| > 0.5:"); display(sig.round(4).reset_index(drop=True))
    tidy = pbs.melt(id_vars=["sample_name", "cell_type", "condition", "model", "n_cells"], value_vars=list(P_eae.columns), var_name="pathway", value_name="score")
    g = sns.catplot(data=tidy, x="cell_type", y="score", hue="condition", col="pathway", col_wrap=2, kind="box", sharey=False, height=3.2, aspect=2.4, showfliers=False, order=cts)
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=45)
    plt.show()
    del sub

# %% [markdown]
# ### 2.3 Lesion distance (EAE samples only)
#
# Pathway scores and focus genes as a function of distance from lesions, within oligodendrocyte-lineage cells, myeloid cells, astrocytes and all cells.

# %%
if ad_eae is not None:
    m_eae = ((ad_eae.obs["condition"] == "EAE") & ad_eae.obs["lesion_distance_bin"].notna()).values
    bins = bin_order(ad_eae.obs.loc[m_eae, "lesion_distance_bin"].astype(str).unique())
    print("bins:", bins)
    groups = {"all cells": np.ones(ad_eae.n_obs, bool),
              "oligodendrocyte lineage": ad_eae.obs["cell_type"].isin(["Oligodendrocytes", "DA Oligodendrocytes"]).values,
              "myeloid": ad_eae.obs["cell_type"].isin(["Microglia", "Macrophages", "Foamy Mic_Mac", "Activate Mic_Mac 1", "Activated Mic_Mac 2", "Efflux Mic_Mac", "Myeloid cells", "Proliferating microglia"]).values,
              "astrocytes": ad_eae.obs["cell_type"].isin(["Astrocytes", "DA astrocytes"]).values}
    ld_rows = []
    for label, mask in groups.items():
        s = ad_eae[m_eae & mask]
        sb = score_by(s, P_eae.loc[s.obs_names], "lesion_distance_bin", min_cells=50).reindex([b for b in bins])
        results[ds][f"lesion_distance_scores_{label}"] = sb
        mt = om.mean_table(s, [g for g in FOCUS_ENERGY + FOCUS_LIPID if g in genes_eae] + ["C4b"], "lesion_distance_bin", min_cells=50).reindex(bins)
        results[ds][f"lesion_distance_genes_{label}"] = mt
        del s
        print(f"--- {label}: mean pathway score by lesion distance"); display(sb.round(3))
        t = sb.reset_index().melt(id_vars="index", var_name="pathway", value_name="score").rename(columns={"index": "bin"}); t["population"] = label; ld_rows.append(t)
    ld = pd.concat(ld_rows)
    g = sns.relplot(data=ld, x="bin", y="score", hue="population", col="pathway", col_wrap=4, kind="line", marker="o", height=2.4, aspect=1.2, facet_kws={"sharey": False})
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=45)
    g.set_titles("{col_name}"); plt.show()
    # every cell type separately: cell type x lesion-distance bin, one heatmap per pathway (value = mean score minus the >500 µm value of that cell type)
    s_all = ad_eae[m_eae]
    ct_bin = P_eae.loc[s_all.obs_names].groupby([s_all.obs["cell_type"].astype(str).values, s_all.obs["lesion_distance_bin"].astype(str).values]).mean()
    n_ct_bin = s_all.obs.groupby([s_all.obs["cell_type"].astype(str).values, s_all.obs["lesion_distance_bin"].astype(str).values]).size()
    ct_bin = ct_bin[n_ct_bin >= 50]
    results[ds]["lesion_distance_scores_by_celltype"] = ct_bin
    for pw in P_eae.columns:
        m = ct_bin[pw].unstack(1).reindex(columns=bins)
        m = m.loc[[c for c in cts if c in m.index]]
        far = m[">500µm"] if ">500µm" in m.columns else m.iloc[:, -1]
        heat(m.sub(far, axis=0), f"Xenium EAE, {pw}: mean pathway score by cell type and lesion distance, relative to the >500 µm value of the same cell type", "score minus far-from-lesion score", fmt=".2f", vmin=-0.4, vmax=0.4, row_labels=list(m.index))
    del s_all, ct_bin
    print("oligodendrocyte lineage: mean log-expression by lesion distance"); display(results[ds]["lesion_distance_genes_oligodendrocyte lineage"].round(3))

# %% [markdown]
# ### 2.4 Co-expression with C4b inside DA oligodendrocytes
#
# Detection of each metabolic gene in C4b⁺ vs C4b⁻ DA oligodendrocytes, the rank of every metabolic gene among the ~5,000 panel genes correlated with C4b, and a
# pathway-level test of whether metabolic pathways sit systematically high or low in that ranking (Mann–Whitney of the pathway genes' rho against all other genes).

# %%
if ad_eae is not None:
    dao = ad_eae[ad_eae.obs["cell_type"] == "DA Oligodendrocytes"].copy()
    print("DA oligodendrocytes:", dao.n_obs, f"({(om.expr(dao, 'C4b') > 0).mean():.1%} C4b+)")
    co = om.coexpression_with_anchor(dao, "C4b", genes_eae + ["Serpina3n"])
    rho = om.genome_wide_rho(dao, "C4b", min_frac=0.005, max_cells=40000)
    rk = om.rank_panel(rho, dao, genes_eae + ["Serpina3n"])
    en = om.pathway_rank_enrichment(rho, dao, PANEL, min_genes=3)
    results[ds]["coexpr_DAoligo"] = co.join(rk[["rank", "n_ranked"]]); results[ds]["C4b_pathway_enrichment_DAoligo"] = en
    print(f"pathway-level enrichment among {len(rho)} genes ranked by correlation with C4b:"); display(en.round(4))
    print("20 best-ranked metabolic genes:"); display(results[ds]["coexpr_DAoligo"].sort_values("rank").head(20)[["pathway", "frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p", "rank"]].round(4))
    print("20 worst-ranked (anti-correlated) metabolic genes:"); display(results[ds]["coexpr_DAoligo"].sort_values("rank").tail(20)[["pathway", "frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p", "rank"]].round(4))
    print("top 25 C4b-correlated panel genes (any gene):", ", ".join(rho.head(25).index))
    # same in all oligodendrocyte-lineage cells of EAE animals
    ol_eae = ad_eae[ad_eae.obs["cell_type"].isin(["Oligodendrocytes", "DA Oligodendrocytes"]) & (ad_eae.obs["condition"] == "EAE")].copy()
    rho2 = om.genome_wide_rho(ol_eae, "C4b", min_frac=0.005, max_cells=40000)
    en2 = om.pathway_rank_enrichment(rho2, ol_eae, PANEL, min_genes=3)
    results[ds]["C4b_pathway_enrichment_oligo_EAE"] = en2
    print("\nAll oligodendrocyte-lineage cells in EAE animals, pathway enrichment:"); display(en2.round(4))
    del dao, ol_eae

# %% [markdown]
# ### 2.5 Disease-associated versus homeostatic oligodendrocytes: the metabolic profile of the DA state
#
# Within each sample that has ≥20 cells of both types, the mean pathway score (and gene mean) of DA oligodendrocytes is compared with that of homeostatic
# oligodendrocytes; paired Wilcoxon across samples. Also C4b-high vs C4b-negative cells within the oligodendrocyte lineage.

# %%
if ad_eae is not None:
    is_da = (ad_eae.obs["cell_type"] == "DA Oligodendrocytes").values
    is_ho = (ad_eae.obs["cell_type"] == "Oligodendrocytes").values
    t = paired_pathway_test(ad_eae, P_eae, "sample_name", is_da, is_ho, "DA", "homeostatic")
    results[ds]["DA_vs_homeostatic_scores"] = t
    print("pathway scores, DA vs homeostatic oligodendrocytes:"); display(t.round(4))
    tg = paired_gene_test(ad_eae, genes_eae, "sample_name", is_da, is_ho, "DA", "homeostatic")
    results[ds]["DA_vs_homeostatic_genes"] = tg
    print("genes most increased in DA oligodendrocytes:"); display(tg.sort_values("log2FC", ascending=False).head(20).round(4))
    print("genes most decreased in DA oligodendrocytes:"); display(tg.sort_values("log2FC").head(20).round(4))
    # C4b-high vs C4b-negative within the lineage
    c4b = om.expr(ad_eae, "C4b"); is_ol = is_da | is_ho
    thr = np.quantile(c4b[is_ol & (c4b > 0)], 0.75)
    t2 = paired_pathway_test(ad_eae, P_eae, "sample_name", is_ol & (c4b >= thr), is_ol & (c4b == 0), "C4b_high", "C4b_neg")
    results[ds]["C4bhigh_vs_neg_scores"] = t2
    print("pathway scores, C4b-high vs C4b-negative oligodendrocyte-lineage cells:"); display(t2.round(4))
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    for ax, tt, ttl in [(axes[0], t, "DA − homeostatic"), (axes[1], t2, "C4b-high − C4b-neg")]:
        cols = ["tab:red" if (p < 0.05) else "lightgrey" for p in tt["wilcoxon_p"]]
        ax.barh(tt.index, tt["median_paired_diff"], color=cols); ax.axvline(0, c="k", lw=0.8); ax.set_title(f"Xenium EAE oligodendrocytes: {ttl}\n(median paired difference across samples; red = Wilcoxon p < 0.05)", fontsize=9)
    plt.tight_layout(); plt.show()

# %% [markdown]
# ### 2.6 Spatial neighbourhood: metabolic partners of C4b-high oligodendrocytes
#
# Same design as OligoC4b section 2.5: in EAE samples, C4b-high (top quartile of C4b⁺) versus C4b-negative oligodendrocyte-lineage cells as sources, all cells
# within 30 µm as neighbours. Targets: Hcar2⁺ cells (the ketone-body / niacin receptor, myeloid), Hcar2⁺ myeloid cells, Slc16a3⁺ (MCT4, lactate-exporting)
# astrocytes, Hk2⁺ cells and Pdk1⁺ cells (hypoxia / glycolytic switch markers), Plin2⁺ cells (lipid-droplet-laden). A ratio > 1 means enrichment around
# C4b-high oligodendrocytes. Then the lesion-distance-matched control (sources stratified by `lesion_distance_bin`) and the non-lesion-only control.
# Finally, the mean glycolysis / OXPHOS / lipid score of the *neighbours* of C4b-high vs C4b-negative oligodendrocytes (source cells excluded).

# %%
if ad_eae is not None:
    eae = ad_eae[ad_eae.obs["condition"] == "EAE"].copy()
    P_e = P_eae.loc[eae.obs_names]
    c4b = om.expr(eae, "C4b")
    is_ol = eae.obs["cell_type"].isin(["Oligodendrocytes", "DA Oligodendrocytes"]).values
    thr = np.quantile(c4b[is_ol & (c4b > 0)], 0.75)
    src_high, src_neg = is_ol & (c4b >= thr), is_ol & (c4b == 0)
    myeloid = eae.obs["cell_type"].isin(["Microglia", "Macrophages", "Foamy Mic_Mac", "Activate Mic_Mac 1", "Activated Mic_Mac 2", "Efflux Mic_Mac", "Myeloid cells", "Proliferating microglia"]).values
    astro = eae.obs["cell_type"].isin(["Astrocytes", "DA astrocytes"]).values
    targets = {}
    for g, lab, extra in [("Hcar2", "Hcar2+ (any cell)", None), ("Hcar2", "Hcar2+ myeloid", myeloid), ("Slc16a3", "Slc16a3+ astrocytes", astro),
                          ("Hk2", "Hk2+ (any cell)", None), ("Pdk1", "Pdk1+ (any cell)", None), ("Plin2", "Plin2+ (any cell)", None), ("Slc2a1", "Slc2a1+ (any cell)", None)]:
        if g in genes_eae:
            m = om.expr(eae, g) > 0
            targets[lab] = m if extra is None else (m & extra)
    print(f"C4b-high oligos: {src_high.sum()}, C4b-neg oligos: {src_neg.sum()}")
    nb_summary, per_target = [], {}
    for name, tgt in targets.items():
        print(f"\n--- target: {name} (n={tgt.sum()})")
        j = om.compare_neighbourhoods(eae, "sample_name", src_high, src_neg, tgt, "C4bhigh", "C4bneg", radius=30.0)
        per_target[name] = j
        nb_summary.append({"target": name, "n_samples": len(j), "median_ratio": j["ratio"].median(), "frac_samples_ratio_gt1": (j["ratio"] > 1).mean(),
                           "wilcoxon_p": stats.wilcoxon(j.iloc[:, 0], j.iloc[:, 1]).pvalue if len(j) >= 3 else np.nan})
    nb_summary = pd.DataFrame(nb_summary).set_index("target"); results[ds]["neighbourhood_summary"] = nb_summary
    display(nb_summary.round(4))
    long = pd.concat([per_target[k]["ratio"].rename(k) for k in targets], axis=1).melt(var_name="target", value_name="ratio")
    plt.figure(figsize=(8, 3.5)); sns.boxplot(data=long, x="target", y="ratio", showfliers=False, color="lightgrey"); sns.stripplot(data=long, x="target", y="ratio", size=3, color="k", alpha=0.5)
    plt.axhline(1, ls="--", c="r"); plt.ylabel("neighbour fraction ratio\n(C4b-high / C4b-neg oligos)"); plt.xlabel(""); plt.xticks(rotation=30, ha="right"); plt.tight_layout(); plt.show()

# %%
if ad_eae is not None:
    # lesion-distance-matched control and non-lesion-only control
    strata = eae.obs["lesion_distance_bin"].astype(str).values
    per_bin_rows, per_sample = [], {}
    for name, tgt in targets.items():
        ratios = []
        for b in bins:
            m = strata == b
            A = om.neighbourhood_fraction(eae, "sample_name", src_high & m, tgt, radius=30.0, min_source=20)
            B = om.neighbourhood_fraction(eae, "sample_name", src_neg & m, tgt, radius=30.0, min_source=20)
            j = A[["frac_target_neighbours"]].join(B[["frac_target_neighbours"]], lsuffix="_high", rsuffix="_neg", how="inner")
            j = j[j.iloc[:, 1] > 0]; j["ratio"] = j.iloc[:, 0] / j.iloc[:, 1]
            if len(j) >= 3:
                per_bin_rows.append({"target": name, "bin": b, "n_samples": len(j), "median_ratio": j["ratio"].median(), "frac_gt1": (j["ratio"] > 1).mean(),
                                     "wilcoxon_p": stats.wilcoxon(np.log(j["ratio"].clip(lower=1e-3))).pvalue})
            ratios += [(smp, b, r) for smp, r in j["ratio"].items()]
        rdf = pd.DataFrame(ratios, columns=["sample", "bin", "ratio"])
        per_sample[name] = rdf.groupby("sample")["ratio"].median()
    per_bin = pd.DataFrame(per_bin_rows).set_index(["target", "bin"]); results[ds]["neighbourhood_lesion_matched_per_bin"] = per_bin
    print("per lesion-distance bin:"); display(per_bin.round(3))
    overall = pd.DataFrame({k: {"n_samples": len(v), "median_of_per_sample_median_ratio": v.median(), "frac_samples_gt1": (v > 1).mean(),
                                "wilcoxon_p": stats.wilcoxon(np.log(v.clip(lower=1e-3))).pvalue if len(v) >= 3 else np.nan} for k, v in per_sample.items()}).T
    results[ds]["neighbourhood_lesion_matched_overall"] = overall
    print("\nstratified overall estimate:"); display(overall.round(4))
    nl = (eae.obs["lesion_density_call"].astype(str) == "non_lesion").values
    rows = []
    for name, tgt in targets.items():
        A = om.neighbourhood_fraction(eae, "sample_name", src_high & nl, tgt, radius=30.0, min_source=20)
        B = om.neighbourhood_fraction(eae, "sample_name", src_neg & nl, tgt, radius=30.0, min_source=20)
        j = A[["frac_target_neighbours"]].join(B[["frac_target_neighbours"]], lsuffix="_high", rsuffix="_neg", how="inner")
        j = j[j.iloc[:, 1] > 0]; j["ratio"] = j.iloc[:, 0] / j.iloc[:, 1]
        rows.append({"target": name, "n_samples": len(j), "median_ratio": j["ratio"].median(), "frac_gt1": (j["ratio"] > 1).mean(),
                     "wilcoxon_p": stats.wilcoxon(np.log(j["ratio"].clip(lower=1e-3))).pvalue if len(j) >= 3 else np.nan})
    nonlesion = pd.DataFrame(rows).set_index("target"); results[ds]["neighbourhood_nonlesion_only"] = nonlesion
    print("\nnon-lesion tissue only:"); display(nonlesion.round(4))
    pb_plot = per_bin.reset_index()
    plt.figure(figsize=(8, 3.5))
    sns.pointplot(data=pb_plot, x="bin", y="median_ratio", hue="target", order=[b for b in bins if b in pb_plot["bin"].values], dodge=0.4, linestyle="none", errorbar=None)
    plt.axhline(1, ls="--", c="r"); plt.ylabel("median ratio (C4b-high / C4b-neg)"); plt.xlabel("distance of source oligodendrocytes from lesion")
    plt.xticks(rotation=30, ha="right"); plt.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left"); plt.tight_layout(); plt.show()

# %%
if ad_eae is not None:
    # is the neighbourhood of C4b-high oligodendrocytes itself more glycolytic / less oxidative? mean pathway score of neighbours (sources excluded)
    rows = []
    for pw in P_e.columns:
        vals = P_e[pw].values.copy()
        vals[is_ol] = np.nan          # exclude all oligodendrocyte-lineage cells so we read the non-oligodendrocyte micro-environment
        A = om.neighbourhood_mean(eae, "sample_name", src_high, vals, radius=30.0, min_source=30)
        B = om.neighbourhood_mean(eae, "sample_name", src_neg, vals, radius=30.0, min_source=30)
        j = A[["mean_neighbour_value"]].join(B[["mean_neighbour_value"]], lsuffix="_high", rsuffix="_neg", how="inner")
        d = j.iloc[:, 0] - j.iloc[:, 1]
        rows.append({"pathway": pw, "n_samples": len(j), "median_diff_high_minus_neg": d.median(), "frac_samples_positive": (d > 0).mean(),
                     "wilcoxon_p": stats.wilcoxon(j.iloc[:, 0], j.iloc[:, 1]).pvalue if len(j) >= 3 else np.nan})
    nbm = pd.DataFrame(rows).set_index("pathway"); results[ds]["neighbourhood_pathway_scores"] = nbm
    print("mean pathway score of non-oligodendrocyte neighbours (30 µm) around C4b-high vs C4b-neg oligodendrocytes:"); display(nbm.round(4))
    del eae, P_e

# %% [markdown]
# ### 2.7 Spatial maps (one peak-EAE sample and one control sample)

# %%
if ad_eae is not None:
    meta = ad_eae.obs.groupby("sample_name", observed=True).agg(condition=("condition", "first"), course=("course", "first"), model=("model", "first"), n=("condition", "size"))
    peak = meta[(meta.condition == "EAE") & meta.course.astype(str).str.contains("peak", case=False)].sort_values("n", ascending=False).index[0]
    ctrl = meta[meta.condition == "CONTROL"].sort_values("n", ascending=False).index[0]
    show_scores = P_eae[[c for c in ["Glycolysis", "Pyruvate to lactate", "Beta-oxidation", "Lipid transport / storage"] if c in P_eae.columns]]
    for s in [peak, ctrl]:
        spatial_panel(ad_eae, "sample_name", s, ["C4b", "Slc16a1", "Ldha", "Hk2", "Hcar2", "Cpt1a", "Plin2"], spot_size=20,
                      title=f"{s} ({meta.loc[s, 'condition']}, {meta.loc[s, 'course']})", extra_obs=show_scores)

# %%
if ad_eae is not None:
    for k, v in results[ds].items():
        if isinstance(v, pd.DataFrame):
            v.to_csv(f"../../results/xenium_eae_{k.replace(' ', '_').replace('/', '-')}.csv")
    del ad_eae, P_eae, show_scores, meta; gc.collect()

# %% [markdown]
# ---
# ## 3. Visium aging mouse brain (whole transcriptome)
#
# Every panel gene is measurable. Spots are ~55 µm mixtures, so cell-type attribution is indirect; a white-matter score (Plp1, Mbp, Mobp, Mag, Cldn11) selects
# oligodendrocyte-rich spots. Age mapping as in OligoC4b (Young = 6, Mid = 18, Old = 21 months; two sections per group).

# %%
ds = "Visium aging"
ad_vis = None
if PATHS[ds]:
    ad_vis = sc.read_h5ad(PATHS[ds]); ad_vis.var_names_make_unique()
    ad_vis.obs["age_months"] = ad_vis.obs["age_group"].map({"Young": 6, "Mid": 18, "Old": 21}).astype(float)
    ad_vis.obs["age_group"] = pd.Categorical(ad_vis.obs["age_group"].astype(str), categories=["Young", "Mid", "Old"])
    availability[ds] = {g: (v is not None) for g, v in om.resolve(ad_vis, ALL_GENES).items()}
    print(ad_vis.shape, "| missing:", om.missing(ad_vis, ALL_GENES))
    P_vis = om.pathway_scores(ad_vis, SCORE_SETS)
    wm_genes = om.present(ad_vis, ["Plp1", "Mbp", "Mobp", "Mag", "Cldn11"])
    sc.tl.score_genes(ad_vis, wm_genes, score_name="wm_score", use_raw=False)
    ad_vis.obs["wm_rich"] = ad_vis.obs["wm_score"] >= ad_vis.obs["wm_score"].quantile(0.75)
    print("white-matter-rich spots:", int(ad_vis.obs["wm_rich"].sum()), "of", ad_vis.n_obs)
    results[ds] = {}
else:
    print("skipped")

# %% [markdown]
# ### 3.1 Age dependence (pseudobulk per section, n = 6)

# %%
if ad_vis is not None:
    dotplot(ad_vis, FOCUS + ["C4b", "Serpina3n"], "age_group", title="Visium aging: focus genes by age group")
    genes_v = om.present(ad_vis, ALL_GENES)
    for label, wm in [("all spots", None), ("white-matter-rich spots", True)]:
        a_ = ad_vis if wm is None else ad_vis[ad_vis.obs["wm_rich"].values]
        pb = om.pseudobulk(a_, genes_v, "sample", ["age_group"]); pb["age_months"] = pb["age_group"].map({"Young": 6, "Mid": 18, "Old": 21}).astype(float)
        pbs = om.pseudobulk(a_, [], "sample", ["age_group"], values=P_vis.loc[a_.obs_names]); pbs["age_months"] = pbs["age_group"].map({"Young": 6, "Mid": 18, "Old": 21}).astype(float)
        rows = []
        for g in genes_v:
            r = stats.linregress(pb["age_months"], pb[g])
            rows.append({"gene": g, "pathway": PATHWAY_OF[g], "slope_per_month": r.slope, "pearson_r": r.rvalue, "p": r.pvalue, "mean_Young": pb.loc[pb.age_group == "Young", g].mean(),
                         "mean_Old": pb.loc[pb.age_group == "Old", g].mean(), "log2FC_old_vs_young": om.log2fc(pb.loc[pb.age_group == "Old", g].mean(), pb.loc[pb.age_group == "Young", g].mean())})
        age_df = pd.DataFrame(rows).set_index("gene"); results[ds][f"age_slope_genes_{label.split(' ')[0]}"] = age_df
        rows = []
        for pw in P_vis.columns:
            r = stats.linregress(pbs["age_months"], pbs[pw])
            rows.append({"pathway": pw, "slope_per_month": r.slope, "pearson_r": r.rvalue, "p": r.pvalue, "mean_Young": pbs.loc[pbs.age_group == "Young", pw].mean(), "mean_Old": pbs.loc[pbs.age_group == "Old", pw].mean()})
        age_s = pd.DataFrame(rows).set_index("pathway"); results[ds][f"age_slope_scores_{label.split(' ')[0]}"] = age_s
        print(f"=== {label}: pathway scores vs age (linear fit over 6 sections)"); display(age_s.round(4))
        print(f"{label}: focus genes vs age"); display(age_df.loc[[g for g in FOCUS if g in age_df.index]].round(4))
        print(f"{label}: panel genes with p < 0.05 for the age trend"); display(age_df[age_df.p < 0.05].sort_values("pearson_r").round(4))
        tidy = pbs.melt(id_vars=["sample", "age_group", "age_months"], value_vars=list(P_vis.columns), var_name="pathway", value_name="score")
        g = sns.catplot(data=tidy, x="age_group", y="score", col="pathway", col_wrap=6, kind="strip", sharey=False, height=2.2, size=7); g.fig.suptitle(label, y=1.03); plt.show()
        del a_

# %% [markdown]
# ### 3.2 Correlation with C4b across spots and within white-matter-rich spots

# %%
if ad_vis is not None:
    for label, sub in [("all spots", ad_vis), ("white-matter-rich spots", ad_vis[ad_vis.obs["wm_rich"].values].copy())]:
        co = om.coexpression_with_anchor(sub, "C4b", genes_v + ["Serpina3n"])
        rho = om.genome_wide_rho(sub, "C4b", min_frac=0.01, max_cells=60000)
        rk = om.rank_panel(rho, sub, genes_v + ["Serpina3n"]); en = om.pathway_rank_enrichment(rho, sub, PANEL)
        results[ds][f"coexpr_{label.split(' ')[0]}"] = co.join(rk[["rank", "n_ranked"]]); results[ds][f"C4b_pathway_enrichment_{label.split(' ')[0]}"] = en
        print(f"=== {label}: pathway enrichment among {len(rho)} genes ranked by correlation with C4b"); display(en.round(4))
        print("20 best-ranked metabolic genes:"); display(results[ds][f"coexpr_{label.split(' ')[0]}"].sort_values("rank").head(20)[["pathway", "spearman_rho", "rank", "frac_in_C4b-", "frac_in_C4b+"]].round(4))
        print("top 20 C4b-correlated genes (any):", ", ".join(rho.head(20).index))
        del sub
    # pathway scores in C4b-high vs C4b-negative white-matter spots (per section paired)
    wm = ad_vis[ad_vis.obs["wm_rich"].values].copy(); c4b = om.expr(wm, "C4b"); pos = c4b > 0
    thr = np.quantile(c4b[pos], 0.75)
    t = paired_pathway_test(wm, P_vis.loc[wm.obs_names], "sample", c4b >= thr, c4b == 0, "C4b_high", "C4b_neg", min_cells=10)
    results[ds]["C4bhigh_vs_neg_scores_wm"] = t
    print("white-matter-rich spots: pathway scores in C4b-high vs C4b-negative spots (paired per section):"); display(t.round(4))
    del wm

# %% [markdown]
# ### 3.3 Spatial maps (one section per age group)

# %%
if ad_vis is not None:
    picks = ad_vis.obs.groupby("age_group", observed=True)["sample"].first()
    show_scores = P_vis[["Glycolysis", "OXPHOS", "Ketone utilisation", "Cholesterol synthesis"]]
    for age, s in picks.items():
        spatial_panel(ad_vis, "sample", s, ["C4b", "Slc16a1", "Hmgcs2", "Ldha", "Cox4i1", "Apoe"], spot_size=90, title=f"{s} ({age})", extra_obs=show_scores)
    for k, v in results[ds].items():
        if isinstance(v, pd.DataFrame):
            v.to_csv(f"../../results/visium_aging_{k.replace(' ', '_')}.csv")

# %% [markdown]
# ---
# ## 4. Intrinsic expression: Falcão et al. 2018 EAE scRNA-seq (Smart-seq2, sorted cells)
#
# In situ, a transcript inside an oligodendrocyte segment can come from a neighbouring process; sorted single cells do not have that problem. The dataset has
# oligodendrocyte-lineage cells (OPC, COP, NFOL, MOL), microglia and VLMC from control and EAE spinal cord; `Renamed_clusternames` marks the disease-associated
# MOL clusters (`MOL*_EAE`). Cells are the replicates here (few animals), so p-values are descriptive. The matrix is already log-normalised.

# %%
ds = "Falcão EAE sc"
ad_f = None
if PATHS[ds]:
    ad_f = sc.read_h5ad(PATHS[ds]); ad_f.var_names_make_unique()
    cl = ad_f.obs["Renamed_clusternames"].astype(str)
    grp = np.select([cl.str.contains("MOL") & cl.str.contains("EAE"), cl.str.contains("MOL"), cl.str.startswith("MiGl"), cl.str.startswith("VLMC"),
                     cl.isin(["OPC1", "OPC2", "OPC3", "OPC_Cycling", "COP", "NFOL", "PLC"])],
                    ["MOL (EAE clusters)", "MOL (control clusters)", "Microglia", "VLMC", "OPC/COP/NFOL"], default="other")
    ad_f.obs["group_coarse"] = pd.Categorical(grp, categories=["MOL (control clusters)", "MOL (EAE clusters)", "OPC/COP/NFOL", "Microglia", "VLMC", "other"])
    availability[ds] = {g: (v is not None) for g, v in om.resolve(ad_f, ALL_GENES).items()}
    print(ad_f.shape, ad_f.obs["Group"].value_counts().to_dict(), "| missing:", om.missing(ad_f, ALL_GENES))
    print(ad_f.obs["group_coarse"].value_counts().to_dict())
    P_f = om.pathway_scores(ad_f, SCORE_SETS)
    results[ds] = {}
else:
    print("skipped")

# %%
if ad_f is not None:
    genes_f = om.present(ad_f, ALL_GENES)
    dotplot(ad_f, FOCUS + ["C4b", "Serpina3n", "Plp1", "Hexb"], "Renamed_clusternames", title="Falcão EAE scRNA-seq: focus genes by cluster")
    s_cl = score_by(ad_f, P_f, "Renamed_clusternames", min_cells=8); results[ds]["scores_by_cluster"] = s_cl
    heat(s_cl, "Falcão: mean pathway score by cluster (lineage stages, EAE clusters, microglia)", "pathway score", fmt=".2f", row_labels=list(s_cl.index))
    s_g = score_by(ad_f, P_f, "group_coarse", min_cells=10); results[ds]["scores_by_group"] = s_g
    display(s_g.round(3))
    # MOL EAE vs control clusters (cells as replicates; descriptive)
    m_eae = (ad_f.obs["group_coarse"] == "MOL (EAE clusters)").values; m_ctl = (ad_f.obs["group_coarse"] == "MOL (control clusters)").values
    rows = []
    for pw in P_f.columns:
        v = P_f[pw].values
        rows.append({"pathway": pw, "mean_MOL_EAE": v[m_eae].mean(), "mean_MOL_ctrl": v[m_ctl].mean(), "diff": v[m_eae].mean() - v[m_ctl].mean(), "MWU_p_cells": stats.mannwhitneyu(v[m_eae], v[m_ctl]).pvalue})
    t = pd.DataFrame(rows).set_index("pathway"); results[ds]["MOL_EAE_vs_ctrl_scores"] = t
    print("pathway scores, MOL EAE clusters vs MOL control clusters:"); display(t.round(4))
    df = om.expr_df(ad_f, genes_f)
    rows = []
    for g in genes_f:
        x = df[g].values
        rows.append({"gene": g, "pathway": PATHWAY_OF[g], "frac_MOL_EAE": (x[m_eae] > 0).mean(), "frac_MOL_ctrl": (x[m_ctl] > 0).mean(), "mean_MOL_EAE": x[m_eae].mean(), "mean_MOL_ctrl": x[m_ctl].mean(),
                     "log2FC": om.log2fc(x[m_eae].mean(), x[m_ctl].mean()), "MWU_p_cells": stats.mannwhitneyu(x[m_eae], x[m_ctl]).pvalue if x.std() > 0 else np.nan})
    tg = pd.DataFrame(rows).set_index("gene"); results[ds]["MOL_EAE_vs_ctrl_genes"] = tg
    print("genes most increased in EAE MOL clusters:"); display(tg.sort_values("log2FC", ascending=False).head(20).round(4))
    print("genes most decreased in EAE MOL clusters:"); display(tg.sort_values("log2FC").head(20).round(4))

# %%
if ad_f is not None:
    mol = ad_f[ad_f.obs["group_coarse"].isin(["MOL (control clusters)", "MOL (EAE clusters)"])].copy()
    n_pos = int((om.expr(mol, "C4b") > 0).sum())
    print(f"MOL cells: {mol.n_obs}; C4b+: {n_pos}")
    co = om.coexpression_with_anchor(mol, "C4b", genes_f + ["Serpina3n"])
    rho = om.genome_wide_rho(mol, "C4b", min_frac=0.02, max_cells=60000)
    rk = om.rank_panel(rho, mol, genes_f + ["Serpina3n"]); en = om.pathway_rank_enrichment(rho, mol, PANEL)
    results[ds]["coexpr_MOL"] = co.join(rk[["rank", "n_ranked"]]); results[ds]["C4b_pathway_enrichment_MOL"] = en
    print(f"pathway enrichment among {len(rho)} genes ranked by correlation with C4b in MOL:"); display(en.round(4))
    print("20 best-ranked metabolic genes:"); display(results[ds]["coexpr_MOL"].sort_values("rank").head(20)[["pathway", "frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p", "rank"]].round(4))
    print("20 worst-ranked metabolic genes:"); display(results[ds]["coexpr_MOL"].sort_values("rank").tail(20)[["pathway", "frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p", "rank"]].round(4))
    print("top 25 C4b-correlated genes (any):", ", ".join(rho.head(25).index))
    for k, v in results[ds].items():
        if isinstance(v, pd.DataFrame):
            v.to_csv(f"../../results/falcao_{k.replace(' ', '_')}.csv")
    del mol

# %% [markdown]
# ---
# ## 5. Cross-dataset availability and summary

# %%
avail = pd.DataFrame(availability).reindex(ALL_GENES); avail.index.name = "gene"
avail.insert(0, "pathway", [PATHWAY_OF[g] for g in avail.index])
summary = avail.groupby("pathway").agg(["sum", "size"]) if False else pd.concat({c: avail.groupby("pathway")[c].sum() for c in avail.columns[1:]}, axis=1)
summary.insert(0, "n_genes", avail.groupby("pathway").size())
print("panel genes measurable per pathway and dataset:"); display(summary.loc[list(PANEL)])
display(avail.replace({True: "✓", False: "✗"}))
avail.to_csv("../../results/inhouse_panel_availability.csv")
if "Xenium AD" in results and "coexpr_oligo" in results["Xenium AD"]:
    results["Xenium AD"]["coexpr_oligo"].to_csv("../../results/xenium_ad_coexpr_oligo.csv")

# %% [markdown]
# ### Interpretation
#
# **Panel content first.** The 347-gene Xenium AD panel carries three metabolic genes (Apoe, Apod, Acsbg1), so metabolism cannot be read from that dataset; Apod is the only one that tracks C4b in oligodendrocytes (rank 7 of 341 panel genes, in line with OligoC4b). The 5K EAE panel has 100 of 183 panel genes (no ketogenesis enzymes, no Oxct1/Acat1, two OXPHOS subunits, three cholesterol-synthesis genes), Visium and Falcão have everything.
#
# **Baseline division of labour (Xenium EAE, 27 cell types; Falcão lineage stages).** Astrocytes and neurons carry the highest glycolysis and TCA scores; homeostatic oligodendrocytes the highest lipid-synthesis score and newly formed oligodendrocytes the highest cholesterol-synthesis score; microglia, macrophages and especially foamy myeloid cells the highest lipid-transport / storage score together with the lowest cholesterol synthesis. In sorted cells, OPC/COP/NFOL have the highest ketone-utilisation and TCA scores of the lineage, mature oligodendrocytes the highest cholesterol- and myelin-lipid synthesis, microglia the highest glycolysis and pentose-phosphate scores.
#
# **EAE is a tissue-wide metabolic shift, not an oligodendrocyte-specific one.** Comparing 92 EAE with 15 control samples per cell type, cholesterol synthesis falls and lipid transport / storage rises in essentially every cell type with enough control samples (oligodendrocytes, DA oligodendrocytes, OPCs, newly formed oligodendrocytes, astrocytes, DA astrocytes, microglia, activated myeloid cells, neurons, endothelial, vascular, fibroblasts, stromal, ependymal and Schwann cells), and the TCA score falls in most of them. Glycolysis rises only in microglia and DA oligodendrocytes and falls in OPCs, neurons and endothelial cells. At the gene level Hcar2, Cd36, Hk2, Plin2, Plin4, Abca1, Lpl, Slc16a3 and Txnip rise across many cell types; because this is in situ, the myeloid genes among them (Hcar2, Cd36, Slc16a3) can be spill-over from adjacent myeloid processes and are checked below in sorted cells.
#
# **Lesion distance.** Toward lesions, glycolysis and TCA scores fall in neurons, homeostatic oligodendrocytes, endothelial and vascular cells, while lipid transport / storage rises in nearly every cell type (largest in endothelial cells, microglia, OPCs and oligodendrocytes). Cholesterol synthesis falls most steeply in OPCs (−0.33), newly formed oligodendrocytes (−0.29) and DA oligodendrocytes (−0.16), the cells that would need it to remyelinate, and rises in homeostatic oligodendrocytes far from the core. Myeloid cells are glycolytic everywhere and become more lactate-producing (Pyruvate-to-lactate score) toward lesions.
#
# **The disease-associated oligodendrocyte state (paired within 107 samples).** DA oligodendrocytes score lower than homeostatic ones on glycolysis (−0.10), TCA (−0.12), lipid synthesis (−0.15) and cholesterol synthesis (−0.16) and higher on lipid transport / storage (+0.15; 99 % of samples), with β-oxidation unchanged. Genes: Hk2 (3.5×), Ldha, Pdk1, Slc16a3, Slc2a4, Plin2, Abca1, Lpl, Cd36, Hcar2, Mlxipl, Txnip and Fabp7 up; Slc2a3, Eno2, Ppargc1a, Mpc2, Pcx, Got2, Fh1, Ogdh, Mfn2, Dnm1l, Elovl1 and Pfkp down. C4b-high versus C4b-negative lineage cells give the same picture. **Sorted cells (Falcão) confirm the core of this intrinsically:** EAE MOL clusters versus control MOL clusters show cholesterol synthesis −1.0 (Hmgcs1, Fdps, Cyp51, Dhcr7, Dhcr24, Sqle each 20–30 % lower), myelin-lipid synthesis −0.30, lipid transport / storage +0.95 (Plin4 5×, Plin2 3×, Abca1 2.5×), Pdk4 4×, Ldha 2.2×, Cpt1a 2.5×, Bdh1 2.3×, and Slc16a1 (MCT1) −0.24 log2. Where the two datasets disagree is the direction of the glycolysis / TCA / β-oxidation scores (down in Xenium DA oligodendrocytes, up in sorted EAE MOL); the Xenium glycolysis set is 15 genes weighted toward Eno2, Pfkp and Slc2a3, which fall, whereas the whole-transcriptome set is dominated by Aldoa, Tpi1, Gapdh and Ldha, which rise, so the gene-level tables are the safer read-out here. Ketone-body handling in oligodendrocytes is modest: Bdh1 rises, Oxct1 / Acat1 are flat or lower, Hmgcs2 is barely expressed; Hcar2 is a myeloid receptor (8 % of microglia, 0.3 % of oligodendrocytes in EAE) and its 2 % detection in EAE MOL is a handful of cells.
#
# **The C4b axis.** Inside DA oligodendrocytes the metabolic genes that follow C4b are the lipid-droplet / lipid-handling genes (Plin4 rank 2 of 3,673 panel genes, right behind C4a-type immune genes and Serpina3n; Apod 42; Cers2 72; Fasn 86); no metabolic pathway is enriched as a set, and across the whole oligodendrocyte lineage of EAE animals the glycolysis, TCA and mitochondrial-biogenesis sets are significantly *anti*-correlated with C4b. In sorted MOL the same holds with more power: lipid transport / storage is the top pathway (Apod 0.46, Plin4 0.40, Abca1 0.32), cholesterol synthesis is significantly anti-correlated, and glycolysis, TCA and β-oxidation are mildly positive. In Visium the spot-level C4b correlates are white-matter lipid genes (Apoe, Apod, Scd2, Fa2h, Ugt8a), i.e. tissue composition; within white-matter-rich spots, C4b-high spots have higher lipid-synthesis / storage / myelin-lipid and pentose-phosphate scores and lower glycolysis / TCA / OXPHOS in all six sections.
#
# **Spatial partners.** In EAE, Hcar2⁺ cells (1.7×), Hcar2⁺ myeloid cells (1.8×), Slc16a3⁺ (MCT4, lactate-exporting) astrocytes (2.3×), Hk2⁺ (1.8×), Pdk1⁺ (1.3×) and Plin2⁺ (1.4×) cells are enriched within 30 µm of C4b-high compared with C4b-negative oligodendrocytes (paired Wilcoxon p < 10⁻⁵ in 87 samples), whereas Slc2a1⁺ cells are not (0.97). The enrichment persists after matching source oligodendrocytes for lesion distance (1.2–1.75) and in non-lesion tissue (1.3–2.2), so, as with C5aR1 in OligoC4b, C4b-high oligodendrocytes mark glycolytic, lactate-exporting, lipid-loaded micro-niches also outside lesions. The non-oligodendrocyte neighbourhood of C4b-high cells is itself less glycolytic and less oxidative on the score scale and richer in lipid storage, which mostly reflects who the neighbours are (myeloid rather than neuronal / astrocytic) rather than a per-cell change.
#
# **Aging (Visium, 6 sections).** Whole-section glycolysis, TCA and OXPHOS scores rise slightly with age (r ≈ 0.84–0.86, p ≈ 0.03–0.04, slopes of a few hundredths over 15 months), driven by housekeeping glycolytic and OXPHOS genes (Pkm, Aldoa, Gapdh, Hk1, Ldhb, Ndufv1, Uqcrc1, Sdhb). In white-matter-rich spots the lipid-synthesis score falls (p = 0.03) and the myelin-lipid and cholesterol-synthesis genes Fa2h, Ugt8a, Cyp51 and Hmgcs1 decline by 10–19 %, while Apoe (+0.12 log2), Plin2 (+0.48 log2) and Bdh1 rise and Ddit4 (−0.63 log2), Slc2a1 (−0.15 log2) and Sirt1 fall. Aging white matter therefore shows a milder version of the EAE pattern: less lipid / cholesterol synthesis, more lipid storage, with Hcar2 rising from a very low base.
#
# **Caveats.** One Xenium AD section per genotype × age; the EAE 5K panel misses the ketogenesis, ketolysis, most OXPHOS and most cholesterol-synthesis genes, so those scores are absent or thin there; in situ transcripts inside a segment can come from a neighbouring cell (Hcar2, Cd36, Slc16a3 in "oligodendrocytes"); the Visium age trend rests on 6 sections; Falcão cells are the replicates (few animals), so its p-values are descriptive; pathway scores are relative to size-matched control genes and are compared by difference, not fold change.
