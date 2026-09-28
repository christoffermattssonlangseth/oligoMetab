# %% [markdown]
# # Ketone-body metabolism across brain cell types, disease and the *C4b*⁺ oligodendrocyte state
#
# The two main notebooks treat ketone bodies as two pathway scores among twelve. This notebook takes the ketone arm apart gene by gene, in every
# annotated cell type, across the in-house datasets (Xenium EAE, Visium aging, Falcão sorted cells) and the sixteen public datasets. Xenium AD has no
# ketone gene on its panel and is skipped.
#
# **Biology in one paragraph.** Ketone bodies (acetoacetate, β-hydroxybutyrate) are made from fatty-acid-derived acetyl-CoA by *Hmgcs2* (rate-limiting),
# *Hmgcl* and *Bdh1*; in the brain this ketogenic capacity sits in astrocytes, which oxidise fatty acids (*Cpt1a*, *Acadm*, *Hadha*) and are driven by
# PPARα. They are used by *Bdh1* (BHB → acetoacetate), *Oxct1* (SCOT, the committed ketolytic step, absent from liver) and *Acat1*, and they cross
# membranes through the monocarboxylate transporters MCT1 (*Slc16a1*, oligodendrocytes and endothelium), MCT2 (*Slc16a7*, neurons), MCT4 (*Slc16a3*,
# astrocytes), MCT7 (*Slc16a6*) and SMCT1 (*Slc5a8*). BHB also signals: it is the ligand of *Hcar2* (GPR109A) on myeloid cells and an antagonist of the
# SCFA receptor *Ffar3*. Acetate is activated by *Acss1* / *Acss2*.
#
# Questions:
#
# 1. **Who makes, who uses, who transports and who senses ketone bodies?** Detection and mean expression of every ketone gene by cell type, plus the
#    fraction of cells that co-detect *Oxct1* and *Bdh1* (complete ketolytic machinery).
# 2. **How does that change with disease, age, demyelination and lesion distance, in every cell type?** Pseudobulk per sample, gene-wise and as
#    ketogenesis / ketolysis / fatty-acid-supply scores.
# 3. **Is astrocyte ketogenesis coupled to fatty-acid oxidation and PPARα?** Correlation of *Hmgcs2* with *Cpt1a*, *Acadm*, *Ppara* inside astrocytes.
# 4. **What does the disease-associated oligodendrocyte do with ketones?** *Bdh1*, *Oxct1*, *Acat1*, *Slc16a1* in DA vs homeostatic oligodendrocytes (in situ,
#    sorted), and the ketone genes' correlation with *C4b*.
# 5. **Human MS / AD** by lesion type and Braak stage.
#
# Gene set: `oligometab.KETONE` (26 genes in six groups); scores: `oligometab.KETONE_SCORE_SETS`. Pathway scores are `scanpy.tl.score_genes` values compared by
# difference; genes are compared by log2 fold change of pseudobulk means; Mann–Whitney across samples when both groups have ≥3 samples.

# %%
import os, sys, gc, traceback, warnings
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
from oligometab import KETONE, KETONE_GENES, KETONE_SCORE_SETS, heat, short

sc.settings.verbosity = 0
sc.set_figure_params(dpi=80, frameon=False)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 300)

GROUP_OF = {g: k for k, gs in KETONE.items() for g in gs}
KG = KETONE_GENES
CORE = ["Hmgcs2", "Hmgcl", "Bdh1", "Bdh2", "Oxct1", "Acat1", "Acss1", "Acss2", "Slc16a1", "Slc16a7", "Slc16a3", "Slc5a8", "Hcar2", "Ffar3", "Cpt1a", "Acadm", "Ppara", "Fgf21"]
CTS_TEST = [c for c in om.COARSE_TYPES if c != "Other"]

PATHS_INHOUSE = {
    "Xenium EAE":    om.resolve_path("OLIGOMETAB_XENIUM_EAE_H5AD",   "../../data/RREAE_5k_raw_only_integration_processed.h5ad"),
    "Visium aging":  om.resolve_path("OLIGOMETAB_VISIUM_AGING_H5AD", "../../data/visum_aging_brain.h5ad"),
    "Falcão EAE sc": om.resolve_path("OLIGOMETAB_FALCAO_H5AD",       "../../data/falcao_et_al_2018.h5ad"),
}
PATHS_PUBLIC = om.processed_paths()
only = os.getenv("OLIGOMETAB_ONLY_DATASETS")
if only:
    PATHS_PUBLIC = {k: v for k, v in PATHS_PUBLIC.items() if k in only.split(",")}
for k, v in PATHS_INHOUSE.items():
    print(f"{k:13s} -> {v if v else 'NOT FOUND (section will be skipped)'}")
print(len(PATHS_PUBLIC), "public datasets")
display(pd.DataFrame([(k, ", ".join(v)) for k, v in KETONE.items()], columns=["group", "genes"]).set_index("group"))
os.makedirs("../../results", exist_ok=True)
results = {}


def codetect(adata, groupby, a="Oxct1", b="Bdh1", min_cells=30):
    """Fraction of cells per group detecting a, b, and both (complete ketolytic machinery when a=Oxct1, b=Bdh1)."""
    df = om.expr_df(adata, [a, b])
    if a not in df.columns or b not in df.columns:
        return pd.DataFrame()
    g = adata.obs[groupby].astype(str).values
    rows = []
    for grp, sub in df.groupby(g):
        if len(sub) < min_cells:
            continue
        pa, pb = (sub[a] > 0).mean(), (sub[b] > 0).mean()
        both = ((sub[a] > 0) & (sub[b] > 0)).mean()
        rows.append({"group": grp, "n_cells": len(sub), f"{a}+": pa, f"{b}+": pb, "both": both, "both_expected_if_independent": pa * pb})
    return pd.DataFrame(rows).set_index("group")


def astro_coupling(adata, mask, genes=("Cpt1a", "Acadm", "Hadha", "Ppara", "Bdh1", "Hmgcl", "Slc16a3", "Gfap")):
    """Inside astrocytes: Spearman correlation of Hmgcs2 with FA-oxidation / regulator genes and detection in Hmgcs2+ vs Hmgcs2- cells."""
    sub = adata[mask]
    if sub.n_obs < 100 or om.resolve(sub, ["Hmgcs2"])["Hmgcs2"] is None or (om.expr(sub, "Hmgcs2") > 0).sum() < 20:
        return pd.DataFrame()
    return om.coexpression_with_anchor(sub, "Hmgcs2", list(genes), min_pos=20)


def gene_heat(df, title, cbar, fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, annot=None, row_labels=None):
    cols = [g for g in CORE if g in df.columns]
    heat(df[cols], title, cbar, fmt=fmt, cmap=cmap, center=center, vmin=vmin, vmax=vmax, annot=annot, row_labels=row_labels)

# %% [markdown]
# ---
# ## 1. Xenium mouse EAE (5K panel)
#
# Panel content limits this section: of the 26 ketone genes the 5K panel carries *Bdh1*, *Acss1*, *Acss2*, *Slc16a1*, *Slc16a3*, *Hcar2*, *Ffar2*, *Ffar3*,
# *Cpt1a*, *Cpt2*, *Hadha*, *Acaa2*, *Ppara*, *Ppargc1a*. No *Hmgcs2*, *Hmgcl*, *Oxct1*, *Acat1*, *Slc16a7*. The matrix is raw counts and is normalised here.

# %%
ds = "Xenium EAE"
ad_eae = None
if PATHS_INHOUSE[ds]:
    from anndata import AnnData
    from anndata.experimental import read_elem
    with h5py.File(PATHS_INHOUSE[ds], "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        shape = tuple(f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X/data"][:].astype(np.float32), f["X/indices"][:], f["X/indptr"][:]), shape=shape)
        spatial = f["obsm/spatial"][:]
    obs.index = obs.index.astype(str); var.index = var.index.astype(str)
    ad_eae = AnnData(X=X, obs=obs, var=var); ad_eae.obsm["spatial"] = spatial; ad_eae.obs_names_make_unique(); del X, obs, var
    ad_eae.obs["cell_type"] = ad_eae.obs["cell_type"].astype(str).replace({"DA oligodendrocytes": "DA Oligodendrocytes", "Astrocyte": "Astrocytes"})
    sc.pp.normalize_total(ad_eae, target_sum=1e4); sc.pp.log1p(ad_eae)
    kg_eae = om.present(ad_eae, KG)
    print(ad_eae.shape, "| ketone genes present:", kg_eae, "| missing:", om.missing(ad_eae, KG))
    results[ds] = {}
else:
    print("skipped")

# %%
if ad_eae is not None:
    vc = ad_eae.obs["cell_type"].value_counts(); cts = [c for c in vc.index if vc[c] >= 1500 and c != "unclear"]
    det = om.detection_table(ad_eae, kg_eae, "cell_type", min_cells=1500); results[ds]["detection_by_celltype"] = det
    heat(det[[g for g in kg_eae]], "Xenium EAE: fraction of cells with ≥1 transcript, ketone-related genes on the 5K panel", "fraction detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, row_labels=list(det.index))
    mt = om.mean_table(ad_eae, kg_eae, "cell_type", min_cells=1500); results[ds]["mean_by_celltype"] = mt
    # EAE vs control per cell type, and DA vs homeostatic oligodendrocytes / DA vs homeostatic astrocytes
    sub = ad_eae[ad_eae.obs["cell_type"].isin(cts)]
    n = sub.obs.groupby(["sample_name", "cell_type"], observed=True).size().rename("n_cells").reset_index(); n[["sample_name", "cell_type"]] = n[["sample_name", "cell_type"]].astype(str)
    pb = om.pseudobulk(sub, kg_eae, "sample_name", ["cell_type", "condition"]).merge(n, on=["sample_name", "cell_type"]).query("n_cells >= 20")
    rows = []
    for ct in cts:
        p_ct = pb[pb.cell_type == ct]
        for g in kg_eae:
            ma, mr, p, na, nr = om.group_test(p_ct, g, "condition", "EAE", "CONTROL")
            rows.append({"cell_type": ct, "gene": g, "mean_EAE": ma, "mean_CONTROL": mr, "log2FC": om.log2fc(ma, mr), "MWU_p": p, "n_EAE": na, "n_CONTROL": nr})
    de = pd.DataFrame(rows); results[ds]["EAE_vs_control"] = de
    d = de.pivot(index="cell_type", columns="gene", values="log2FC").loc[cts, kg_eae]; d = d[d.notna().any(axis=1)]
    p = de.pivot(index="cell_type", columns="gene", values="MWU_p").loc[d.index, kg_eae]
    heat(d, "Xenium EAE: pseudobulk log2FC EAE vs CONTROL, ketone-related genes (* MWU p < 0.05)", "log2 fold change", fmt=".1f", cmap="RdBu_r", center=0, vmin=-2, vmax=2,
         annot=(d.round(1).astype(str) + np.where(p.values < 0.05, "*", "")).values, row_labels=list(d.index))
    # paired within sample: DA vs homeostatic oligodendrocytes, DA vs homeostatic astrocytes
    def paired_genes(mask_a, mask_b, la, lb, min_cells=20):
        df = om.expr_df(ad_eae, kg_eae); smp = ad_eae.obs["sample_name"].astype(str).values
        A, B = [], []
        for s_ in np.unique(smp):
            m = smp == s_
            if (m & mask_a).sum() >= min_cells and (m & mask_b).sum() >= min_cells:
                A.append(df[m & mask_a].mean()); B.append(df[m & mask_b].mean())
        A, B = pd.DataFrame(A), pd.DataFrame(B)
        return pd.DataFrame([{"gene": g, f"mean_{la}": A[g].mean(), f"mean_{lb}": B[g].mean(), "log2FC": om.log2fc(A[g].mean(), B[g].mean()), "n_samples": len(A),
                              "wilcoxon_p": stats.wilcoxon(A[g], B[g]).pvalue if len(A) >= 3 and not np.allclose(A[g], B[g]) else np.nan} for g in df.columns]).set_index("gene")
    ct_ = ad_eae.obs["cell_type"].values
    t1 = paired_genes(ct_ == "DA Oligodendrocytes", ct_ == "Oligodendrocytes", "DA_oligo", "homeostatic_oligo"); results[ds]["DA_vs_homeostatic_oligo"] = t1
    t2 = paired_genes(ct_ == "DA astrocytes", ct_ == "Astrocytes", "DA_astro", "astro"); results[ds]["DA_vs_homeostatic_astro"] = t2
    print("DA vs homeostatic oligodendrocytes (paired per sample):"); display(t1.round(4))
    print("DA vs homeostatic astrocytes (paired per sample):"); display(t2.round(4))
    # lesion distance
    m_eae = ((ad_eae.obs["condition"] == "EAE") & ad_eae.obs["lesion_distance_bin"].notna()).values
    import re
    bins = sorted(ad_eae.obs.loc[m_eae, "lesion_distance_bin"].astype(str).unique(), key=lambda b: (1, 0) if b.startswith(">") else (0, int(re.search(r"(\d+)", b).group(1))))
    ld_rows = []
    for label, cts_ in [("oligodendrocyte lineage", ["Oligodendrocytes", "DA Oligodendrocytes"]), ("astrocytes", ["Astrocytes", "DA astrocytes"]),
                        ("myeloid", ["Microglia", "Macrophages", "Foamy Mic_Mac", "Activate Mic_Mac 1", "Activated Mic_Mac 2", "Efflux Mic_Mac", "Myeloid cells", "Proliferating microglia"]), ("neurons", ["Neurons"])]:
        s_ = ad_eae[m_eae & ad_eae.obs["cell_type"].isin(cts_).values]
        mt_ = om.mean_table(s_, [g for g in ["Bdh1", "Slc16a1", "Slc16a3", "Hcar2", "Cpt1a", "Ppara", "Acss2"] if g in kg_eae], "lesion_distance_bin", min_cells=50).reindex(bins)
        mt_["population"] = label; ld_rows.append(mt_.reset_index().rename(columns={"index": "bin"})); del s_
    ld = pd.concat(ld_rows); results[ds]["lesion_distance"] = ld
    tidy = ld.melt(id_vars=["bin", "population"], var_name="gene", value_name="mean_expr")
    g = sns.relplot(data=tidy, x="bin", y="mean_expr", hue="population", col="gene", col_wrap=4, kind="line", marker="o", height=2.4, aspect=1.2, facet_kws={"sharey": False})
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=45)
    g.set_titles("{col_name}"); plt.show()
    # C4b co-expression with the ketone genes in DA oligodendrocytes
    dao = ad_eae[ad_eae.obs["cell_type"] == "DA Oligodendrocytes"].copy()
    co = om.coexpression_with_anchor(dao, "C4b", kg_eae); results[ds]["C4b_coexpr_DAoligo"] = co
    print("ketone genes vs C4b in DA oligodendrocytes:"); display(co[["frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
    del dao, sub
    for k, v in results[ds].items():
        v.to_csv(f"../../results/ketone_xenium_eae_{k}.csv")
    del ad_eae; gc.collect()

# %% [markdown]
# ---
# ## 2. Visium aging mouse brain (whole transcriptome)

# %%
ds = "Visium aging"
ad_vis = None
if PATHS_INHOUSE[ds]:
    ad_vis = sc.read_h5ad(PATHS_INHOUSE[ds]); ad_vis.var_names_make_unique()
    ad_vis.obs["age_group"] = pd.Categorical(ad_vis.obs["age_group"].astype(str), categories=["Young", "Mid", "Old"])
    ad_vis.obs["age_months"] = ad_vis.obs["age_group"].map({"Young": 6, "Mid": 18, "Old": 21}).astype(float)
    wm_genes = om.present(ad_vis, ["Plp1", "Mbp", "Mobp", "Mag", "Cldn11"]); sc.tl.score_genes(ad_vis, wm_genes, score_name="wm_score", use_raw=False)
    ad_vis.obs["wm_rich"] = ad_vis.obs["wm_score"] >= ad_vis.obs["wm_score"].quantile(0.75)
    kg_vis = om.present(ad_vis, KG); print(ad_vis.shape, "| missing:", om.missing(ad_vis, KG))
    results[ds] = {}
    P_vis = om.pathway_scores(ad_vis, KETONE_SCORE_SETS)
    # detection / mean by age group and WM vs non-WM
    ad_vis.obs["compartment_age"] = np.where(ad_vis.obs["wm_rich"], "WM-rich ", "other ") + ad_vis.obs["age_group"].astype(str)
    det = om.detection_table(ad_vis, kg_vis, "compartment_age"); results[ds]["detection_by_compartment_age"] = det
    gene_heat(det, "Visium aging: fraction of spots with ≥1 UMI, by compartment and age", "fraction of spots", row_labels=list(det.index))
    mt = om.mean_table(ad_vis, kg_vis, "compartment_age"); results[ds]["mean_by_compartment_age"] = mt
    gene_heat(mt, "Visium aging: mean log-expression by compartment and age", "mean log-expression", vmin=None, row_labels=list(mt.index))
    # age slopes per section, all spots and WM-rich
    rows = []
    for label, m in [("all spots", np.ones(ad_vis.n_obs, bool)), ("WM-rich", ad_vis.obs["wm_rich"].values)]:
        a_ = ad_vis[m]
        pb = om.pseudobulk(a_, kg_vis, "sample", ["age_group"]); pb["age_months"] = pb["age_group"].map({"Young": 6, "Mid": 18, "Old": 21}).astype(float)
        pbs = om.pseudobulk(a_, [], "sample", ["age_group"], values=P_vis.loc[a_.obs_names]); pbs["age_months"] = pb["age_months"].values
        for g in kg_vis:
            r = stats.linregress(pb["age_months"], pb[g])
            rows.append({"spots": label, "feature": g, "kind": "gene", "pearson_r": r.rvalue, "p": r.pvalue, "mean_Young": pb.loc[pb.age_group == "Young", g].mean(), "mean_Old": pb.loc[pb.age_group == "Old", g].mean(),
                         "log2FC_old_vs_young": om.log2fc(pb.loc[pb.age_group == "Old", g].mean(), pb.loc[pb.age_group == "Young", g].mean())})
        for pw in P_vis.columns:
            r = stats.linregress(pbs["age_months"], pbs[pw])
            rows.append({"spots": label, "feature": pw, "kind": "score", "pearson_r": r.rvalue, "p": r.pvalue, "mean_Young": pbs.loc[pbs.age_group == "Young", pw].mean(), "mean_Old": pbs.loc[pbs.age_group == "Old", pw].mean(), "log2FC_old_vs_young": np.nan})
        del a_
    age = pd.DataFrame(rows); results[ds]["age_trends"] = age
    print("age trends (linear fit over 6 sections):"); display(age.round(4))
    # C4b co-expression within WM-rich spots
    wm = ad_vis[ad_vis.obs["wm_rich"].values].copy()
    co = om.coexpression_with_anchor(wm, "C4b", kg_vis); results[ds]["C4b_coexpr_WM"] = co
    print("ketone genes vs C4b in white-matter-rich spots:"); display(co[["frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
    del wm
    # maps: one section per age group
    picks = ad_vis.obs.groupby("age_group", observed=True)["sample"].first()
    for age_, s_ in picks.items():
        sub = ad_vis[ad_vis.obs["sample"] == s_].copy()
        hv = om.present(sub, ["Hmgcs2", "Bdh1", "Oxct1", "Slc16a7", "Slc16a1", "Hcar2"])
        sc.pl.spatial(sub, color=[om.resolve(sub, [g])[g] for g in hv], spot_size=90, cmap="magma", vmax="p99", ncols=len(hv), show=False); plt.suptitle(f"{s_} ({age_})", y=1.02); plt.show()
        del sub
    for k, v in results[ds].items():
        v.to_csv(f"../../results/ketone_visium_{k}.csv")
    del ad_vis, P_vis; gc.collect()

# %% [markdown]
# ---
# ## 3. Falcão et al. 2018 (sorted single cells): intrinsic ketone handling along the oligodendrocyte lineage

# %%
ds = "Falcão EAE sc"
ad_f = None
if PATHS_INHOUSE[ds]:
    ad_f = sc.read_h5ad(PATHS_INHOUSE[ds]); ad_f.var_names_make_unique()
    cl = ad_f.obs["Renamed_clusternames"].astype(str)
    grp = np.select([cl.str.contains("MOL") & cl.str.contains("EAE"), cl.str.contains("MOL"), cl.str.startswith("MiGl"), cl.str.startswith("VLMC"),
                     cl.isin(["OPC1", "OPC2", "OPC3", "OPC_Cycling", "COP", "NFOL", "PLC"])], ["MOL (EAE clusters)", "MOL (control clusters)", "Microglia", "VLMC", "OPC/COP/NFOL"], default="other")
    ad_f.obs["group_coarse"] = pd.Categorical(grp, categories=["OPC/COP/NFOL", "MOL (control clusters)", "MOL (EAE clusters)", "Microglia", "VLMC", "other"])
    kg_f = om.present(ad_f, KG); print(ad_f.shape, "| missing:", om.missing(ad_f, KG))
    results[ds] = {}
    r = om.resolve(ad_f, kg_f)
    sc.pl.dotplot(ad_f, var_names=[r[g] for g in kg_f], groupby="Renamed_clusternames", standard_scale="var", color_map="Reds", figsize=(0.4 * len(kg_f) + 2, 6), show=False, title="Falcão: ketone-related genes by cluster"); plt.show()
    det = om.detection_table(ad_f, kg_f, "group_coarse", min_cells=10); results[ds]["detection_by_group"] = det
    gene_heat(det, "Falcão sorted cells: fraction detecting each ketone-related gene", "fraction of cells", row_labels=list(det.index))
    cd = codetect(ad_f, "group_coarse", min_cells=10); results[ds]["codetect_Oxct1_Bdh1"] = cd
    print("complete ketolytic machinery (Oxct1 and Bdh1 in the same cell):"); display(cd.round(3))
    cd_cl = codetect(ad_f, "Renamed_clusternames", min_cells=8); results[ds]["codetect_by_cluster"] = cd_cl; display(cd_cl.round(3).sort_values("both", ascending=False))
    # EAE-MOL vs control-MOL, gene-wise
    m_eae = (ad_f.obs["group_coarse"] == "MOL (EAE clusters)").values; m_ctl = (ad_f.obs["group_coarse"] == "MOL (control clusters)").values
    df = om.expr_df(ad_f, kg_f); rows = []
    for g in df.columns:
        x = df[g].values
        rows.append({"gene": g, "group": GROUP_OF.get(g), "frac_MOL_EAE": (x[m_eae] > 0).mean(), "frac_MOL_ctrl": (x[m_ctl] > 0).mean(), "mean_MOL_EAE": x[m_eae].mean(), "mean_MOL_ctrl": x[m_ctl].mean(),
                     "log2FC": om.log2fc(x[m_eae].mean(), x[m_ctl].mean()), "MWU_p_cells": stats.mannwhitneyu(x[m_eae], x[m_ctl]).pvalue if x.std() > 0 else np.nan})
    t = pd.DataFrame(rows).set_index("gene"); results[ds]["MOL_EAE_vs_ctrl"] = t
    print("EAE-MOL vs control-MOL (cells as replicates):"); display(t.round(4))
    # C4b co-expression in MOL
    mol = ad_f[ad_f.obs["group_coarse"].isin(["MOL (control clusters)", "MOL (EAE clusters)"])].copy()
    co = om.coexpression_with_anchor(mol, "C4b", kg_f); results[ds]["C4b_coexpr_MOL"] = co
    print("ketone genes vs C4b in MOL:"); display(co[["frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
    del mol
    for k, v in results[ds].items():
        v.to_csv(f"../../results/ketone_falcao_{k}.csv")
    del ad_f; gc.collect()

# %% [markdown]
# ---
# ## 4. Public datasets: per-dataset pass
#
# For each dataset: detection and mean expression of every ketone gene by coarse cell type; *Oxct1* + *Bdh1* co-detection; ketogenesis / ketolysis /
# fatty-acid-supply scores by cell type; pseudobulk contrasts for every cell type, gene-wise and score-wise; *Hmgcs2* coupling inside astrocytes; *C4b*
# co-expression inside oligodendrocytes (mouse); by-condition tables (human). Cross-dataset summaries follow in section 5.

# %%
DET, MEAN, CODET, SCORE_CT, DE, DE_S, ASTRO, CO, COND, INFO, ERR = [], [], [], [], [], [], {}, {}, {}, {}, {}
for name, path in PATHS_PUBLIC.items():
    try:
        a = sc.read_h5ad(path); a.obs["cell_type_coarse"] = a.obs["cell_type_coarse"].astype(str)
        species, modality = str(a.obs["species"].iloc[0]), str(a.obs["modality"].iloc[0]); ref, alts = om.alt_groups(a)
        INFO[name] = {"species": species, "modality": modality, "n_cells": a.n_obs, "n_samples": a.obs["sample"].nunique()}
        kg = om.present(a, KG)
        print(f"\n{'=' * 110}\n{name}: {a.n_obs} cells | {species} {modality} | missing ketone genes: {om.missing(a, KG)}")
        P = om.pathway_scores(a, KETONE_SCORE_SETS)
        key_ct = "cell_type_coarse"
        if modality == "spatial":
            key_ct = "group_age" if "age" in a.obs else ("condition_original" if "condition_original" in a.obs else "group")
            if key_ct == "group_age":
                a.obs["group_age"] = a.obs["group"].astype(str) + "_" + a.obs["age"].astype(str)
        det = om.detection_table(a, kg, key_ct, min_cells=30); det.insert(0, "dataset", name); DET.append(det.reset_index().rename(columns={"group": "cell_type"}))
        mt = om.mean_table(a, kg, key_ct, min_cells=30); mt.insert(0, "dataset", name); MEAN.append(mt.reset_index().rename(columns={"index": "cell_type"}))
        cd = codetect(a, key_ct)
        if len(cd):
            cd.insert(0, "dataset", name); CODET.append(cd.reset_index().rename(columns={"group": "cell_type"}))
        s_ct = P.groupby(a.obs[key_ct].astype(str).values).mean(); n_ct = a.obs[key_ct].astype(str).value_counts(); s_ct = s_ct.loc[n_ct[n_ct >= 30].index.intersection(s_ct.index)]
        s_ct.insert(0, "dataset", name); SCORE_CT.append(s_ct.reset_index().rename(columns={"index": "cell_type"}))
        gene_heat(det.drop(columns="dataset"), f"{short(name)}: fraction of cells with ≥1 UMI", "fraction detected", row_labels=list(det.index))
        # contrasts per cell type (or per spot)
        groups_ct = [c for c in CTS_TEST if c in a.obs[key_ct].unique()] if modality != "spatial" else ["spot"]
        for ct in groups_ct:
            m_ct = (a.obs[key_ct] == ct).values if modality != "spatial" else np.ones(a.n_obs, bool)
            if m_ct.sum() < 50 or not alts:
                continue
            sub = a[m_ct]; n = sub.obs.groupby("sample", observed=True).size(); ok = n[n >= 20].index
            pb = om.pseudobulk(sub, kg); pb = pb[pb["sample"].isin(ok)]
            pbs = om.pseudobulk(sub, [], values=P.loc[sub.obs_names]); pbs = pbs[pbs["sample"].isin(ok)]
            for alt in alts:
                for g in [c for c in pb.columns if c in kg]:
                    ma, mr, p, na, nr = om.group_test(pb, g, "group", alt, ref)
                    DE.append({"dataset": name, "cell_type": ct, "gene": g, "group_of": GROUP_OF.get(g), "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr, "log2FC": om.log2fc(ma, mr), "MWU_p": p, "n_alt": na, "n_ref": nr})
                for pw in P.columns:
                    ma, mr, p, na, nr = om.group_test(pbs, pw, "group", alt, ref)
                    DE_S.append({"dataset": name, "cell_type": ct, "score": pw, "alt": alt, "ref": ref, "mean_alt": ma, "mean_ref": mr, "diff": ma - mr, "MWU_p": p, "n_alt": na, "n_ref": nr})
            del sub
        if modality != "spatial":
            ac = astro_coupling(a, (a.obs["cell_type_coarse"] == "Astrocyte").values)
            if len(ac):
                ASTRO[name] = ac; print("astrocytes: genes vs Hmgcs2"); display(ac[["frac_in_Hmgcs2-", "frac_in_Hmgcs2+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
            ol = a[(a.obs["cell_type_coarse"] == "Oligodendrocyte").values].copy()
            if species == "mouse" and ol.n_obs >= 100 and om.resolve(ol, ["C4b"])["C4b"] is not None and (om.expr(ol, "C4b") > 0).sum() >= 50:
                CO[name] = om.coexpression_with_anchor(ol, "C4b", kg); print("oligodendrocytes: ketone genes vs C4b"); display(CO[name][["frac_in_C4b-", "frac_in_C4b+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
            del ol
            key = "condition_original" if ("condition_original" in a.obs and a.obs["condition_original"].nunique() > 1) else ("braak" if "braak" in a.obs else None)
            if key:
                for ct in ["Oligodendrocyte", "Astrocyte", "Microglia", "OPC", "Neuron"]:
                    m_ct = (a.obs["cell_type_coarse"] == ct).values
                    if m_ct.sum() >= 50:
                        sub = a[m_ct]; t = om.detection_table(sub, kg, key, min_cells=30); t2 = om.mean_table(sub, kg, key, min_cells=30)
                        COND[(name, ct, "detection")] = t; COND[(name, ct, "mean")] = t2; del sub
        else:
            print(f"mean log-expression by {key_ct}:"); display(mt.drop(columns="dataset")[[g for g in CORE if g in mt.columns]].round(3))
    except Exception as exc:
        ERR[name] = traceback.format_exc(); print(f"!! {name} failed: {type(exc).__name__}: {exc}")
    finally:
        for v in ["a", "P", "sub", "ol"]:
            if v in globals():
                del globals()[v]
        gc.collect()
DET, MEAN, SCORE_CT, DE, DE_S = pd.concat(DET, ignore_index=True), pd.concat(MEAN, ignore_index=True), pd.concat(SCORE_CT, ignore_index=True), pd.DataFrame(DE), pd.DataFrame(DE_S)
CODET = pd.concat(CODET, ignore_index=True) if CODET else pd.DataFrame()
print("\ndone:", len(INFO), "datasets;", len(ERR), "failed", list(ERR))
for k, v in ERR.items():
    print(k, v)

# %% [markdown]
# ## 5. Who makes, who uses, who transports and who senses ketone bodies?
#
# Detection of each ketone gene by coarse cell type, averaged over the mouse whole-cell datasets and over the mouse nucleus datasets separately (nuclei detect
# cytoplasmic / mitochondrial transcripts less), and over the human nucleus datasets.

# %%
def avg_detection(datasets, title):
    sub = DET[DET.dataset.isin(datasets) & DET.cell_type.isin(CTS_TEST)]
    if not len(sub):
        return None
    m = sub.groupby("cell_type")[[g for g in CORE if g in sub.columns]].mean().reindex([c for c in CTS_TEST if c in sub.cell_type.unique()])
    heat(m, title, "mean fraction of cells detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, row_labels=list(m.index))
    return m
mouse_sc = [k for k, v in INFO.items() if v["species"] == "mouse" and v["modality"] == "scRNA"]
mouse_sn = [k for k, v in INFO.items() if v["species"] == "mouse" and v["modality"] == "snRNA"]
human_sn = [k for k, v in INFO.items() if v["species"] == "human" and v["modality"] == "snRNA"]
WHO = {}
WHO["mouse whole cells"] = avg_detection(mouse_sc, f"Mouse whole-cell datasets ({len(mouse_sc)}): fraction of cells detecting each ketone gene, by cell type")
WHO["mouse nuclei"] = avg_detection(mouse_sn, f"Mouse nucleus datasets ({len(mouse_sn)}): fraction of nuclei detecting each ketone gene, by cell type")
WHO["human nuclei"] = avg_detection(human_sn, f"Human nucleus datasets ({len(human_sn)}): fraction of nuclei detecting each ketone gene, by cell type")
if len(CODET):
    c = CODET[CODET.cell_type.isin(CTS_TEST)].groupby("cell_type")[["Oxct1+", "Bdh1+", "both", "both_expected_if_independent"]].mean().reindex([c for c in CTS_TEST if c in CODET.cell_type.unique()])
    print("complete ketolytic machinery (Oxct1 and Bdh1 in the same cell), mean over datasets:"); display(c.round(3))
    per = CODET.pivot_table(index="dataset", columns="cell_type", values="both")[[c for c in CTS_TEST if c in CODET.cell_type.unique()]]
    heat(per, "Fraction of cells co-detecting Oxct1 and Bdh1, per dataset and cell type", "fraction with both", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None)
s = SCORE_CT[SCORE_CT.cell_type.isin(CTS_TEST)]
z = []
for dsn in s.dataset.unique():
    b = s[s.dataset == dsn].set_index("cell_type").drop(columns="dataset"); z.append((b - b.mean()) / b.std(ddof=0))
z = pd.concat(z).groupby(level=0).mean().reindex([c for c in CTS_TEST if c in s.cell_type.unique()])
heat(z, "Ketogenesis / ketolysis / FA-supply scores by cell type (z-scored within dataset, averaged over datasets)", "z-score across cell types", fmt=".1f", cmap="RdBu_r", center=0, row_labels=list(z.index))

# %% [markdown]
# ## 6. Disease, age and demyelination effects on the ketone genes, per cell type

# %%
DE["contrast"] = DE["dataset"].map(short) + "  [" + DE["alt"] + " vs " + DE["ref"] + "]"
DE_S["contrast"] = DE_S["dataset"].map(short) + "  [" + DE_S["alt"] + " vs " + DE_S["ref"] + "]"
def pv(df, index, columns, values):
    return df.drop_duplicates([index, columns]).set_index([index, columns])[values].unstack(columns)
for ct in [c for c in CTS_TEST + ["spot"] if c in DE.cell_type.unique()]:
    sub = DE[DE.cell_type == ct]
    d = pv(sub, "contrast", "gene", "log2FC"); d = d[[g for g in CORE if g in d.columns]]
    p = pv(sub, "contrast", "gene", "MWU_p").reindex(index=d.index, columns=d.columns)
    heat(d, f"{ct}: pseudobulk log2 fold change of ketone-related genes (* Mann–Whitney p < 0.05)", "log2 fold change", fmt=".1f", cmap="RdBu_r", center=0, vmin=-2, vmax=2,
         annot=(d.round(1).astype(str).replace("nan", "") + np.where(p.values < 0.05, "*", "")).values, row_labels=list(d.index))
for ct in ["Oligodendrocyte", "Astrocyte", "Microglia", "OPC", "Neuron"]:
    sub = DE_S[DE_S.cell_type == ct]
    if not len(sub):
        continue
    d = pv(sub, "contrast", "score", "diff"); p = pv(sub, "contrast", "score", "MWU_p").reindex(index=d.index, columns=d.columns)
    heat(d, f"{ct}: ketogenesis / ketolysis / FA-supply score difference (* MWU p < 0.05)", "score difference", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.2, vmax=0.2,
         annot=(d.round(2).astype(str).replace("nan", "") + np.where(p.values < 0.05, "*", "")).values, row_labels=list(d.index))
# consistency per cell type
cons_rows = []
for ct in [c for c in CTS_TEST if c in DE.cell_type.unique()]:
    sub = DE[DE.cell_type == ct]
    for g, gg in sub.groupby("gene"):
        cons_rows.append({"cell_type": ct, "gene": g, "n_contrasts": len(gg), "n_up": int((gg.log2FC > 0.5).sum()), "n_down": int((gg.log2FC < -0.5).sum()),
                          "median_log2FC": gg.log2FC.median(), "n_sig_up": int(((gg.MWU_p < 0.05) & (gg.log2FC > 0)).sum()), "n_sig_down": int(((gg.MWU_p < 0.05) & (gg.log2FC < 0)).sum())})
CONS = pd.DataFrame(cons_rows); CONS["net"] = CONS.n_up - CONS.n_down
print("consistency of direction across contrasts (|log2FC| > 0.5), all cell types:")
display(CONS.pivot(index="gene", columns="cell_type", values="net").reindex([g for g in KG if g in CONS.gene.unique()]).fillna(0).astype(int))
print("median log2FC across contrasts:"); display(CONS.pivot(index="gene", columns="cell_type", values="median_log2FC").reindex([g for g in KG if g in CONS.gene.unique()]).round(2))

# %% [markdown]
# ## 7. Astrocyte ketogenesis: is *Hmgcs2* coupled to fatty-acid oxidation and PPARα?

# %%
if ASTRO:
    R = pd.DataFrame({k: v["spearman_rho"] for k, v in ASTRO.items()}).T
    heat(R, "Inside astrocytes: Spearman correlation of each gene with Hmgcs2, per dataset", "Spearman rho with Hmgcs2", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.2, vmax=0.2)
    E = pd.DataFrame({k: v["enrichment(+/-)"] for k, v in ASTRO.items()}).T
    print("detection enrichment in Hmgcs2+ vs Hmgcs2- astrocytes:"); display(E.round(2))
    frac = DET[DET.cell_type == "Astrocyte"].set_index("dataset")[["n_cells"] + [g for g in ["Hmgcs2", "Hmgcl", "Bdh1", "Oxct1", "Cpt1a", "Acadm", "Ppara", "Slc16a3", "Fgf21"] if g in DET.columns]]
    print("astrocytes: fraction detecting the ketogenic / FA-oxidation genes per dataset:"); display(frac.round(3))

# %% [markdown]
# ## 8. Ketone genes and the *C4b*⁺ oligodendrocyte state (mouse datasets)

# %%
if CO:
    R = pd.DataFrame({k: v["spearman_rho"] for k, v in CO.items()}).T; R = R[[g for g in CORE if g in R.columns]]
    heat(R, "Spearman correlation of ketone genes with C4b inside oligodendrocytes, per mouse dataset", "Spearman rho with C4b", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.15, vmax=0.15)
    E = pd.DataFrame({k: v["enrichment(+/-)"] for k, v in CO.items()}).T; E = E[[g for g in CORE if g in E.columns]]
    heat(np.log2(E.replace(0, np.nan)), "Detection enrichment in C4b+ vs C4b- oligodendrocytes (log2), per mouse dataset", "log2 detection ratio", fmt=".1f", cmap="RdBu_r", center=0, vmin=-2, vmax=2)

# %% [markdown]
# ## 9. Human MS and AD by the authors' condition labels

# %%
for (name, ct, kind), t in COND.items():
    if kind != "detection" or ct not in ("Oligodendrocyte", "Astrocyte", "Microglia"):
        continue
    gene_heat(t, f"{short(name)} — {ct}: fraction detected by condition", "fraction detected", row_labels=list(t.index))

# %% [markdown]
# ## 10. Summary tables

# %%
display(pd.DataFrame(INFO).T)
DET.to_csv("../../results/ketone_public_detection_by_celltype.csv", index=False); MEAN.to_csv("../../results/ketone_public_mean_by_celltype.csv", index=False)
if len(CODET): CODET.to_csv("../../results/ketone_public_codetection_Oxct1_Bdh1.csv", index=False)
SCORE_CT.to_csv("../../results/ketone_public_scores_by_celltype.csv", index=False)
DE.to_csv("../../results/ketone_public_gene_contrasts.csv", index=False); DE_S.to_csv("../../results/ketone_public_score_contrasts.csv", index=False)
CONS.to_csv("../../results/ketone_public_gene_consistency.csv", index=False)
if ASTRO: pd.concat(ASTRO, names=["dataset", "gene"]).to_csv("../../results/ketone_public_astrocyte_Hmgcs2_coupling.csv")
if CO: pd.concat(CO, names=["dataset", "gene"]).to_csv("../../results/ketone_public_c4b_coexpression_oligodendrocytes.csv")
if COND: pd.concat({f"{k[0]}|{k[1]}|{k[2]}": v for k, v in COND.items()}, names=["dataset|cell_type|kind", "condition"]).to_csv("../../results/ketone_public_by_condition.csv")
for k, v in WHO.items():
    if v is not None:
        v.to_csv(f"../../results/ketone_who_{k.replace(' ', '_')}.csv")
print("results written")

# %% [markdown]
# ### Interpretation
#
# **Who makes ketone bodies.** Almost nobody in brain parenchyma at the transcript level. *Hmgcs2*, the rate-limiting ketogenic enzyme, is detected in 3 % of astrocytes in whole cells and nuclei alike, in 0 % of oligodendrocytes and microglia, and in 0 % of any human nucleus type; where it does appear it is vascular: 12–22 % of endothelial cells, 7–14 % of mural / fibroblast cells, 53 % of the sorted VLMC in Falcão, and a subset of neurons in the whole-cell datasets (15 %). *Hmgcl* is broad and *Bdh1* (which works in both directions) sits in astrocytes, OPCs, ependymal cells and neurons (16–27 %), not in oligodendrocytes (5–7 %) or microglia (1–2 %). Astrocyte ketogenesis therefore rests on a small *Hmgcs2*⁺ subset; inside astrocytes *Hmgcs2*⁺ cells are 1.5–3.8× more likely to express *Acadm*, *Cpt1a* and *Hmgcl* in the demyelination and aged-white-matter datasets, i.e. the ketogenic astrocytes are the fatty-acid-oxidising ones, but the coupling to *Ppara* is weak (ρ ≤ 0.05 everywhere).
#
# **Who uses them.** *Oxct1* (SCOT), the committed ketolytic step, is the most widely expressed ketone gene: 65 % of neuronal, 48 % of astrocytic, 33–39 % of oligodendrocyte and 43 % of microglial nuclei / cells, 82 % of OPC / COP / NFOL in Falcão. The *complete* machinery (*Oxct1* and *Bdh1* in the same cell) is rarer and follows *Bdh1*: ependymal cells 18 %, neurons 13 %, OPCs 7 %, astrocytes 6 %, oligodendrocytes 2 %, microglia < 1 %; in sorted cells 45 % of OPC / COP / NFOL, 9 % of control MOL, 23 % of EAE MOL, 3 % of microglia. Co-detection is at the level expected from the two marginals, so there is no dedicated "ketolytic" subpopulation beyond what *Bdh1* defines. Along the oligodendrocyte lineage, ketolytic capacity is a progenitor feature that mature oligodendrocytes lose and partly regain in the disease-associated state.
#
# **Who transports and who senses.** MCT1 (*Slc16a1*) is endothelial (42 % of whole cells) and oligodendroglial (34 %), MCT2 (*Slc16a7*) neuronal (50 % of nuclei), MCT4 (*Slc16a3*) microglial and, in EAE, astrocytic (84 % of DA astrocytes in situ), SMCT1 (*Slc5a8*) essentially absent. *Hcar2* is a myeloid receptor everywhere (6–9 % of microglia / immune cells in whole cells, 8 % of microglia and 21 % of proliferating microglia in EAE in situ, 59 % of sorted microglia in Falcão), *Ffar3* is not detected in brain, *Ffar2* barely.
#
# **Disease.** The ketone genes move modestly and in a reproducible direction. Across the public contrasts, *Bdh1* falls in oligodendrocytes (down in 3 more contrasts than up; median log2FC −0.41; Jäkel MS −1.6, Ximerakis aging −0.6) and in immune cells, *Acss1* / *Acss2* fall in most cell types, *Ppara* falls in oligodendrocytes, microglia and astrocytes, MCT1 falls in oligodendrocytes (−0.34) but rises in astrocytes and immune cells, and *Hcar2* rises in microglia (median +0.8 log2; +1.3 in aged brain, +2.6 in 5XFAD) and immune cells. The ketogenesis and ketolysis *scores* are slightly lower in MS oligodendrocytes, OPCs, microglia and astrocytes in Lerma-Martin 2024 (−0.02 to −0.05, all p < 0.05), in Jäkel oligodendrocytes and in Leng Braak-6 OPCs; nothing rises. In EAE in situ the tissue-wide changes are *Acss2* down (13 cell types), *Ppargc1a* down (10), *Ppara* down (6), MCT1 up (12) and MCT4 up (13), *Hcar2* "up" in 15 cell types (spill-over from myeloid cells, which are the only cells that express it in sorted data), and *Bdh1* up in homeostatic and DA oligodendrocytes (+0.2 log2) but down in microglia, endothelial and fibroblast populations.
#
# **The disease-associated oligodendrocyte.** In sorted EAE MOL versus control MOL, *Bdh1* (2.3×), *Hmgcll1* (2.3×), *Acat1* (1.2×), *Cpt1a* (2.5×), *Cpt2*, *Hadha* and MCT4 rise while MCT1 (−0.24 log2) and *Acss2* fall; *Oxct1* and *Hmgcs2* do not change. In situ the picture is subtler: EAE raises *Bdh1* in both homeostatic and DA oligodendrocytes, but within the same sample DA oligodendrocytes have *less* *Bdh1*, *Acss2*, *Cpt2*, *Ppara* and *Ppargc1a* than homeostatic ones and more MCT4 and *Hcar2* (the latter myeloid spill-over). So the DA state is not a ketogenic or strongly ketolytic state; it is a fatty-acid-importing, lactate-exporting (MCT4 up, MCT1 down) state in which *Bdh1* is the one ketone enzyme that is induced by the disease environment. In white-matter-rich Visium spots *Bdh1* rises with age (r = 0.86, p = 0.03) alongside *Slc16a3*, *Hcar2* and *Bsg*.
#
# **The *C4b* axis.** No ketone gene correlates with *C4b* inside oligodendrocytes beyond |ρ| ≈ 0.1 in any dataset. The detection enrichments in *C4b*⁺ cells that do appear, *Hcar2* (7–43×), *Hmgcs2* (2–7×) and *Cpt1a* (2–6×) in the LPC / cuprizone and *Serpina3n*-cKO nuclei, are the microglial ambient-RNA signature of demyelinating-lesion nuclei documented in OligoC4b, and *Bdh1* is ~2× enriched in *C4b*⁺ sorted MOL and in *C4b*⁺ oligodendrocyte nuclei after LPC / cuprizone. Ketone bodies are therefore not part of the *C4b*⁺ program; what the program shares with the ketone arm is the fatty-acid-oxidation entry (*Cpt1a*, *Hadha*) and the lactate / ketone exporter MCT4.
#
# **Caveats.** *Hmgcs2*, *Hcar2*, *Ffar3*, *Slc5a8* and *Fgf21* are near the detection floor in most datasets, so their fold changes are unstable and nuclei under-detect them further; the EAE 5K panel lacks *Hmgcs2*, *Hmgcl*, *Oxct1*, *Acat1* and MCT2, so ketogenesis / ketolysis scores cannot be computed there; in situ *Hcar2* inside non-myeloid segments is spill-over; and transcript abundance says nothing about BHB flux, which in brain is dominated by import from the circulation rather than local synthesis.
