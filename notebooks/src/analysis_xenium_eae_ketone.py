# %% [markdown]
# # Ketone-body handling in Xenium EAE: the genes that are on the 5K panel
#
# The central question of the project is ketone-body metabolism in the EAE spinal cord. The 5K panel has none of the committed ketone enzymes
# (*Hmgcs2*, *Hmgcl*, *Oxct1*, *Acat1*, *Aacs*), so this notebook works with what the panel does carry and asks the most that can be asked of it:
#
# | Group | Genes on the panel | What they can tell us |
# | --- | --- | --- |
# | Ketone enzyme | *Bdh1* | interconverts acetoacetate and BHB in both directions: the only ketone-body enzyme measurable in situ. Its *company* (fatty-acid oxidation genes → ketogenic side; MCT1 → export side) is the clue to its direction |
# | Transport | *Slc16a1* (MCT1), *Slc16a3* (MCT4), *Slc5a12* (SMCT2) | who can import / export BHB and lactate |
# | Sensing | *Hcar2* (GPR109A, BHB agonist), *Ffar3* (GPR41, BHB antagonist), *Ffar2*, *Nlrp3* (inflammasome inhibited by BHB) | who can respond to BHB |
# | Fatty-acid supply | *Cpt1a*, *Cpt1c*, *Cpt2*, *Slc25a20*, *Hadhb*, *Hadh*, *Echs1*, *Acaa2* | the mitochondrial β-oxidation chain that feeds ketogenesis |
# | Regulators | *Ppara*, *Ppargc1a*, *Pdk4*, *Fgf21*, *Klb* | the PPARα / FGF21 axis that switches ketogenesis on |
# | Acetate | *Acss1*, *Acss2* | the parallel two-carbon fuel |
#
# Design: 107 samples (92 EAE, 15 control), 891k cells, 27 annotated cell types, lesion-distance bins, RR and chronic EAE models with disease-course
# labels. Every comparison is pseudobulk per sample (≥20 cells per sample × cell type), Mann–Whitney across samples or paired Wilcoxon within samples.
# Sections: (1) who expresses what, (2) EAE vs control in every cell type, by model and by disease course, (3) lesion distance per cell type,
# (4) the disease-associated oligodendrocyte and astrocyte states, (5) what *Bdh1* keeps company with (direction), (6) myeloid ketone sensing
# (*Hcar2* / *Nlrp3*), (7) spatial neighbourhoods: are *Bdh1*⁺ / *Cpt1a*⁺ astrocytes and *Hcar2*⁺ myeloid cells enriched around DA and *C4b*-high
# oligodendrocytes, (8) maps, (9) interpretation.

# %%
import os, sys, gc, re, warnings
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
from oligometab import EAE5K_KETONE, EAE5K_KETONE_GENES, heat

sc.settings.verbosity = 0
sc.set_figure_params(dpi=80, frameon=False)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 300)

PATH = om.resolve_path("OLIGOMETAB_XENIUM_EAE_H5AD", "../../data/RREAE_5k_raw_only_integration_processed.h5ad")
print("Xenium EAE ->", PATH)
GROUP_OF = {g: k for k, gs in EAE5K_KETONE.items() for g in gs}
KG = EAE5K_KETONE_GENES
KEY = ["Bdh1", "Slc16a1", "Slc16a3", "Slc5a12", "Hcar2", "Ffar3", "Nlrp3", "Cpt1a", "Cpt2", "Hadhb", "Echs1", "Acaa2", "Ppara", "Ppargc1a", "Pdk4", "Fgf21", "Acss2"]
MYELOID = ["Microglia", "Macrophages", "Foamy Mic_Mac", "Activate Mic_Mac 1", "Activated Mic_Mac 2", "Efflux Mic_Mac", "Myeloid cells", "Proliferating microglia", "Mic_AST"]
ASTRO = ["Astrocytes", "DA astrocytes"]
OLIGO = ["Oligodendrocytes", "DA Oligodendrocytes", "Newly formed oligodendrocytes"]
os.makedirs("../../results", exist_ok=True)
R = {}


def save(name, df):
    R[name] = df; df.to_csv(f"../../results/eae5k_ketone_{name}.csv")


def pseudobulk_test(adata, genes, sample_col, group_col, a, b, ct_col="cell_type", cts=None, min_cells=20):
    """Per cell type: pseudobulk per sample, Mann-Whitney a vs b across samples, log2FC of means."""
    cts = cts or list(adata.obs[ct_col].unique())
    sub = adata[adata.obs[ct_col].isin(cts)]
    n = sub.obs.groupby([sample_col, ct_col], observed=True).size().rename("n_cells").reset_index(); n[[sample_col, ct_col]] = n[[sample_col, ct_col]].astype(str)
    pb = om.pseudobulk(sub, genes, sample_col, [ct_col, group_col]).merge(n, on=[sample_col, ct_col]).query("n_cells >= @min_cells")
    rows = []
    for ct in cts:
        p_ct = pb[pb[ct_col] == ct]
        for g in genes:
            ma, mr, p, na, nr = om.group_test(p_ct, g, group_col, a, b)
            rows.append({"cell_type": ct, "gene": g, "group": GROUP_OF.get(g), f"mean_{a}": ma, f"mean_{b}": mr, "log2FC": om.log2fc(ma, mr), "MWU_p": p, f"n_{a}": na, f"n_{b}": nr})
    return pd.DataFrame(rows), pb


def fc_heat(de, title, cts, genes, a, b, vmin=-2, vmax=2):
    d = de.pivot(index="cell_type", columns="gene", values="log2FC").reindex(index=cts, columns=genes); d = d[d.notna().any(axis=1)]
    p = de.pivot(index="cell_type", columns="gene", values="MWU_p").reindex(index=d.index, columns=d.columns)
    heat(d, title, "log2 fold change (* MWU p < 0.05)", fmt=".1f", cmap="RdBu_r", center=0, vmin=vmin, vmax=vmax,
         annot=(d.round(1).astype(str).replace("nan", "") + np.where(p.values < 0.05, "*", "")).values, row_labels=list(d.index))
    return d


def paired_genes(adata, genes, sample_col, mask_a, mask_b, la, lb, min_cells=20):
    df = om.expr_df(adata, genes); smp = adata.obs[sample_col].astype(str).values
    A, B = [], []
    for s_ in np.unique(smp):
        m = smp == s_
        if (m & mask_a).sum() >= min_cells and (m & mask_b).sum() >= min_cells:
            A.append(df[m & mask_a].mean()); B.append(df[m & mask_b].mean())
    A, B = pd.DataFrame(A), pd.DataFrame(B)
    return pd.DataFrame([{"gene": g, "group": GROUP_OF.get(g), f"mean_{la}": A[g].mean(), f"mean_{lb}": B[g].mean(), "log2FC": om.log2fc(A[g].mean(), B[g].mean()),
                          "n_samples": len(A), "frac_samples_up": float((A[g] > B[g]).mean()) if len(A) else np.nan,
                          "wilcoxon_p": stats.wilcoxon(A[g], B[g]).pvalue if len(A) >= 3 and not np.allclose(A[g], B[g]) else np.nan} for g in df.columns]).set_index("gene")


def bin_order(values):
    return sorted([v for v in values if str(v) not in ("nan", "None")], key=lambda b: (1, 0) if str(b).startswith(">") else (0, int(re.search(r"(\d+)", str(b)).group(1))))

# %%
from anndata import AnnData
from anndata.experimental import read_elem
with h5py.File(PATH, "r") as f:
    obs = read_elem(f["obs"]); var = read_elem(f["var"])
    shape = tuple(f["X"].attrs["shape"])
    X = sp.csr_matrix((f["X/data"][:].astype(np.float32), f["X/indices"][:], f["X/indptr"][:]), shape=shape)
    spatial = f["obsm/spatial"][:]
obs.index = obs.index.astype(str); var.index = var.index.astype(str)
ad = AnnData(X=X, obs=obs, var=var); ad.obsm["spatial"] = spatial; ad.obs_names_make_unique(); del X, obs, var
ad.obs["cell_type"] = ad.obs["cell_type"].astype(str).replace({"DA oligodendrocytes": "DA Oligodendrocytes", "Astrocyte": "Astrocytes"})
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
present = om.present(ad, KG); print("on panel:", present, "| not on panel:", om.missing(ad, KG))
KG = present; KEY = [g for g in KEY if g in KG]
for c in ["condition", "model", "course", "lesion_density_call"]:
    print(c, ad.obs[c].astype(str).value_counts().to_dict())
meta = ad.obs.groupby("sample_name", observed=True).agg(condition=("condition", "first"), model=("model", "first"), course=("course", "first"), n_cells=("condition", "size"))
meta["frac_lesion"] = ad.obs.groupby("sample_name", observed=True)["lesion_density_call"].apply(lambda s: (s.astype(str) != "non_lesion").mean())
print(meta.groupby(["condition", "model", "course"], observed=True).size())
vc = ad.obs["cell_type"].value_counts(); CTS = [c for c in vc.index if vc[c] >= 1500 and c != "unclear"]
P = om.pathway_scores(ad, {"FA supply for ketogenesis": ["Cpt1a", "Cpt2", "Slc25a20", "Hadhb", "Hadh", "Echs1", "Acaa2"]}, min_genes=3)
ad.obs["FA_supply_score"] = P["FA supply for ketogenesis"].values

# %% [markdown]
# ## 1. Who expresses what

# %%
det = om.detection_table(ad, KG, "cell_type", min_cells=1500).loc[CTS]; save("detection_by_celltype", det)
heat(det[KG], "Xenium EAE: fraction of cells with ≥1 transcript (all samples)", "fraction detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, row_labels=list(det.index))
mt = om.mean_table(ad, KG, "cell_type", min_cells=1500).loc[CTS]; save("mean_by_celltype", mt)
ctx = om.present(ad, ["C4b", "Plp1", "Aqp4", "Hexb", "Gfap", "Mbp"]); r = om.resolve(ad, KG + ctx)
sc.pl.dotplot(ad[ad.obs["cell_type"].isin(CTS)], var_names=[r[g] for g in KG + ctx], groupby="cell_type", standard_scale="var", color_map="Reds", figsize=(0.4 * (len(KG) + 4) + 2, 7), show=False, title="Xenium EAE: ketone-related genes by cell type"); plt.show()
# control tissue only: the healthy baseline
det_ctrl = om.detection_table(ad[ad.obs["condition"] == "CONTROL"], KG, "cell_type", min_cells=300); det_ctrl = det_ctrl.loc[[c for c in CTS if c in det_ctrl.index]]; save("detection_by_celltype_control_only", det_ctrl)
heat(det_ctrl[KG], "Xenium EAE, CONTROL samples only: fraction of cells with ≥1 transcript (healthy baseline)", "fraction detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, row_labels=list(det_ctrl.index))
fa = ad.obs.groupby("cell_type")["FA_supply_score"].mean().loc[CTS].sort_values(ascending=False); save("FA_supply_score_by_celltype", fa.to_frame())
print("fatty-acid-supply score by cell type:"); display(fa.round(3).to_frame().T)

# %% [markdown]
# ## 2. EAE vs control in every cell type; by model (RR / chronic) and by disease course

# %%
de, pb = pseudobulk_test(ad, KG, "sample_name", "condition", "EAE", "CONTROL", cts=CTS); save("EAE_vs_control", de)
fc_heat(de, "Xenium EAE: pseudobulk log2FC EAE vs CONTROL per cell type (* MWU p < 0.05 across samples)", CTS, KG, "EAE", "CONTROL")
sig = de[(de.MWU_p < 0.05)].sort_values(["gene", "log2FC"]); print(f"{len(sig)} significant cell type × gene contrasts:"); display(sig.round(3).reset_index(drop=True))
# by model
for model in sorted(ad.obs["model"].astype(str).unique()):
    m = ad[(ad.obs["model"].astype(str) == model) | (ad.obs["condition"] == "CONTROL")]
    de_m, _ = pseudobulk_test(m, KEY, "sample_name", "condition", "EAE", "CONTROL", cts=["Oligodendrocytes", "DA Oligodendrocytes", "OPCs", "Astrocytes", "DA astrocytes", "Microglia", "Neurons", "Endothelial cells"])
    de_m["model"] = model; save(f"EAE_vs_control_model_{model}", de_m)
    fc_heat(de_m, f"{model} EAE vs CONTROL, key genes (* MWU p < 0.05)", ["Oligodendrocytes", "DA Oligodendrocytes", "OPCs", "Astrocytes", "DA astrocytes", "Microglia", "Neurons", "Endothelial cells"], KEY, "EAE", "CONTROL")
    del m
# by disease course: pseudobulk per sample for the main cell types, plotted along course labels
courses = [c for c in ad.obs["course"].astype(str).unique() if c not in ("nan", "None")]
pb_c = pb[pb.cell_type.isin(["Oligodendrocytes", "DA Oligodendrocytes", "Astrocytes", "DA astrocytes", "Microglia", "Neurons"])].merge(meta[["course", "model", "frac_lesion"]], left_on="sample_name", right_index=True, suffixes=("", "_meta"))
tidy = pb_c.melt(id_vars=["sample_name", "cell_type", "condition", "course", "model", "frac_lesion"], value_vars=KEY, var_name="gene", value_name="mean_expr")
order_c = list(meta.groupby("course", observed=True)["frac_lesion"].mean().sort_values().index)
g = sns.catplot(data=tidy[tidy.gene.isin(["Bdh1", "Slc16a1", "Slc16a3", "Hcar2", "Cpt1a", "Ppara", "Acss2", "Nlrp3"])], x="course", y="mean_expr", hue="cell_type", col="gene", col_wrap=4, kind="point", order=order_c, dodge=0.4, errorbar=("ci", 95), linestyle="none", height=3.2, aspect=1.6, sharey=False)
for ax in g.axes.flat:
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=8)
g.set_titles("{col_name}"); plt.show()
course_tab = pb_c.groupby(["cell_type", "course"], observed=True)[KEY].mean().round(3); save("pseudobulk_by_course", course_tab); display(course_tab)
# sample-level: does Bdh1 / Slc16a1 in oligodendrocytes track lesion burden?
rows = []
for ct in ["Oligodendrocytes", "DA Oligodendrocytes", "Astrocytes", "DA astrocytes", "Microglia"]:
    p_ct = pb_c[(pb_c.cell_type == ct) & (pb_c.condition == "EAE")]
    for gname in KEY:
        if p_ct[gname].std() > 0:
            rho, p = stats.spearmanr(p_ct["frac_lesion"], p_ct[gname]); rows.append({"cell_type": ct, "gene": gname, "n_samples": len(p_ct), "spearman_rho_vs_lesion_fraction": rho, "p": p})
burden = pd.DataFrame(rows); save("gene_vs_lesion_burden_across_samples", burden)
print("Spearman correlation across EAE samples between per-sample mean expression and the sample's lesion fraction:"); display(burden.pivot(index="gene", columns="cell_type", values="spearman_rho_vs_lesion_fraction").round(2))

# %% [markdown]
# ## 3. Lesion distance, per cell type

# %%
m_eae = ((ad.obs["condition"] == "EAE") & ad.obs["lesion_distance_bin"].notna()).values
bins = bin_order(ad.obs.loc[m_eae, "lesion_distance_bin"].astype(str).unique())
s_all = ad[m_eae]
df = om.expr_df(s_all, KEY); ct_bin = df.groupby([s_all.obs["cell_type"].astype(str).values, s_all.obs["lesion_distance_bin"].astype(str).values]).mean()
n_ct_bin = s_all.obs.groupby([s_all.obs["cell_type"].astype(str).values, s_all.obs["lesion_distance_bin"].astype(str).values]).size(); ct_bin = ct_bin[n_ct_bin >= 50]
ct_bin.index.names = ["cell_type", "bin"]; save("lesion_distance_by_celltype", ct_bin)
for gname in ["Bdh1", "Slc16a1", "Slc16a3", "Hcar2", "Cpt1a", "Ppara", "Nlrp3", "Acss2"]:
    m = ct_bin[gname].unstack("bin").reindex(columns=bins); m = m.loc[[c for c in CTS if c in m.index]]
    far = m[bins[-1]]
    heat(m.sub(far, axis=0), f"{gname}: mean log-expression by cell type and lesion distance, relative to the >500 µm value of the same cell type", "difference from far-from-lesion", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6, row_labels=list(m.index))
fa_bin = s_all.obs.groupby([s_all.obs["cell_type"].astype(str).values, s_all.obs["lesion_distance_bin"].astype(str).values])["FA_supply_score"].mean(); fa_bin = fa_bin[n_ct_bin >= 50].unstack(1).reindex(columns=bins).loc[[c for c in CTS if c in fa_bin.index.get_level_values(0)]]
save("FA_supply_score_lesion_distance", fa_bin)
heat(fa_bin.sub(fa_bin[bins[-1]], axis=0), "Fatty-acid-supply score by cell type and lesion distance, relative to >500 µm", "score difference", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.3, vmax=0.3, row_labels=list(fa_bin.index))
del s_all, df

# %% [markdown]
# ## 4. The disease-associated oligodendrocyte and astrocyte states (paired within sample)

# %%
ct_ = ad.obs["cell_type"].values
c4b = om.expr(ad, "C4b"); is_ol = np.isin(ct_, ["Oligodendrocytes", "DA Oligodendrocytes"]); thr = np.quantile(c4b[is_ol & (c4b > 0)], 0.75)
tests = {"DA_vs_homeostatic_oligodendrocytes": (ct_ == "DA Oligodendrocytes", ct_ == "Oligodendrocytes", "DA", "homeo"),
         "C4bhigh_vs_C4bneg_oligodendrocytes": (is_ol & (c4b >= thr), is_ol & (c4b == 0), "C4b_high", "C4b_neg"),
         "DA_vs_homeostatic_astrocytes": (ct_ == "DA astrocytes", ct_ == "Astrocytes", "DA", "homeo"),
         "NFOL_vs_homeostatic_oligodendrocytes": (ct_ == "Newly formed oligodendrocytes", ct_ == "Oligodendrocytes", "NFOL", "homeo")}
fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
for ax, (name, (ma, mb, la, lb)) in zip(axes, tests.items()):
    t = paired_genes(ad, KG, "sample_name", ma, mb, la, lb); save(name, t)
    print(f"=== {name}"); display(t.round(4))
    cols = ["tab:red" if (p < 0.05) else "lightgrey" for p in t["wilcoxon_p"]]
    ax.barh(t.index, t["log2FC"], color=cols); ax.axvline(0, c="k", lw=0.8); ax.set_title(name.replace("_", " "), fontsize=9); ax.set_xlabel("log2FC (mean of per-sample means)")
plt.tight_layout(); plt.show()

# %% [markdown]
# ## 5. What does *Bdh1* keep company with? (direction of the reaction)
#
# *Bdh1* is bidirectional. If *Bdh1*⁺ cells of a type are enriched for the β-oxidation chain (*Cpt1a*, *Hadhb*, *Echs1*, *Acaa2*) and *Ppara*, the cell is set
# up to *make* BHB from fatty acids; if they are enriched for MCT1 and the acetate / glycolytic side, the picture is export or consumption. Tested inside
# astrocytes, oligodendrocytes, OPCs, microglia and neurons, in control and EAE separately.

# %%
partners = [g for g in ["Cpt1a", "Cpt2", "Slc25a20", "Hadhb", "Hadh", "Echs1", "Acaa2", "Ppara", "Ppargc1a", "Pdk4", "Slc16a1", "Slc16a3", "Acss2", "Hcar2", "Nlrp3", "C4b", "Gfap", "Plp1"] if g in ad.var_names]
rows = []
for ct in ["Astrocytes", "DA astrocytes", "Oligodendrocytes", "DA Oligodendrocytes", "Newly formed oligodendrocytes", "OPCs", "Microglia", "Neurons", "Endothelial cells"]:
    for cond in ["CONTROL", "EAE"]:
        sub = ad[(ad.obs["cell_type"] == ct) & (ad.obs["condition"] == cond)]
        if sub.n_obs < 300 or (om.expr(sub, "Bdh1") > 0).sum() < 50:
            continue
        co = om.coexpression_with_anchor(sub, "Bdh1", partners, min_pos=50); co["cell_type"] = ct; co["condition"] = cond; co["n_cells"] = sub.n_obs
        rows.append(co.reset_index()); del sub
BDH = pd.concat(rows, ignore_index=True); save("Bdh1_company", BDH)
for cond in ["CONTROL", "EAE"]:
    piv = BDH[BDH.condition == cond].pivot(index="cell_type", columns="gene", values="spearman_rho")[partners]
    heat(piv, f"{cond}: Spearman correlation of each gene with Bdh1 inside the cell type", "rho with Bdh1", fmt=".2f", cmap="RdBu_r", center=0, vmin=-0.2, vmax=0.2, row_labels=list(piv.index))
    piv2 = BDH[BDH.condition == cond].pivot(index="cell_type", columns="gene", values="enrichment(+/-)")[partners]
    heat(np.log2(piv2.replace(0, np.nan)), f"{cond}: detection enrichment in Bdh1+ vs Bdh1- cells (log2)", "log2 detection ratio", fmt=".1f", cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, row_labels=list(piv2.index))
fa_by_bdh = []
for ct in ["Astrocytes", "DA astrocytes", "Oligodendrocytes", "DA Oligodendrocytes", "OPCs", "Microglia", "Neurons"]:
    m = (ad.obs["cell_type"] == ct).values; b = om.expr(ad, "Bdh1") > 0
    fa_by_bdh.append({"cell_type": ct, "FA_supply_score_Bdh1+": ad.obs.loc[m & b, "FA_supply_score"].mean(), "FA_supply_score_Bdh1-": ad.obs.loc[m & ~b, "FA_supply_score"].mean(), "n_Bdh1+": int((m & b).sum())})
fa_by_bdh = pd.DataFrame(fa_by_bdh).set_index("cell_type"); fa_by_bdh["diff"] = fa_by_bdh.iloc[:, 0] - fa_by_bdh.iloc[:, 1]; save("FA_supply_score_by_Bdh1", fa_by_bdh); display(fa_by_bdh.round(3))

# %% [markdown]
# ## 6. Myeloid ketone sensing: *Hcar2* and *Nlrp3*

# %%
my = ad[ad.obs["cell_type"].isin(MYELOID)].copy()
d = om.detection_table(my, ["Hcar2", "Nlrp3", "Ffar2", "Ffar3", "Slc16a3", "Slc16a1", "Cpt1a", "Ppara", "Acss1", "Acss2", "Bdh1"], "cell_type", min_cells=300); save("myeloid_detection", d)
heat(d.drop(columns="n_cells"), "Myeloid populations: fraction detecting the ketone-sensing and transport genes", "fraction detected", fmt=".2f", cmap="Reds", center=None, vmin=0, vmax=None, row_labels=list(d.index))
cd = pd.DataFrame([{"cell_type": ct, "n": int((my.obs["cell_type"] == ct).sum()), "Hcar2+": (om.expr(my, "Hcar2")[(my.obs["cell_type"] == ct).values] > 0).mean(),
                    "Nlrp3+": (om.expr(my, "Nlrp3")[(my.obs["cell_type"] == ct).values] > 0).mean(),
                    "both": ((om.expr(my, "Hcar2") > 0) & (om.expr(my, "Nlrp3") > 0))[(my.obs["cell_type"] == ct).values].mean()} for ct in d.index]).set_index("cell_type")
cd["both_expected"] = cd["Hcar2+"] * cd["Nlrp3+"]; save("myeloid_Hcar2_Nlrp3_codetection", cd); display(cd.round(3))
de_my, pb_my = pseudobulk_test(my, ["Hcar2", "Nlrp3", "Ffar2", "Slc16a3", "Slc16a1", "Cpt1a", "Cpt2", "Hadhb", "Ppara", "Pdk4", "Acss1", "Acss2", "Bdh1"], "sample_name", "condition", "EAE", "CONTROL", cts=list(d.index)); save("myeloid_EAE_vs_control", de_my)
fc_heat(de_my, "Myeloid populations: pseudobulk log2FC EAE vs CONTROL (* MWU p < 0.05)", list(d.index), ["Hcar2", "Nlrp3", "Ffar2", "Slc16a3", "Slc16a1", "Cpt1a", "Cpt2", "Hadhb", "Ppara", "Pdk4", "Acss1", "Acss2", "Bdh1"], "EAE", "CONTROL")
m_e = ((my.obs["condition"] == "EAE") & my.obs["lesion_distance_bin"].notna()).values
ld = om.mean_table(my[m_e], ["Hcar2", "Nlrp3", "Slc16a3", "Cpt1a", "Ppara"], "lesion_distance_bin", min_cells=100).reindex(bins); save("myeloid_lesion_distance", ld)
print("all myeloid cells in EAE: mean log-expression by lesion distance"); display(ld.round(3))
hc = om.coexpression_with_anchor(my, "Hcar2", ["Nlrp3", "Slc16a3", "Cpt1a", "Ppara", "Pdk4", "Bdh1", "Acss2", "Cd36", "Lpl", "Plin2", "C4b"], min_pos=50); save("myeloid_Hcar2_company", hc)
print("inside myeloid cells: genes vs Hcar2"); display(hc[["frac_in_Hcar2-", "frac_in_Hcar2+", "enrichment(+/-)", "spearman_rho", "fisher_p"]].round(4))
del my

# %% [markdown]
# ## 7. Spatial neighbourhoods
#
# Are ketogenic-competent astrocytes (*Bdh1*⁺, *Cpt1a*⁺, *Bdh1*⁺*Cpt1a*⁺), MCT4⁺ astrocytes, *Hcar2*⁺ and *Nlrp3*⁺ myeloid cells enriched within 30 µm of
# (a) DA versus homeostatic oligodendrocytes and (b) *C4b*-high versus *C4b*-negative oligodendrocytes, in EAE samples, with the lesion-distance-matched control.

# %%
eae = ad[ad.obs["condition"] == "EAE"].copy()
ct_e = eae.obs["cell_type"].values; c4b_e = om.expr(eae, "C4b"); is_ol_e = np.isin(ct_e, ["Oligodendrocytes", "DA Oligodendrocytes"]); thr_e = np.quantile(c4b_e[is_ol_e & (c4b_e > 0)], 0.75)
sources = {"DA vs homeostatic oligodendrocytes": (ct_e == "DA Oligodendrocytes", ct_e == "Oligodendrocytes"),
           "C4b-high vs C4b-neg oligodendrocytes": (is_ol_e & (c4b_e >= thr_e), is_ol_e & (c4b_e == 0))}
astro_m = np.isin(ct_e, ASTRO); my_m = np.isin(ct_e, MYELOID)
targets = {"Bdh1+ astrocytes": (om.expr(eae, "Bdh1") > 0) & astro_m, "Cpt1a+ astrocytes": (om.expr(eae, "Cpt1a") > 0) & astro_m,
           "Bdh1+Cpt1a+ astrocytes": (om.expr(eae, "Bdh1") > 0) & (om.expr(eae, "Cpt1a") > 0) & astro_m, "Slc16a3+ astrocytes": (om.expr(eae, "Slc16a3") > 0) & astro_m,
           "Hcar2+ myeloid": (om.expr(eae, "Hcar2") > 0) & my_m, "Nlrp3+ myeloid": (om.expr(eae, "Nlrp3") > 0) & my_m, "Slc16a1+ endothelial": (om.expr(eae, "Slc16a1") > 0) & (ct_e == "Endothelial cells")}
strata = eae.obs["lesion_distance_bin"].astype(str).values
rows = []
for sname, (sa, sb) in sources.items():
    for tname, tgt in targets.items():
        j = om.compare_neighbourhoods(eae, "sample_name", sa, sb, tgt, "A", "B", radius=30.0, verbose=False)
        # lesion-distance-matched: per bin, per sample; median over bins within sample
        ratios = []
        for b in bins:
            m = strata == b
            A = om.neighbourhood_fraction(eae, "sample_name", sa & m, tgt, 30.0, min_source=20); B = om.neighbourhood_fraction(eae, "sample_name", sb & m, tgt, 30.0, min_source=20)
            jj = A[["frac_target_neighbours"]].join(B[["frac_target_neighbours"]], lsuffix="_A", rsuffix="_B", how="inner"); jj = jj[jj.iloc[:, 1] > 0]
            ratios += [(smp, r_ / rb) for smp, (r_, rb) in zip(jj.index, jj.values)]
        ps = pd.DataFrame(ratios, columns=["sample", "ratio"]).groupby("sample")["ratio"].median() if ratios else pd.Series(dtype=float)
        rows.append({"source": sname, "target": tname, "n_samples": len(j), "median_ratio": j["ratio"].median(), "frac_samples_gt1": (j["ratio"] > 1).mean(),
                     "wilcoxon_p": stats.wilcoxon(j.iloc[:, 0], j.iloc[:, 1]).pvalue if len(j) >= 3 else np.nan,
                     "lesion_matched_n": len(ps), "lesion_matched_median_ratio": ps.median() if len(ps) else np.nan,
                     "lesion_matched_wilcoxon_p": stats.wilcoxon(np.log(ps.clip(lower=1e-3))).pvalue if len(ps) >= 3 else np.nan})
NB = pd.DataFrame(rows); save("neighbourhoods", NB); display(NB.round(4))
plt.figure(figsize=(9, 3.8)); sns.barplot(data=NB, x="target", y="median_ratio", hue="source"); plt.axhline(1, ls="--", c="r"); plt.xticks(rotation=30, ha="right"); plt.ylabel("median neighbour-fraction ratio (30 µm)"); plt.xlabel(""); plt.tight_layout(); plt.show()
del eae

# %% [markdown]
# ## 8. Maps (one peak-EAE and one control sample)

# %%
peak = meta[(meta.condition == "EAE") & meta.course.astype(str).str.contains("peak", case=False)].sort_values("n_cells", ascending=False).index[0]
ctrl = meta[meta.condition == "CONTROL"].sort_values("n_cells", ascending=False).index[0]
for s_ in [peak, ctrl]:
    sub = ad[ad.obs["sample_name"].astype(str) == str(s_)].copy()
    genes = [g for g in ["C4b", "Bdh1", "Slc16a1", "Slc16a3", "Hcar2", "Nlrp3", "Cpt1a", "Ppara"] if g in sub.var_names]
    sc.pl.spatial(sub, color=genes, spot_size=20, cmap="magma", vmax="p99", ncols=4, show=False); plt.suptitle(f"{s_} ({meta.loc[s_, 'condition']}, {meta.loc[s_, 'course']})", y=1.02); plt.show()
    sc.pl.spatial(sub, color=["FA_supply_score"], spot_size=20, cmap="RdBu_r", vcenter=0, show=False); plt.suptitle(f"{s_}: fatty-acid-supply score", y=1.02); plt.show()
    del sub

# %% [markdown]
# ### Interpretation
#
# **What the panel allows.** The 5K panel has no committed ketone enzyme, so this is a study of *Bdh1* (bidirectional), the monocarboxylate transporters, the BHB sensors and the fatty-acid-oxidation supply chain, in 107 samples and 27 cell types. Read with that limit in mind, the in-situ data give five clear answers.
#
# **1. Ketogenic competence is astrocytic at baseline, not oligodendroglial or microglial.** In control tissue *Bdh1* is detected in 51 % of astrocytes, 45 % of newly formed oligodendrocytes, 29 % of ependymal cells, 27 % of neurons, 22 % of OPCs, 16 % of mature oligodendrocytes and endothelial cells and 8 % of microglia, and astrocytes also carry the highest β-oxidation chain (*Cpt1a* 41 %, *Hadhb* 50 %, *Echs1* 55 %, *Acaa2* 24 %) and the highest *Ppara* / *Ppargc1a*; their fatty-acid-supply score (0.31) is the highest of all cell types, followed by foamy myeloid cells (0.25) and DA astrocytes (0.18), with oligodendrocytes at zero and neurons at −0.27. Inside every cell type *Bdh1*⁺ cells are 1.3–3.8× enriched for the β-oxidation genes and for *Ppara* / *Ppargc1a* (oligodendrocytes: *Acaa2* 2.3×, *Cpt1a* 2.1×, *Ppara* 2.2×, *Ppargc1a* 2.4×; astrocytes ≈ 1.5× across the chain) and not for MCT4 or *Hcar2*, so *Bdh1* sits on the fatty-acid-oxidation, i.e. ketogenic, side of the reaction wherever it is expressed. MCT1 (*Slc16a1*) is oligodendroglial (40 %), endothelial (32 %) and astrocytic (31 %; 70 % of DA astrocytes); MCT4 (*Slc16a3*) is microglial (13 %); *Hcar2* (2 % of control microglia) and *Nlrp3* (20 %) are myeloid; SMCT2, *Ffar3*, *Ffar2*, *Fgf21* and *Klb* are at the detection floor.
#
# **2. EAE does not switch on a ketogenic program in any cell type; it dismantles it near lesions.** *Bdh1* falls toward lesions in every cell type (astrocytes −0.33, DA astrocytes −0.49, endothelial −0.28, neurons −0.23, newly formed oligodendrocytes −0.22, OPCs −0.13 log-expression units between the 0–10 µm and >500 µm bins), *Acss2* falls toward lesions and in 13 cell types overall, *Ppargc1a* falls in 10 and *Ppara* in 5 (astrocytes −0.35, OPCs −0.75 log2), and the astrocytic β-oxidation chain thins toward lesions (*Cpt1a* −0.29, fatty-acid-supply score −0.19). Across EAE samples, astrocytic *Acss2* (ρ −0.49) and *Ppargc1a* (−0.38) fall with lesion burden. The exceptions are oligodendrocytes, in which EAE raises *Bdh1* modestly (+0.21 log2 in homeostatic, +0.20 in DA oligodendrocytes; per-sample *Bdh1* tracks lesion burden, ρ 0.39), and myeloid cells and DA oligodendrocytes, in which *Cpt1a* / *Cpt2* rise (microglia *Cpt1a* +0.73, DA oligodendrocytes +0.55 / +0.71, both increasing toward lesions and with lesion burden, ρ 0.4–0.6). *Bdh1* is flat across the disease course in astrocytes (1.2–1.4) and oligodendrocytes (0.35–0.50). The RR model shows every effect more strongly than the chronic model (e.g. astrocyte MCT1 +1.0 vs +0.3 log2, microglial *Bdh1* −0.72 vs −0.32).
#
# **3. What EAE turns on is monocarboxylate transport.** MCT4 rises in 14 cell types (astrocytes +1.3, DA astrocytes +2.5, endothelial +1.7, microglia +1.0 log2) and toward lesions in all of them (DA oligodendrocytes +0.30, OPCs +0.32, endothelial +0.29, macrophages +0.29), and MCT1 rises in 10 (astrocytes +0.8, DA oligodendrocytes +0.6, newly formed oligodendrocytes +0.7, OPCs +0.4). Astrocytic MCT1 is the single ketone-related gene that best tracks the disease: it climbs from 0.7 in CFA controls to 1.5–2.1 at the RR peaks and returns to 0.9–1.3 in remission, and its per-sample level correlates with lesion fraction at ρ 0.68 (DA astrocytes 0.59, microglia 0.45). Endothelial MCT1 goes the other way (−0.33 toward lesions). So EAE re-wires lactate / ketone shuttling (astrocytes and myeloid cells export through MCT4, astrocytes and oligodendrocytes gain MCT1) without changing local ketone synthesis.
#
# **4. The disease-associated oligodendrocyte is not a ketone-using state; the newly formed oligodendrocyte is.** Paired within 107 samples, DA oligodendrocytes have less *Bdh1* (−0.42 log2; higher in only 4 % of samples), less MCT1 (−0.21), less of the whole β-oxidation chain (*Cpt1c* −0.60, *Cpt2* −0.70, *Hadhb* −0.50, *Hadh* −0.41, *Echs1* −0.35, *Slc25a20* −0.44), less *Ppara* (−0.60), *Ppargc1a* (−1.3), *Pdk4* (−0.53) and *Acss2* (−0.69) than homeostatic oligodendrocytes, and more MCT4 (+1.2). *C4b*-high lineage cells show the same *Bdh1* / β-oxidation deficit with higher *Pdk4* (+0.77, 87 % of samples), MCT1 (+0.25) and *Acss2* (+0.46). Newly formed oligodendrocytes are the opposite: *Bdh1* +1.4 log2 (higher in 100 % of samples), *Cpt1a* +0.9, *Cpt1c* +0.7, *Slc25a20* +0.6, *Acss1* +0.6, MCT4 +1.8, MCT1 −1.0 and *Pdk4* −2.2 relative to mature oligodendrocytes: the remyelinating stage is the lineage's fatty-acid- and ketone-handling stage, in line with the progenitor ketolysis seen in the sorted Falcão cells. DA astrocytes likewise lose *Bdh1* (−0.34), *Acss1* (−1.2), *Hadhb* (−0.67), *Acaa2* (−0.85), *Ppara* (−1.0) and *Ppargc1a* (−0.93) and gain MCT1 (+1.2, 100 % of samples) and MCT4 (+2.5).
#
# **5. BHB sensing belongs to myeloid cells, and those cells surround DA oligodendrocytes.** *Hcar2* is in 8 % of microglia, 13 % of activated myeloid cells and 20–21 % of proliferating microglia and infiltrating myeloid cells; *Nlrp3* in 22 % of microglia, 28 % of macrophages and 45 % of infiltrating myeloid cells; the two co-occur 1.5–1.9× more often than expected, and *Hcar2*⁺ myeloid cells are 1.9× enriched for *Nlrp3*, 2.4× for *Cd36* and 1.4× for *Lpl* and MCT4, i.e. the lipid-laden phagocyte. EAE raises microglial *Hcar2* four-fold (+2.1 log2; from 0.05 in controls to 0.27–0.47 at onset and peak, 0.16 in remission), *Ffar2* (+1.5), *Cpt1a* (+0.7) and MCT4 (+1.0) and lowers *Bdh1* (−0.6) and *Acss2* (−0.6); *Nlrp3* per cell is unchanged, its tissue-level rise is recruitment. Toward lesions, microglial *Hcar2* falls (0.26 → 0.06) while macrophage *Hcar2* rises (0.08 → 0.17), so the BHB-responsive cells at the lesion are infiltrating and activated myeloid cells, not resident microglia. The apparent *Hcar2* and *Nlrp3* "induction" in 12–15 non-myeloid cell types (up to +5 log2 in newly formed oligodendrocytes) is spill-over from these cells into neighbouring segments and should not be read as intrinsic. Spatially, DA oligodendrocytes have *fewer* *Bdh1*⁺ (0.82), *Cpt1a*⁺ (0.77) and *Bdh1*⁺*Cpt1a*⁺ (0.69) astrocytes and fewer MCT1⁺ endothelial cells (0.77) within 30 µm than homeostatic oligodendrocytes of the same sample, and *more* MCT4⁺ astrocytes (1.9), *Hcar2*⁺ (1.6) and *Nlrp3*⁺ (1.6) myeloid cells; all of this survives lesion-distance matching (0.82–0.90 and 1.2–1.5). Around *C4b*-high versus *C4b*-negative oligodendrocytes the astrocyte targets are neutral (0.9–1.0) while MCT4⁺ astrocytes (2.3; matched 1.75) and *Hcar2*⁺ / *Nlrp3*⁺ myeloid cells (1.7–1.8; matched 1.7) remain enriched. The DA oligodendrocyte niche is therefore lactate-rich and ketogenesis-poor.
#
# **What this means for the central question.** At the transcript level the EAE spinal cord does not up-regulate local ketone production; the ketogenic-competent cells (astrocytes) lose *Bdh1*, PPARα / PGC1α and β-oxidation near lesions, and the DA oligodendrocyte is a lactate-exporting, lipid-importing cell rather than a ketone-consuming one. If BHB is beneficial in EAE, these data point to exogenous supply (diet, supplementation) rather than induced synthesis, and to two cellular targets: *Hcar2*⁺ / *Nlrp3*⁺ myeloid cells that accumulate around DA oligodendrocytes (anti-inflammatory signalling), and cells that gain MCT1 in disease (astrocytes, DA and newly formed oligodendrocytes) for import as fuel, with the newly formed oligodendrocyte as the lineage stage best equipped to use it. Confirming utilisation needs *Oxct1* / *Aacs* (a custom Xenium add-on or the nuclei object with whole-transcriptome data) and BHB measurements.
#
# **Caveats.** No *Hmgcs2*, *Oxct1*, *Acat1*, *Aacs* or MCT2 on the panel; *Hcar2*, *Nlrp3*, *Fgf21*, *Klb*, SMCT2, *Ffar2* / *Ffar3* are near the detection floor and *Hcar2* / *Nlrp3* in non-myeloid segments are spill-over; *Bdh1* direction is inferred from its company, not measured; controls are 15 CFA-immunised animals; transcript abundance is not flux.
