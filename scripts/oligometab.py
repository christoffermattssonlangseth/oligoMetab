"""Gene panel and shared helpers for the OligoMetab analyses.

OligoMetab asks the questions that OligoC4b asked about the complement system of the *metabolic* machinery of
oligodendrocytes: ketone-body metabolism, glycolysis and lactate handling, the TCA cycle and oxidative
phosphorylation, fatty-acid oxidation, lipid and cholesterol synthesis, and the nutrient-sensing regulators that
sit on top of them. The datasets are the ones OligoC4b already harmonised (``OLIGOC4B_PUBLIC_PROCESSED_DIR``), so
this module contains no loaders, only:

- ``PANEL`` / ``ALL_GENES`` / ``FOCUS``   the metabolic gene panel, grouped by pathway (mouse symbols)
- ``MOUSE_TO_HUMAN`` / ``ALIASES``         symbol handling across species and annotation versions
- expression / detection / pseudobulk / co-expression / pathway-score / spatial-neighbourhood helpers shared by the
  two analysis notebooks

Every object produced by OligoC4b has raw counts in ``layers["counts"]``, log-normalised ``X`` and a harmonised
``obs`` (``dataset``, ``species``, ``modality``, ``sample``, ``group``, ``group_ref``, ``cell_type_coarse``,
``cell_type_original``); see ``docs/public_datasets.md`` in OligoC4b.
"""
from __future__ import annotations

import glob
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

try:  # scanpy is only needed for pathway scores and plots
    import scanpy as sc
except Exception:  # pragma: no cover
    sc = None

# --------------------------------------------------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------------------------------------------------

PROCESSED_DIR = os.getenv("OLIGOMETAB_PUBLIC_PROCESSED_DIR") or os.getenv("OLIGOC4B_PUBLIC_PROCESSED_DIR", "../../data/public/processed")

COARSE_TYPES = ["Oligodendrocyte", "OPC", "Microglia", "Astrocyte", "Neuron", "Endothelial", "Vascular/Fibroblast",
                "Immune (lymphoid/myeloid)", "Ependymal", "Other"]

# --------------------------------------------------------------------------------------------------------------------
# the metabolic panel (mouse symbols)
# --------------------------------------------------------------------------------------------------------------------

PANEL: Dict[str, List[str]] = {
    # ketone bodies: synthesis (Hmgcs2 is rate-limiting), utilisation (Oxct1 = SCOT, Bdh1, Acat1) and the ketone / SCFA receptors
    "Ketone body synthesis":            ["Hmgcs2", "Hmgcl", "Bdh1", "Bdh2"],
    "Ketone body utilisation":          ["Oxct1", "Acat1", "Acat2", "Acss1", "Acss2"],
    "Ketone / SCFA receptors":          ["Hcar2", "Ffar3", "Ffar2"],
    # monocarboxylate transport: Slc16a1 (MCT1, oligodendrocyte lactate/ketone export), Slc16a7 (MCT2, neurons),
    # Slc16a3 (MCT4, astrocytes), Slc16a6 (MCT7, ketone body exporter)
    "Monocarboxylate transporters":     ["Slc16a1", "Slc16a7", "Slc16a3", "Slc16a6", "Bsg"],
    # glucose uptake and glycolysis
    "Glucose transporters":             ["Slc2a1", "Slc2a3", "Slc2a4", "Slc2a8"],
    "Glycolysis":                       ["Hk1", "Hk2", "Hk3", "Gpi1", "Pfkl", "Pfkm", "Pfkp", "Pfkfb2", "Pfkfb3", "Aldoa", "Aldoc",
                                         "Tpi1", "Gapdh", "Pgk1", "Pgam1", "Eno1", "Eno2", "Pkm", "Ldha", "Ldhb"],
    "Pyruvate fate":                    ["Pdha1", "Pdhb", "Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Mpc1", "Mpc2", "Pcx", "Pck2"],
    "Pentose phosphate pathway":        ["G6pdx", "Pgd", "Tkt", "Taldo1"],
    "Glycogen":                         ["Gys1", "Pygb", "Gbe1", "Agl"],
    # mitochondrial oxidation
    "TCA cycle":                        ["Cs", "Aco2", "Idh2", "Idh3a", "Ogdh", "Dlst", "Sucla2", "Suclg1", "Sdha", "Sdhb", "Fh1", "Mdh1", "Mdh2"],
    "Anaplerosis / amino acids":        ["Got1", "Got2", "Glud1", "Gls", "Glul", "Idh1", "Slc25a1", "Slc25a11"],
    "Oxidative phosphorylation":        ["Ndufa4", "Ndufs1", "Ndufv1", "Ndufb8", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Cycs", "Cox4i1", "Cox5a",
                                         "Cox6c", "Cox7c", "Atp5a1", "Atp5b", "Atp5o"],
    "Mitochondrial biogenesis / dynamics": ["Ppargc1a", "Tfam", "Nrf1", "Opa1", "Mfn2", "Dnm1l", "Pink1", "Prkn"],
    # fatty acids
    "Fatty acid uptake / activation":   ["Cd36", "Slc27a1", "Fabp5", "Fabp7", "Acsl1", "Acsl3", "Acsl6", "Acsbg1"],
    "Fatty acid beta-oxidation":        ["Cpt1a", "Cpt1c", "Cpt2", "Slc25a20", "Acadm", "Acadl", "Acadvl", "Acads", "Hadha", "Hadhb",
                                         "Echs1", "Acaa2", "Decr1", "Eci1", "Acox1", "Ehhadh"],
    "Lipid synthesis":                  ["Acly", "Acaca", "Fasn", "Elovl1", "Elovl5", "Elovl6", "Scd1", "Scd2", "Srebf1", "Mlxipl"],
    "Cholesterol synthesis":            ["Hmgcs1", "Hmgcr", "Mvk", "Fdps", "Fdft1", "Sqle", "Lss", "Cyp51", "Dhcr24", "Dhcr7", "Srebf2", "Insig1"],
    "Cholesterol / lipid transport":    ["Ldlr", "Abca1", "Abca2", "Abcg1", "Apoe", "Apod", "Lpl", "Nr1h2", "Nr1h3", "Plin2", "Plin3", "Plin4"],
    "Myelin lipid synthesis":           ["Ugt8a", "Gal3st1", "Cers2", "Sptlc1", "Pigt", "Fa2h"],
    # regulators
    "Nutrient sensing / regulators":    ["Hif1a", "Epas1", "Prkaa1", "Prkaa2", "Mtor", "Rptor", "Tsc2", "Ddit4", "Txnip", "Sirt1", "Sirt3",
                                         "Ppara", "Ppard", "Pparg", "Nfe2l2"],
}
ALL_GENES: List[str] = [g for gs in PANEL.values() for g in gs]
PATHWAY_OF: Dict[str, str] = {g: k for k, gs in PANEL.items() for g in gs}

# the genes the question is about: one or two per arm, shown in every summary figure
FOCUS: List[str] = [
    "Hmgcs2", "Bdh1", "Oxct1", "Acat1", "Hcar2",           # ketone bodies
    "Slc16a1", "Slc16a7", "Slc16a3",                       # lactate / ketone shuttles
    "Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkp", "Aldoc", "Pkm", "Ldha", "Ldhb", "Pdk1",   # glycolysis and pyruvate fate
    "Cs", "Sdha", "Cox4i1", "Atp5a1", "Ppargc1a",          # TCA / OXPHOS
    "Cpt1a", "Acadm", "Hadha",                             # beta-oxidation
    "Fasn", "Hmgcr", "Srebf2", "Plin2", "Apoe",            # lipid / cholesterol synthesis and storage
    "Hif1a", "Txnip", "Ddit4",                             # sensing
]

# short pathway groups for pathway-score plots (the full panel is too fine for one figure)
SCORE_SETS: Dict[str, List[str]] = {
    "Ketone synthesis":     PANEL["Ketone body synthesis"],
    "Ketone utilisation":   PANEL["Ketone body utilisation"],
    "Glycolysis":           PANEL["Glucose transporters"] + PANEL["Glycolysis"],
    "Pyruvate to lactate":  ["Ldha", "Ldhb", "Pdk1", "Pdk2", "Pdk3", "Pdk4", "Slc16a1", "Slc16a3", "Slc16a7"],
    "Pentose phosphate":    PANEL["Pentose phosphate pathway"],
    "TCA cycle":            PANEL["TCA cycle"],
    "OXPHOS":               PANEL["Oxidative phosphorylation"],
    "Beta-oxidation":       PANEL["Fatty acid beta-oxidation"],
    "Lipid synthesis":      PANEL["Lipid synthesis"],
    "Cholesterol synthesis": PANEL["Cholesterol synthesis"],
    "Lipid transport / storage": PANEL["Cholesterol / lipid transport"],
    "Myelin lipids":        PANEL["Myelin lipid synthesis"],
}

# ketone-body metabolism in depth (used by analysis_ketone_metabolism.ipynb); a superset of the two ketone pathway groups above
KETONE: Dict[str, List[str]] = {
    "Synthesis (ketogenesis)":        ["Hmgcs2", "Hmgcl", "Hmgcll1", "Bdh1", "Bdh2"],          # Hmgcs2 is rate-limiting; Bdh1 interconverts AcAc <-> BHB (both directions)
    "Utilisation (ketolysis)":        ["Oxct1", "Acat1", "Bdh1"],                              # Oxct1 (SCOT) is the committed ketolytic step, absent from liver
    "Acetate activation":             ["Acss1", "Acss2"],
    "Ketone / lactate transport":     ["Slc16a1", "Slc16a7", "Slc16a3", "Slc16a6", "Slc5a8", "Bsg"],   # MCT1 (oligodendrocytes, endothelium), MCT2 (neurons), MCT4 (astrocytes), MCT7, SMCT1
    "Receptors":                      ["Hcar2", "Ffar3", "Ffar2"],                             # BHB agonist Hcar2 (GPR109A); Ffar3 (GPR41) is antagonised by BHB
    "Fatty-acid supply / regulators": ["Cpt1a", "Cpt2", "Acadm", "Hadha", "Acaa2", "Ppara", "Ppargc1a", "Fgf21"],
}
KETONE_GENES: List[str] = list(dict.fromkeys(g for gs in KETONE.values() for g in gs))
KETONE_SCORE_SETS: Dict[str, List[str]] = {
    "Ketogenesis": ["Hmgcs2", "Hmgcl", "Hmgcll1", "Bdh1", "Bdh2"],
    "Ketolysis": ["Oxct1", "Acat1", "Bdh1", "Slc16a1", "Slc16a7"],
    "FA supply for ketogenesis": ["Cpt1a", "Cpt2", "Acadm", "Hadha", "Acaa2", "Ppara"],
}
for _g in ["Hmgcll1", "Slc5a8", "Fgf21"]:
    PATHWAY_OF.setdefault(_g, "ketone extras")

# context genes used to check cell identity and to anchor on the C4b program described in OligoC4b
CONTEXT: List[str] = ["C4b", "Serpina3n", "Plp1", "Mbp", "Hexb", "Aqp4", "Gfap", "Itgam"]

# mouse -> human symbols where upper-casing is wrong or ambiguous. None = no one-to-one human ortholog (gene skipped).
MOUSE_TO_HUMAN: Dict[str, Optional[str]] = {
    "G6pdx": "G6PD", "Gpi1": "GPI", "Fh1": "FH", "Cyp51": "CYP51A1", "Scd1": "SCD", "Scd2": None, "Ugt8a": "UGT8", "Pcx": "PC",
    "Atp5a1": "ATP5F1A", "Atp5b": "ATP5F1B", "Atp5o": "ATP5PO", "Hcar2": "HCAR2", "Ffar3": "FFAR3",
    "Acsbg1": "ACSBG1", "Slc16a6": "SLC16A6", "Nr1h2": "NR1H2", "Nr1h3": "NR1H3", "Prkn": "PRKN",
    "C4b": "C4B", "Serpina3n": "SERPINA3", "Plp1": "PLP1", "Mbp": "MBP", "Hexb": "HEXB", "Aqp4": "AQP4", "Gfap": "GFAP", "Itgam": "ITGAM",
    "H2-D1": "HLA-A", "H2-K1": "HLA-B", "B2m": "B2M",
}

# alternative symbols for the same gene (annotation-version differences); tried in order after the primary symbol
ALIASES: Dict[str, List[str]] = {
    "Atp5a1": ["Atp5f1a"], "Atp5b": ["Atp5f1b"], "Atp5o": ["Atp5po"], "Hcar2": ["Gpr109a", "Niacr1"],
    "ATP5F1A": ["ATP5A1"], "ATP5F1B": ["ATP5B"], "ATP5PO": ["ATP5O"], "Slc16a6": ["Mct7"], "Acsbg1": ["Bgm"],
    "Ppargc1a": ["Pgc1a"], "Eci1": ["Dci"], "Ehhadh": ["Lbp"], "Idh3a": ["Idh3a"],
}


def sym(species: str, g: str) -> Optional[str]:
    """Mouse symbol -> symbol for ``species`` (None when the gene has no human ortholog)."""
    if species == "human":
        return MOUSE_TO_HUMAN.get(g, g.upper()) if g in MOUSE_TO_HUMAN else g.upper()
    return g


def resolve(adata, genes: Sequence[str]) -> Dict[str, Optional[str]]:
    """Map mouse panel symbols to the var_names actually present in ``adata`` (species from obs, aliases tried)."""
    species = str(adata.obs["species"].iloc[0]) if "species" in adata.obs else "mouse"
    names = set(adata.var_names)
    out: Dict[str, Optional[str]] = {}
    for g in genes:
        s = sym(species, g)
        cands = [] if s is None else [s] + ALIASES.get(s, []) + ALIASES.get(g, [])
        if species == "human":
            cands += [a.upper() for a in ALIASES.get(g, [])]
        out[g] = next((c for c in cands if c in names), None)
    return out


def present(adata, genes: Sequence[str]) -> List[str]:
    r = resolve(adata, genes)
    return [g for g in genes if r[g] is not None]


def missing(adata, genes: Sequence[str]) -> List[str]:
    r = resolve(adata, genes)
    return [g for g in genes if r[g] is None]


# --------------------------------------------------------------------------------------------------------------------
# expression access
# --------------------------------------------------------------------------------------------------------------------

def expr(adata, gene: str) -> np.ndarray:
    """1-D array of ``X`` for one gene given as a MOUSE symbol (resolved to the dataset's symbol)."""
    v = resolve(adata, [gene])[gene]
    if v is None:
        return np.full(adata.n_obs, np.nan, dtype=np.float32)
    x = adata[:, v].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()


def expr_df(adata, genes: Sequence[str]) -> pd.DataFrame:
    """cells x genes DataFrame of ``X``; columns are MOUSE symbols, unavailable genes are dropped."""
    r = resolve(adata, genes)
    keep = [g for g in genes if r[g] is not None]
    if not keep:
        return pd.DataFrame(index=adata.obs_names)
    X = adata[:, [r[g] for g in keep]].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return pd.DataFrame(X, index=adata.obs_names, columns=keep)


def detection_table(adata, genes: Sequence[str], groupby: str, min_cells: int = 30) -> pd.DataFrame:
    """Fraction of cells with X > 0 per group (rows) and gene (columns), plus n_cells."""
    df = expr_df(adata, genes)
    g = adata.obs[groupby].astype(str).values
    rows = []
    for grp, sub in df.groupby(g):
        if len(sub) < min_cells:
            continue
        row = {"group": grp, "n_cells": len(sub)}
        row.update({gene: float((sub[gene] > 0).mean()) for gene in df.columns})
        rows.append(row)
    return pd.DataFrame(rows).set_index("group") if rows else pd.DataFrame()


def mean_table(adata, genes: Sequence[str], groupby: str, min_cells: int = 30) -> pd.DataFrame:
    """Mean log-expression per group and gene."""
    df = expr_df(adata, genes)
    g = adata.obs[groupby].astype(str).values
    out = df.groupby(g).mean()
    n = pd.Series(g).value_counts()
    return out.loc[n[n >= min_cells].index.intersection(out.index)]


def counts_check(adata, genes: Sequence[str]) -> pd.DataFrame:
    """Total UMIs and detecting cells per gene from layers['counts'] (or X)."""
    r = resolve(adata, genes)
    rows = []
    for g in genes:
        v = r[g]
        if v is None:
            rows.append({"gene": g, "symbol": None, "total_UMIs": np.nan, "cells_detected": np.nan})
            continue
        L = adata[:, v].layers["counts"] if "counts" in adata.layers else adata[:, v].X
        L = L.toarray().ravel() if sp.issparse(L) else np.asarray(L).ravel()
        rows.append({"gene": g, "symbol": v, "total_UMIs": int(L.sum()), "cells_detected": int((L > 0).sum())})
    return pd.DataFrame(rows).set_index("gene")


# --------------------------------------------------------------------------------------------------------------------
# pseudobulk and group tests
# --------------------------------------------------------------------------------------------------------------------

def pseudobulk(adata, genes: Sequence[str], sample_col: str = "sample", extra: Sequence[str] = ("group",),
               values: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """One row per sample (x extra keys) with the mean of each gene (or of ``values`` columns) over its cells."""
    df = expr_df(adata, genes) if values is None else values.loc[adata.obs_names]
    keys = adata.obs[[sample_col] + list(extra)].astype(str)
    return df.join(keys).groupby([sample_col] + list(extra), observed=True).mean().reset_index()


def group_test(pb: pd.DataFrame, col: str, group_col: str, a: str, b: str) -> Tuple[float, float, float, int, int]:
    """Mean in a, mean in b, Mann-Whitney p (needs >= 3 samples per side, else nan), n_a, n_b."""
    xa, xb = pb.loc[pb[group_col] == a, col].dropna(), pb.loc[pb[group_col] == b, col].dropna()
    p = stats.mannwhitneyu(xa, xb, alternative="two-sided").pvalue if (len(xa) >= 3 and len(xb) >= 3) else np.nan
    return float(xa.mean()), float(xb.mean()), float(p), len(xa), len(xb)


def alt_groups(adata) -> Tuple[str, List[str]]:
    """Reference level and the alternative levels of ``obs['group']``."""
    ref = str(adata.obs["group_ref"].iloc[0])
    return ref, [g for g in adata.obs["group"].astype(str).unique() if g not in (ref, "unknown", "intermediate", "nan")]


def log2fc(ma: float, mr: float, eps: float = 1e-3) -> float:
    return float(np.log2((ma + eps) / (mr + eps)))


# --------------------------------------------------------------------------------------------------------------------
# pathway scores
# --------------------------------------------------------------------------------------------------------------------

def pathway_scores(adata, sets: Optional[Dict[str, List[str]]] = None, min_genes: int = 3, ctrl_size: int = 50) -> pd.DataFrame:
    """scanpy ``score_genes`` for each pathway (mean of the set minus a size-matched control set), cells x pathways.

    Scores are on the log-expression scale and can be negative; compare them by difference, not fold change.
    """
    if sc is None:
        raise ImportError("scanpy is required for pathway_scores")
    sets = sets or SCORE_SETS
    out = {}
    for name, genes in sets.items():
        r = resolve(adata, genes)
        vs = [r[g] for g in genes if r[g] is not None]
        if len(vs) < min_genes:
            continue
        key = f"score::{name}"
        sc.tl.score_genes(adata, vs, score_name=key, use_raw=False, ctrl_size=max(ctrl_size, len(vs)))
        out[name] = adata.obs.pop(key).astype(float)
    return pd.DataFrame(out, index=adata.obs_names)


# --------------------------------------------------------------------------------------------------------------------
# co-expression with an anchor gene (C4b) and genome-wide ranks
# --------------------------------------------------------------------------------------------------------------------

def coexpression_with_anchor(adata, anchor: str, genes: Sequence[str], high_q: float = 0.75, min_pos: int = 20) -> pd.DataFrame:
    """Inside one population: detection of each gene in anchor+ vs anchor- cells, Fisher p, Spearman rho with the anchor."""
    df = expr_df(adata, [anchor] + list(genes))
    if anchor not in df.columns or (df[anchor] > 0).sum() < min_pos:
        return pd.DataFrame()
    a = df[anchor]
    pos = a > 0
    thr = np.quantile(a[pos], high_q)
    rows = []
    for g in df.columns:
        if g == anchor:
            continue
        x = df[g]
        f_neg, f_pos = float((x[~pos] > 0).mean()), float((x[pos] > 0).mean())
        rho = stats.spearmanr(a, x).correlation if x.std() > 0 else np.nan
        tab = [[int((x[pos] > 0).sum()), int((x[pos] <= 0).sum())], [int((x[~pos] > 0).sum()), int((x[~pos] <= 0).sum())]]
        try:
            p = stats.fisher_exact(tab)[1]
        except Exception:
            p = np.nan
        rows.append({"gene": g, "pathway": PATHWAY_OF.get(g, "context"), f"frac_in_{anchor}-": f_neg, f"frac_in_{anchor}+": f_pos,
                     f"frac_in_{anchor}-high": float((x[a >= thr] > 0).mean()),
                     f"mean_in_{anchor}-": float(x[~pos].mean()), f"mean_in_{anchor}+": float(x[pos].mean()),
                     "enrichment(+/-)": f_pos / f_neg if f_neg > 0 else np.nan, "spearman_rho": rho, "fisher_p": p})
    return pd.DataFrame(rows).set_index("gene")


def genome_wide_rho(adata, anchor_var: str, min_frac: float = 0.01, max_cells: int = 40000, seed: int = 0, chunk: int = 2000) -> pd.Series:
    """Spearman correlation of every detected gene with ``anchor_var`` (a var_name), sorted descending.

    Genes detected in more than ``min_frac`` of cells are ranked column-wise in chunks (average ranks for ties, as
    ``pandas.rank``), standardised and correlated with the ranked anchor by a matrix product, so memory stays at one
    chunk (cells x ``chunk`` genes) instead of the full dense matrix. Cells are subsampled to ``max_cells``.
    """
    from scipy.stats import rankdata
    if adata.n_obs > max_cells:
        adata = adata[np.random.RandomState(seed).choice(adata.n_obs, max_cells, replace=False)]
    X = adata.X
    n = X.shape[0]
    if sp.issparse(X):
        X = X.tocsc()
        frac = np.asarray((X > 0).sum(axis=0)).ravel() / n
    else:
        X = np.asarray(X)
        frac = (X > 0).mean(axis=0)
    names = np.asarray(adata.var_names)
    ai = int(np.where(names == anchor_var)[0][0])
    keep = np.where((frac > min_frac) | (np.arange(len(names)) == ai))[0]
    a = X[:, ai].toarray().ravel() if sp.issparse(X) else X[:, ai]
    a = rankdata(a).astype(np.float64)
    a = (a - a.mean()) / (a.std() + 1e-12)
    out = np.empty(len(keep), dtype=np.float64)
    for start in range(0, len(keep), chunk):
        cols = keep[start:start + chunk]
        D = X[:, cols].toarray() if sp.issparse(X) else X[:, cols]
        R = rankdata(D, axis=0).astype(np.float64)
        R -= R.mean(axis=0)
        sd = R.std(axis=0)
        sd[sd == 0] = np.nan
        out[start:start + len(cols)] = (R.T @ a) / n / sd
    rho = pd.Series(out, index=names[keep]).drop(anchor_var, errors="ignore").dropna().sort_values(ascending=False)
    return rho


def rank_panel(rho: pd.Series, adata, genes: Sequence[str]) -> pd.DataFrame:
    """Rank (1 = most anchor-correlated) and rho for each panel gene (mouse symbols) among all ranked genes."""
    r = resolve(adata, genes)
    ranks = pd.Series(np.arange(1, len(rho) + 1), index=rho.index)
    rows = []
    for g in genes:
        v = r[g]
        if v is not None and v in rho.index:
            rows.append({"gene": g, "pathway": PATHWAY_OF.get(g, "context"), "symbol": v, "spearman_rho": rho[v],
                         "rank": int(ranks[v]), "n_ranked": len(rho)})
    return pd.DataFrame(rows).set_index("gene")


def _mwu_p(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Mann-Whitney p with the normal approximation.

    SciPy's default ``method='auto'`` switches to the exact distribution when one sample has <= 8 values, which for a
    small gene set against ~15,000 other genes takes minutes per test; the asymptotic p is what we want here.
    """
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").pvalue)
    except TypeError:  # scipy < 1.7 has no `method` and is asymptotic by default
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def pathway_rank_enrichment(rho: pd.Series, adata, sets: Optional[Dict[str, List[str]]] = None, min_genes: int = 4) -> pd.DataFrame:
    """Does a pathway sit systematically high or low in the anchor-correlation ranking?

    For each pathway: mean rho of its genes, median rank, and a two-sided Mann-Whitney test of the pathway genes' rho
    against all other ranked genes (a rank-based, set-level test; no assumption on the rho distribution).
    """
    sets = sets or PANEL
    ranks = pd.Series(np.arange(1, len(rho) + 1), index=rho.index)
    rows = []
    for name, genes in sets.items():
        r = resolve(adata, genes)
        vs = [r[g] for g in genes if r[g] is not None and r[g] in rho.index]
        if len(vs) < min_genes:
            continue
        inset = rho.loc[vs]
        rest = rho.drop(vs)
        p = _mwu_p(inset.values, rest.values)
        rows.append({"pathway": name, "n_genes": len(vs), "mean_rho": float(inset.mean()), "median_rank": float(ranks[vs].median()),
                     "n_ranked": len(rho), "frac_top_500": float((ranks[vs] <= 500).mean()), "MWU_p": p,
                     "top_gene": inset.idxmax(), "top_rho": float(inset.max())})
    return pd.DataFrame(rows).set_index("pathway").sort_values("mean_rho", ascending=False)


# --------------------------------------------------------------------------------------------------------------------
# spatial neighbourhoods (Xenium)
# --------------------------------------------------------------------------------------------------------------------

def neighbourhood_fraction(adata, sample_col: str, source_mask: np.ndarray, target_mask: np.ndarray, radius: float = 30.0,
                           min_source: int = 30) -> pd.DataFrame:
    """Per sample: among cells within ``radius`` of a source cell, what fraction are target cells?"""
    from sklearn.neighbors import KDTree
    rows = []
    coords_all = np.asarray(adata.obsm["spatial"], dtype=float)
    samples = adata.obs[sample_col].astype(str).values
    for s in np.unique(samples):
        m = samples == s
        xy, src, tgt = coords_all[m], source_mask[m], target_mask[m]
        if src.sum() < min_source or tgt.sum() < 10:
            continue
        tree = KDTree(xy)
        ind = tree.query_radius(xy[src], r=radius)
        src_idx = np.where(src)[0]
        n_nb = np.array([len(i) - 1 for i in ind])
        n_tg = np.array([tgt[i].sum() - tgt[j] for i, j in zip(ind, src_idx)])
        rows.append({"sample": s, "n_source": int(src.sum()), "n_neighbours": int(n_nb.sum()),
                     "frac_target_neighbours": n_tg.sum() / max(n_nb.sum(), 1), "baseline_frac_target": float(tgt.mean())})
    cols = ["sample", "n_source", "n_neighbours", "frac_target_neighbours", "baseline_frac_target"]
    return pd.DataFrame(rows, columns=cols).set_index("sample")


def compare_neighbourhoods(adata, sample_col: str, src_a: np.ndarray, src_b: np.ndarray, target_mask: np.ndarray,
                           label_a: str, label_b: str, radius: float = 30.0, min_source: int = 30, verbose: bool = True) -> pd.DataFrame:
    """Paired per-sample comparison of target-neighbour fractions around two source populations."""
    A = neighbourhood_fraction(adata, sample_col, src_a, target_mask, radius, min_source)
    B = neighbourhood_fraction(adata, sample_col, src_b, target_mask, radius, min_source)
    j = A[["frac_target_neighbours"]].join(B[["frac_target_neighbours"]], lsuffix=f"_{label_a}", rsuffix=f"_{label_b}", how="inner")
    j["ratio"] = j.iloc[:, 0] / j.iloc[:, 1].replace(0, np.nan)
    if verbose and len(j) >= 3:
        w = stats.wilcoxon(j.iloc[:, 0], j.iloc[:, 1])
        print(f"paired Wilcoxon across {len(j)} samples: p = {w.pvalue:.3g}; median ratio = {j['ratio'].median():.2f}")
    return j


def neighbourhood_mean(adata, sample_col: str, source_mask: np.ndarray, values: np.ndarray, radius: float = 30.0,
                       min_source: int = 30, exclude_self: bool = True) -> pd.DataFrame:
    """Per sample: mean of ``values`` over the neighbours (within ``radius``) of source cells, e.g. a pathway score."""
    from sklearn.neighbors import KDTree
    rows = []
    coords_all = np.asarray(adata.obsm["spatial"], dtype=float)
    samples = adata.obs[sample_col].astype(str).values
    values = np.asarray(values, dtype=float)
    for s in np.unique(samples):
        m = samples == s
        xy, src, val = coords_all[m], source_mask[m], values[m]
        if src.sum() < min_source:
            continue
        tree = KDTree(xy)
        ind = tree.query_radius(xy[src], r=radius)
        src_idx = np.where(src)[0]
        tot, n = 0.0, 0
        for i, j in zip(ind, src_idx):
            nb = i[i != j] if exclude_self else i
            v = val[nb]
            v = v[~np.isnan(v)]
            tot += v.sum(); n += len(v)
        rows.append({"sample": s, "n_source": int(src.sum()), "n_neighbours": n, "mean_neighbour_value": tot / max(n, 1)})
    return pd.DataFrame(rows, columns=["sample", "n_source", "n_neighbours", "mean_neighbour_value"]).set_index("sample")


# --------------------------------------------------------------------------------------------------------------------
# dataset registry and naming
# --------------------------------------------------------------------------------------------------------------------

SHORT: Dict[str, str] = {
    "Park2023_AD_hippocampus_scRNA": "Park 2023 AD hippocampus (m, sc)",
    "Aging_snRNA_HIP_CP_mouse": "Aging HIP+CP (m, sn)",
    "Ximerakis2019_aging_brain_scRNA": "Ximerakis 2019 aging brain (m, sc)",
    "Kaya2022_aged_WM_vs_GM_scRNA": "Kaya 2022 aged WM vs GM (m, sc)",
    "Zhou2020_5XFAD_snRNA": "Zhou 2020 5XFAD (m, sn)",
    "Chen2020_ST_AppNLGF_mouse": "Chen 2020 ST AD (m, spatial)",
    "LPC_Cuprizone_CC_snRNA_mouse": "LPC + cuprizone CC (m, sn)",
    "Serpina3n_Cuprizone_snRNA_mouse": "Serpina3n-cKO cuprizone (m, sn)",
    "Jakel2019_MS_human_snRNA": "Jäkel 2019 MS (h, sn)",
    "Absinta2021_MS_human_snRNA": "Absinta 2021 MS (h, sn)",
    "Schirmer2019_MS_human_snRNA": "Schirmer 2019 MS (h, sn)",
    "LermaMartin2024_MS_human_snRNA": "Lerma-Martin 2024 MS (h, sn)",
    "LermaMartin2024_MS_human_Visium": "Lerma-Martin 2024 MS (h, Visium)",
    "SenescentGlia2025_MS_human_Visium": "Senescent-glia 2025 MS (h, Visium)",
    "Leng2021_AD_human_snRNA": "Leng 2021 AD (h, sn)",
    "Sadick2022_AD_human_astro_oligo_snRNA": "Sadick 2022 AD (h, sn)",
}
# the integrated atlas written by OligoC4b is not a primary dataset and is skipped by processed_paths()
SKIP_FILES = {"oligodendrocyte_atlas"}


def short(name: str) -> str:
    return SHORT.get(name, name)


def processed_paths(processed_dir: Optional[str] = None) -> Dict[str, str]:
    """name -> path for every harmonised public .h5ad written by OligoC4b (atlas excluded)."""
    d = processed_dir or PROCESSED_DIR
    paths = {os.path.basename(p)[:-5]: p for p in sorted(glob.glob(os.path.join(d, "*.h5ad")))}
    return {k: v for k, v in paths.items() if k not in SKIP_FILES}


def resolve_path(env_var: str, fallback: str) -> Optional[str]:
    p = os.getenv(env_var) or fallback
    return p if p and os.path.exists(p) else None


# --------------------------------------------------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------------------------------------------------

def heat(df: pd.DataFrame, title: str, cbar_label: str, fmt: str = ".1f", cmap: str = "RdBu_r", center: Optional[float] = 0,
         vmin=None, vmax=None, annot=None, row_labels=None, col_labels=None, annot_size: int = 7, show: bool = True):
    """Annotated heatmap with readable labels; one hue for magnitudes ('Reds'), diverging for fold changes ('RdBu_r')."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    df = df.copy()
    df.index = row_labels if row_labels is not None else [short(str(i)) for i in df.index]
    if col_labels is not None:
        df.columns = col_labels
    h, w = df.shape
    fig, ax = plt.subplots(figsize=(min(0.55 * w + 4.2, 30), 0.42 * h + 1.9))
    sns.heatmap(df.astype(float), annot=True if annot is None else annot, fmt=fmt if annot is None else "", cmap=cmap, center=center,
                vmin=vmin, vmax=vmax, linewidths=0.5, linecolor="white", annot_kws={"size": annot_size},
                cbar_kws={"label": cbar_label, "shrink": 0.7}, ax=ax)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title(title, loc="left", fontsize=11, pad=10)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8); plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout()
    if show:
        plt.show()
    return ax


def pathway_colour_bar(genes: Sequence[str]):
    """Return a list of pathway names aligned with ``genes`` (for colouring gene axes by pathway)."""
    return [PATHWAY_OF.get(g, "context") for g in genes]
