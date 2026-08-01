
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Post-processing TensorQTL outputs for the following file layout:

cis/
  sv_eqtl.cis_qtl.txt.gz

cis_nominal/
  sv_eqtl.cis_qtl_pairs.W82_01.txt.gz
  ...
  sv_eqtl.cis_qtl_pairs.W82_20.txt.gz

trans/
  sv_eqtl.trans_qtl_pairs.txt.gz

Main tasks
----------
1) Read cis permutation results and define significant eGenes by q-value.
2) Parse all cis_nominal chromosome files and recover all significant cis variant-gene
   pairs using each gene's own pval_nominal_threshold from cis_qtl.
3) Summarize cis results:
   - significant eGenes
   - significant cis pairs
   - counts per gene / per variant
   - lead variant distribution
4) Summarize sparse trans results as candidate pairs only:
   - optional BH q-values within sparse output only (exploratory, not genome-wide)
   - hotspot counts per variant / per phenotype
5) Write a machine-readable summary file.

Notes
-----
- For cis pairs, this script follows the logic used by TensorQTL post.get_significant_pairs():
  keep genes with qval <= FDR, then keep nominal pairs with
  pval_nominal < phenotype-specific pval_nominal_threshold.
- trans_qtl_pairs from TensorQTL default shell mode are sparse output candidates
  (e.g. p < 1e-5 after cis filtering), so any BH correction done here is only on
  the sparse candidate set and must not be over-interpreted as full trans FDR.

Author: ChatGPT
"""

import argparse
import gzip
import json
import math
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


def eprint(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def mkdir(path: str):
    os.makedirs(path, exist_ok=True)


def natural_key(s: str):
    # robust sorting for W82_01 ... W82_20 and generic mixed names
    parts = re.split(r'(\d+)', os.path.basename(s))
    out = []
    for p in parts:
        if p.isdigit():
            out.append(int(p))
        else:
            out.append(p)
    return out


def detect_sep(path: str) -> str:
    # TensorQTL text outputs are tab-delimited
    return "\t"


def read_table(path: str, usecols=None, dtype=None, nrows=None):
    return pd.read_csv(
        path,
        sep=detect_sep(path),
        compression="infer",
        usecols=usecols,
        dtype=dtype,
        low_memory=False,
        nrows=nrows,
    )


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1, dtype=float)
    q = np.empty(n, dtype=float)
    q[order] = pvals[order] * n / ranks
    # enforce monotonicity from largest p to smallest p
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def safe_to_numeric(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def write_df(df: pd.DataFrame, path: str, index: bool = False):
    compression = "gzip" if path.endswith(".gz") else None
    df.to_csv(path, sep="\t", index=index, compression=compression)


def init_tsv_writer(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def stream_write_df(df: pd.DataFrame, handle, header: bool):
    df.to_csv(handle, sep="\t", index=False, header=header)


def collect_nominal_files(nominal_dir: str):
    files = []
    for fn in os.listdir(nominal_dir):
        if fn.endswith(".txt.gz") and ".cis_qtl_pairs." in fn:
            files.append(os.path.join(nominal_dir, fn))
    if not files:
        raise FileNotFoundError(f"No cis nominal files found in: {nominal_dir}")
    files.sort(key=natural_key)
    return files


def summarize_series_count(series: pd.Series, top_n=20, col_name="item", count_name="n"):
    vc = series.value_counts(dropna=False).rename_axis(col_name).reset_index(name=count_name)
    return vc.head(top_n), vc


def parse_args():
    ap = argparse.ArgumentParser(
        description="Post-process TensorQTL cis / cis_nominal / trans outputs."
    )
    ap.add_argument("--cis", required=True, help="Path to sv_eqtl.cis_qtl.txt.gz")
    ap.add_argument("--cis-nominal-dir", required=True, help="Directory containing sv_eqtl.cis_qtl_pairs.*.txt.gz")
    ap.add_argument("--trans", required=True, help="Path to sv_eqtl.trans_qtl_pairs.txt.gz")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for significant eGenes (default: 0.05)")
    ap.add_argument("--chunksize", type=int, default=500000, help="Chunk size when scanning cis_nominal files (default: 500000)")
    ap.add_argument("--trans-max-rows", type=int, default=0, help="Optional hard cap when reading trans sparse output; 0 means read all")
    return ap.parse_args()


def main():
    args = parse_args()
    mkdir(args.outdir)

    summary = {
        "inputs": {
            "cis": os.path.abspath(args.cis),
            "cis_nominal_dir": os.path.abspath(args.cis_nominal_dir),
            "trans": os.path.abspath(args.trans),
        },
        "parameters": {
            "fdr": args.fdr,
            "chunksize": args.chunksize,
            "trans_max_rows": args.trans_max_rows,
        },
    }

    # ------------------------------------------------------------------
    # 1) Read cis results and define significant eGenes
    # ------------------------------------------------------------------
    eprint("[1/4] Reading cis results:", args.cis)
    cis_df = read_table(args.cis)
    if "phenotype_id" not in cis_df.columns:
        raise ValueError("cis result must contain a 'phenotype_id' column.")
    if "qval" not in cis_df.columns:
        raise ValueError("cis result must contain a 'qval' column.")
    if "pval_nominal_threshold" not in cis_df.columns:
        raise ValueError("cis result must contain a 'pval_nominal_threshold' column.")
    if "variant_id" not in cis_df.columns:
        raise ValueError("cis result must contain a 'variant_id' column.")

    numeric_cols = [
        "qval", "pval_nominal_threshold", "pval_nominal", "pval_beta", "pval_perm",
        "slope", "slope_se", "tss_distance", "af", "ma_count", "ma_samples"
    ]
    cis_df = safe_to_numeric(cis_df, numeric_cols)

    write_df(cis_df, os.path.join(args.outdir, "cis.all.tsv.gz"))
    sig_cis = cis_df.loc[cis_df["qval"] <= args.fdr].copy()
    sig_cis = sig_cis.sort_values(["qval", "pval_beta" if "pval_beta" in sig_cis.columns else "pval_nominal"], na_position="last")
    write_df(sig_cis, os.path.join(args.outdir, f"cis.sig_eGenes.FDR{args.fdr}.tsv.gz"))

    lead_cols = [c for c in [
        "phenotype_id", "variant_id", "qval", "pval_nominal_threshold", "pval_nominal",
        "pval_beta", "pval_perm", "slope", "slope_se", "tss_distance", "af", "ma_count", "ma_samples"
    ] if c in sig_cis.columns]
    write_df(sig_cis[lead_cols], os.path.join(args.outdir, f"cis.sig_eGenes.lead_only.FDR{args.fdr}.tsv.gz"))

    summary["cis"] = {
        "n_all_genes": int(cis_df.shape[0]),
        "n_sig_eGenes": int(sig_cis.shape[0]),
    }

    # basic lead summaries
    if not sig_cis.empty:
        top_lead_var, all_lead_var = summarize_series_count(sig_cis["variant_id"], col_name="variant_id", count_name="n_eGenes")
        write_df(all_lead_var, os.path.join(args.outdir, f"cis.lead_variant_counts.FDR{args.fdr}.tsv.gz"))
        summary["cis"]["n_unique_lead_variants"] = int(sig_cis["variant_id"].nunique())

        if "tss_distance" in sig_cis.columns:
            td = pd.to_numeric(sig_cis["tss_distance"], errors="coerce").dropna()
            if len(td) > 0:
                summary["cis"]["lead_tss_distance_abs_median"] = float(np.median(np.abs(td)))
                summary["cis"]["lead_tss_distance_abs_mean"] = float(np.mean(np.abs(td)))

        if "slope" in sig_cis.columns:
            slope = pd.to_numeric(sig_cis["slope"], errors="coerce").dropna()
            if len(slope) > 0:
                summary["cis"]["lead_slope_positive"] = int((slope > 0).sum())
                summary["cis"]["lead_slope_negative"] = int((slope < 0).sum())

    # phenotype -> threshold map
    threshold_map = sig_cis.set_index("phenotype_id")["pval_nominal_threshold"].to_dict()
    phenotype_set = set(threshold_map.keys())

    # phenotype-level lookup columns to append back to significant nominal pairs
    append_cols = [c for c in [
        "phenotype_id", "variant_id", "qval", "pval_nominal_threshold", "pval_nominal",
        "pval_beta", "pval_perm", "slope", "slope_se", "tss_distance", "af", "ma_count", "ma_samples"
    ] if c in sig_cis.columns]
    sig_lookup = sig_cis[append_cols].copy()
    sig_lookup = sig_lookup.rename(columns={
        "variant_id": "lead_variant_id",
        "pval_nominal": "lead_pval_nominal",
        "pval_beta": "lead_pval_beta",
        "pval_perm": "lead_pval_perm",
        "slope": "lead_slope",
        "slope_se": "lead_slope_se",
        "tss_distance": "lead_tss_distance",
        "af": "lead_af",
        "ma_count": "lead_ma_count",
        "ma_samples": "lead_ma_samples",
    })

    # ------------------------------------------------------------------
    # 2) Scan all cis_nominal files and recover all significant cis pairs
    # ------------------------------------------------------------------
    eprint("[2/4] Scanning cis_nominal files and extracting all significant cis pairs")
    nominal_files = collect_nominal_files(args.cis_nominal_dir)
    out_pairs_path = os.path.join(args.outdir, f"cis.sig_pairs.FDR{args.fdr}.tsv.gz")
    out_gene_counts_path = os.path.join(args.outdir, f"cis.sig_pairs.by_gene.FDR{args.fdr}.tsv.gz")
    out_variant_counts_path = os.path.join(args.outdir, f"cis.sig_pairs.by_variant.FDR{args.fdr}.tsv.gz")

    gene_pair_counter = defaultdict(int)
    variant_pair_counter = defaultdict(int)
    chr_pair_counter = defaultdict(int)

    total_nominal_rows = 0
    total_filtered_by_gene = 0
    total_sig_pairs = 0
    first_write = True

    with init_tsv_writer(out_pairs_path) as out_handle:
        for i, path in enumerate(nominal_files, 1):
            eprint(f"  - [{i}/{len(nominal_files)}] {os.path.basename(path)}")
            reader = pd.read_csv(
                path,
                sep="\t",
                compression="infer",
                chunksize=args.chunksize,
                low_memory=False,
            )
            for chunk in reader:
                total_nominal_rows += chunk.shape[0]
                if "phenotype_id" not in chunk.columns or "pval_nominal" not in chunk.columns or "variant_id" not in chunk.columns:
                    raise ValueError(f"{path} must contain phenotype_id, variant_id, pval_nominal columns.")

                chunk = safe_to_numeric(chunk, ["pval_nominal", "slope", "slope_se", "af", "ma_count", "ma_samples"])
                chunk = chunk[chunk["phenotype_id"].isin(phenotype_set)].copy()
                if chunk.empty:
                    continue
                total_filtered_by_gene += chunk.shape[0]

                chunk["pval_nominal_threshold"] = chunk["phenotype_id"].map(threshold_map)
                chunk = chunk[chunk["pval_nominal"] < chunk["pval_nominal_threshold"]].copy()
                if chunk.empty:
                    continue

                # append lead / phenotype-level info
                chunk = chunk.merge(sig_lookup, on="phenotype_id", how="left", suffixes=("", "_from_cis"))
                # is this nominal pair itself the lead pair?
                chunk["is_lead_pair"] = (
                    (chunk["variant_id"].astype(str) == chunk["lead_variant_id"].astype(str))
                )

                # write out
                stream_write_df(chunk, out_handle, header=first_write)
                first_write = False

                total_sig_pairs += chunk.shape[0]
                for phen, cnt in chunk["phenotype_id"].value_counts().items():
                    gene_pair_counter[phen] += int(cnt)
                for var, cnt in chunk["variant_id"].value_counts().items():
                    variant_pair_counter[var] += int(cnt)
                chr_name = os.path.basename(path)
                chr_pair_counter[chr_name] += int(chunk.shape[0])

    by_gene = (
        pd.DataFrame({
            "phenotype_id": list(gene_pair_counter.keys()),
            "n_sig_pairs": list(gene_pair_counter.values())
        })
        .sort_values(["n_sig_pairs", "phenotype_id"], ascending=[False, True])
    )
    if not by_gene.empty:
        by_gene = by_gene.merge(
            sig_cis[[c for c in ["phenotype_id", "variant_id", "qval", "pval_nominal_threshold", "tss_distance", "slope"] if c in sig_cis.columns]]
            .rename(columns={"variant_id": "lead_variant_id", "tss_distance": "lead_tss_distance", "slope": "lead_slope"}),
            on="phenotype_id",
            how="left",
        )
    write_df(by_gene, out_gene_counts_path)

    by_variant = (
        pd.DataFrame({
            "variant_id": list(variant_pair_counter.keys()),
            "n_sig_pairs": list(variant_pair_counter.values())
        })
        .sort_values(["n_sig_pairs", "variant_id"], ascending=[False, True])
    )
    write_df(by_variant, out_variant_counts_path)

    chr_counts = (
        pd.DataFrame({
            "cis_nominal_file": list(chr_pair_counter.keys()),
            "n_sig_pairs": list(chr_pair_counter.values())
        })
        .sort_values("cis_nominal_file")
    )
    write_df(chr_counts, os.path.join(args.outdir, f"cis.sig_pairs.by_nominal_file.FDR{args.fdr}.tsv.gz"))

    summary["cis"].update({
        "n_nominal_files": int(len(nominal_files)),
        "n_nominal_rows_scanned": int(total_nominal_rows),
        "n_nominal_rows_in_sig_eGenes": int(total_filtered_by_gene),
        "n_sig_cis_pairs": int(total_sig_pairs),
        "n_sig_cis_variants": int(by_variant["variant_id"].nunique()) if not by_variant.empty else 0,
    })

    # ------------------------------------------------------------------
    # 3) Read trans sparse output and summarize candidate hotspots
    # ------------------------------------------------------------------
    eprint("[3/4] Reading trans sparse output:", args.trans)
    trans_read_nrows = None if args.trans_max_rows == 0 else args.trans_max_rows
    trans_df = read_table(args.trans, nrows=trans_read_nrows)
    required_trans_cols = {"phenotype_id", "variant_id"}
    if not required_trans_cols.issubset(set(trans_df.columns)):
        raise ValueError("trans result must contain at least phenotype_id and variant_id columns.")
    trans_df = safe_to_numeric(trans_df, [c for c in ["pval", "maf"] if c in trans_df.columns])

    # exploratory BH on sparse output only
    if "pval" in trans_df.columns:
        valid = trans_df["pval"].notna().values
        sparse_q = np.full(trans_df.shape[0], np.nan, dtype=float)
        if valid.sum() > 0:
            sparse_q[valid] = bh_fdr(trans_df.loc[valid, "pval"].values)
        trans_df["sparse_bh_qval"] = sparse_q

    write_df(trans_df, os.path.join(args.outdir, "trans.sparse_pairs.tsv.gz"))

    trans_by_variant = (
        trans_df.groupby("variant_id", dropna=False)
        .agg(
            n_pairs=("phenotype_id", "size"),
            n_unique_phenotypes=("phenotype_id", "nunique"),
            min_pval=("pval", "min") if "pval" in trans_df.columns else ("phenotype_id", "size"),
            min_sparse_bh_qval=("sparse_bh_qval", "min") if "sparse_bh_qval" in trans_df.columns else ("phenotype_id", "size"),
        )
        .reset_index()
        .sort_values(["n_unique_phenotypes", "n_pairs"], ascending=[False, False])
    )
    write_df(trans_by_variant, os.path.join(args.outdir, "trans.hotspot_by_variant.tsv.gz"))

    trans_by_gene = (
        trans_df.groupby("phenotype_id", dropna=False)
        .agg(
            n_pairs=("variant_id", "size"),
            n_unique_variants=("variant_id", "nunique"),
            min_pval=("pval", "min") if "pval" in trans_df.columns else ("variant_id", "size"),
            min_sparse_bh_qval=("sparse_bh_qval", "min") if "sparse_bh_qval" in trans_df.columns else ("variant_id", "size"),
        )
        .reset_index()
        .sort_values(["n_unique_variants", "n_pairs"], ascending=[False, False])
    )
    write_df(trans_by_gene, os.path.join(args.outdir, "trans.hotspot_by_gene.tsv.gz"))

    # overlap with cis eGenes
    if not sig_cis.empty:
        trans_by_gene = trans_by_gene.merge(
            sig_cis[["phenotype_id", "variant_id", "qval"]].rename(columns={
                "variant_id": "cis_lead_variant_id",
                "qval": "cis_qval"
            }),
            on="phenotype_id",
            how="left",
        )
        trans_by_gene["is_sig_cis_eGene"] = trans_by_gene["cis_qval"].notna()
        write_df(trans_by_gene, os.path.join(args.outdir, "trans.hotspot_by_gene.with_cis_overlap.tsv.gz"))

    summary["trans"] = {
        "n_sparse_pairs": int(trans_df.shape[0]),
        "n_variants": int(trans_df["variant_id"].nunique()),
        "n_phenotypes": int(trans_df["phenotype_id"].nunique()),
    }
    if "pval" in trans_df.columns and not trans_df["pval"].dropna().empty:
        summary["trans"]["min_pval"] = float(trans_df["pval"].min())
    if "sparse_bh_qval" in trans_df.columns and not trans_df["sparse_bh_qval"].dropna().empty:
        summary["trans"]["n_sparse_bh_qval_le_0_05"] = int((trans_df["sparse_bh_qval"] <= 0.05).sum())

    # ------------------------------------------------------------------
    # 4) Final summary
    # ------------------------------------------------------------------
    eprint("[4/4] Writing summary")
    summary["outputs"] = {
        "cis_all": "cis.all.tsv.gz",
        "cis_sig_eGenes": f"cis.sig_eGenes.FDR{args.fdr}.tsv.gz",
        "cis_sig_eGenes_lead_only": f"cis.sig_eGenes.lead_only.FDR{args.fdr}.tsv.gz",
        "cis_sig_pairs": f"cis.sig_pairs.FDR{args.fdr}.tsv.gz",
        "cis_sig_pairs_by_gene": f"cis.sig_pairs.by_gene.FDR{args.fdr}.tsv.gz",
        "cis_sig_pairs_by_variant": f"cis.sig_pairs.by_variant.FDR{args.fdr}.tsv.gz",
        "cis_sig_pairs_by_nominal_file": f"cis.sig_pairs.by_nominal_file.FDR{args.fdr}.tsv.gz",
        "cis_lead_variant_counts": f"cis.lead_variant_counts.FDR{args.fdr}.tsv.gz",
        "trans_sparse_pairs": "trans.sparse_pairs.tsv.gz",
        "trans_hotspot_by_variant": "trans.hotspot_by_variant.tsv.gz",
        "trans_hotspot_by_gene": "trans.hotspot_by_gene.tsv.gz",
        "trans_hotspot_by_gene_with_cis_overlap": "trans.hotspot_by_gene.with_cis_overlap.tsv.gz",
        "summary_json": "summary.json",
    }

    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # human-readable notes
    notes = [
        "1) cis.sig_eGenes.* : qval <= FDR 的显著 eGenes（来自 cis permutation 结果）",
        "2) cis.sig_pairs.*   : 对显著 eGenes，在各 cis_nominal 文件中按各自 pval_nominal_threshold 提取的全部显著 cis 配对",
        "3) trans.sparse_pairs.tsv.gz : TensorQTL trans 默认稀疏输出候选集；不是完整 trans 检验全集",
        "4) sparse_bh_qval    : 仅对稀疏 trans 候选集内部做 BH，不能等同于全基因组 trans FDR",
        "5) trans.hotspot_*   : 仅用于候选 hotspot 排序与后续人工审阅/验证",
    ]
    with open(os.path.join(args.outdir, "README.postprocess.txt"), "w", encoding="utf-8") as f:
        for line in notes:
            f.write(line + "\n")

    eprint("Done.")
    eprint(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
