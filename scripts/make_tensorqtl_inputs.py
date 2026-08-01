#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd


def read_vcf_samples(path):
    with open(path) as f:
        samples = [line.strip() for line in f if line.strip()]
    if not samples:
        raise ValueError("VCF sample list is empty.")
    return samples


def read_gene_bed(path):
    # 假定为 4 列：chr start end gene
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError("Gene BED must have at least 4 columns: chr start end gene")
    df = df.iloc[:, :4].copy()
    df.columns = ["chr", "start", "end", "gene"]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def read_expr(path):
    df = pd.read_csv(path, sep="\t", header=0)
    if df.shape[1] < 2:
        raise ValueError("Expression file must have >=2 columns: gene + samples")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "gene"})
    return df


def read_cov_auto(path, vcf_samples):
    """
    自动判断 cov.txt 是：
    1) covariates x samples   （TensorQTL 推荐/CLI 直接需要）
    2) samples x covariates   （常见原始格式）
    返回统一后的 covariates x samples
    """
    cov = pd.read_csv(path, sep="\t", header=0, index_col=0)

    vcf_set = set(vcf_samples)
    col_overlap = len(vcf_set.intersection(set(cov.columns)))
    idx_overlap = len(vcf_set.intersection(set(cov.index.astype(str))))

    if col_overlap > 0 and idx_overlap == 0:
        # 已经是 covariates x samples
        cov_cs = cov.copy()
    elif idx_overlap > 0 and col_overlap == 0:
        # 是 samples x covariates，需要转置
        cov_cs = cov.T
    elif col_overlap > 0 and idx_overlap > 0:
        # 两边都像，优先认为列是样本
        cov_cs = cov.copy()
    else:
        raise ValueError(
            "Cannot detect covariate orientation. "
            "No sample IDs from VCF were found in either rows or columns of cov.txt"
        )

    return cov_cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-bed", required=True, help="chr start end gene")
    ap.add_argument("--expr", required=True, help="gene expression matrix: gene sample1 sample2 ...")
    ap.add_argument("--cov", required=True, help="covariates file")
    ap.add_argument("--vcf-samples", required=True, help="sample list from bcftools query -l")
    ap.add_argument("--out-bed", required=True, help="output phenotype BED(.gz)")
    ap.add_argument("--out-cov", required=True, help="output covariates txt")
    args = ap.parse_args()

    vcf_samples = read_vcf_samples(args.vcf_samples)
    pos = read_gene_bed(args.gene_bed)
    expr = read_expr(args.expr)
    cov_cs = read_cov_auto(args.cov, vcf_samples)  # covariates x samples

    # 检查重复 gene
    if pos["gene"].duplicated().any():
        dup = pos.loc[pos["gene"].duplicated(), "gene"].head(10).tolist()
        raise ValueError(f"Duplicated gene IDs in gene BED, e.g. {dup}")
    if expr["gene"].duplicated().any():
        dup = expr.loc[expr["gene"].duplicated(), "gene"].head(10).tolist()
        raise ValueError(f"Duplicated gene IDs in expression file, e.g. {dup}")

    expr_samples = list(expr.columns[1:])
    cov_samples = list(cov_cs.columns)

    common_samples = [s for s in vcf_samples if s in expr_samples and s in cov_samples]
    if len(common_samples) == 0:
        raise ValueError("No overlapping samples among VCF, expression and covariates.")

    print(f"Common samples kept: {len(common_samples)}")

    # 统一样本顺序
    expr = expr[["gene"] + common_samples].copy()
    cov_cs = cov_cs[common_samples].copy()

    # 按 gene 合并位置信息与表达矩阵
    merged = pos.merge(expr, on="gene", how="inner")
    print(f"Genes in BED: {pos.shape[0]}")
    print(f"Genes in expression: {expr.shape[0]}")
    print(f"Genes kept after merge: {merged.shape[0]}")

    if merged.shape[0] == 0:
        raise ValueError("No overlapping gene IDs between gene BED and expression matrix.")

    # TensorQTL phenotype BED
    merged = merged.rename(columns={"chr": "#chr", "gene": "phenotype_id"})
    out_cols = ["#chr", "start", "end", "phenotype_id"] + common_samples
    merged = merged[out_cols].copy()

    # 可选：按染色体和起点排序
    merged = merged.sort_values(by=["#chr", "start", "end", "phenotype_id"])

    # 写出 phenotype BED
    compression = "gzip" if args.out_bed.endswith(".gz") else None
    merged.to_csv(args.out_bed, sep="\t", index=False, compression=compression)

    # 写出 covariates（covariates x samples）
    cov_cs.to_csv(args.out_cov, sep="\t")

    print(f"Wrote phenotype BED: {args.out_bed}")
    print(f"Wrote covariates: {args.out_cov}")


if __name__ == "__main__":
    main()