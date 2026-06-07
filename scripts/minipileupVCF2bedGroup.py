#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse, gzip
from collections import defaultdict

def smart_open(path, mode='rt'):
    if path == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)

def load_groups(group_tsv):
    """
    读取分组表：第一列样本名、第二列组名（制表符分隔）
    返回 sample2grp 和 grp_order（按出现顺序）
    """
    sample2grp = {}
    grp_order, seen = [], set()
    with smart_open(group_tsv, 'rt') as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = ln.split('\t')
            if len(p) < 2:
                continue
            s, g = p[0].strip(), p[1].strip()
            sample2grp[s] = g
            if g not in seen:
                seen.add(g); grp_order.append(g)
    if not grp_order:
        sys.stderr.write("[warn] 分组文件里没有有效的组名；第6列将为空。\n")
    return sample2grp, grp_order

def infer_len_and_type(ref, alt):
    """
    根据 REF 与第一个 ALT 的长度推断长度与类型（INS/DEL）。
    多等位取第一个 ALT；符号等位(<DEL>)返回 (0, 'NA')。
    """
    if ',' in alt:
        alt = alt.split(',', 1)[0]
    if alt.startswith('<') and alt.endswith('>'):
        return 0, 'NA'
    lr, la = len(ref), len(alt)
    if la > lr:
        return la - lr, 'INS'
    elif lr > la:
        return lr - la, 'DEL'
    else:
        return 0, 'NA'

def is_hom_alt(gt_raw, any_hom=False):
    """
    判断基因型是否为纯合 ALT。
    - any_hom=False: 仅 1/1 记为纯合 ALT
    - any_hom=True : 任何非0的纯合（1/1, 2/2, 3/3, ...）都记为纯合 ALT
    兼容相位写法 1|1
    """
    if gt_raw in ('.', './.', '.|.'):
        return False
    gt = gt_raw.split(':', 1)[0].replace('|', '/')
    alleles = gt.split('/')
    # 排除缺失
    if any(a == '.' or a == '' for a in alleles):
        return False
    # 所有等位完全相同
    uniq = set(alleles)
    if len(uniq) != 1:
        return False
    a = next(iter(uniq))
    if any_hom:
        return a != '0'
    else:
        return a == '1'

def main():
    ap = argparse.ArgumentParser(
        description="VCF to TSV：CHROM, POS, LEN, TYPE, ALT sample list, ALT counts"
    )
    ap.add_argument("-i", "--vcf", required=True, help="input VCF（.vcf / .vcf.gz）")
    ap.add_argument("-g", "--groups", required=True, help="sample group（sample<TAB>group）")
    ap.add_argument("-o", "--out", default="-", help="输出 TSV（默认stdout）")
    ap.add_argument("--any-hom", action="store_true",
                    help="ALT（1/1, 2/2, 3/3, ...）")
    args = ap.parse_args()

    sample2grp, grp_order = load_groups(args.groups)

    with smart_open(args.vcf, 'rt') as fin, smart_open(args.out, 'wt') as fout:
        samples = []
        # 读表头
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip('\n').split('\t')
                samples = header[9:]
                group_col_name = "/".join(grp_order) if grp_order else "group_counts"
                out_header = ["CHROM", "POS", "LEN", "TYPE", "HOM_ALT_SAMPLES", group_col_name]
                fout.write("\t".join(out_header) + "\n")
                break

        if not samples:
            sys.stderr.write("[error] 未找到 #CHROM 表头。确认是标准 VCF？\n")
            sys.exit(1)

        for line in fin:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = parts[:8]

            vlen, vtype = infer_len_and_type(ref, alt)
            hom_alt_samples = []
            grp_counts = defaultdict(int)

            fmt = parts[8] if len(parts) > 8 else "GT"
            for i, sname in enumerate(samples):
                idx = 9 + i
                if idx >= len(parts):
                    break
                field = parts[idx]
                if is_hom_alt(field, any_hom=args.any_hom):
                    hom_alt_samples.append(sname)
                    g = sample2grp.get(sname)
                    if g is not None:
                        grp_counts[g] += 1

            group_counts_str = "/".join(str(grp_counts.get(g, 0)) for g in grp_order) if grp_order else ""
            hom_alt_str = ",".join(hom_alt_samples)
            row = [chrom, pos, str(vlen), vtype, hom_alt_str, group_counts_str]
            fout.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
