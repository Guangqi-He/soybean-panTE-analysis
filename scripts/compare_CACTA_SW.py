#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
import math
import pandas as pd
from scipy.stats import wilcoxon


def broad_class(te_class):
    x = str(te_class).upper()

    if "CACTA" in x or x == "DNA/DTC":
        return "CACTA"
    if "GYPSY" in x:
        return "LTR/Gypsy"
    if "COPIA" in x:
        return "LTR/Copia"
    if x.startswith("LTR"):
        return "Other LTR"
    if x.startswith("DNA") or x.startswith("TIR"):
        return "Other DNA transposon"
    if x.startswith("LINE"):
        return "LINE"
    if x.startswith("SINE"):
        return "SINE"
    if "HELITRON" in x or x.startswith("RC/"):
        return "Helitron"
    if x in {"UNKNOWN", "UNCLASSIFIED"}:
        return "Unknown"

    return "Other TE"


def run_intersect(region_bed, annotation_bed):

    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        suffix=".tsv"
    )
    tmp.close()

    cmd = [
        "bedtools", "intersect",
        "-a", region_bed,
        "-b", annotation_bed,
        "-wao"
    ]

    with open(tmp.name, "w") as out:
        subprocess.run(
            cmd,
            stdout=out,
            check=True
        )

    return tmp.name


def read_intersection(
    filename,
    min_overlap_bp=80,
    min_region_cov=0.5,
    non_cacta_only=False
):

    result = {}

    with open(filename) as f:

        for line in f:

            x = line.rstrip().split("\t")

            # A BED: chr start end
            chrom = x[0]
            start = int(x[1])
            end = int(x[2])

            region_id = f"{chrom}:{start}-{end}"
            region_len = end - start

            if region_id not in result:
                result[region_id] = {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "region_len": region_len,
                    "hits": []
                }

            # No annotation
            if x[3] == ".":
                continue

            # B BED is expected to have 10 columns
            bchrom = x[3]
            bstart = int(x[4])
            bend = int(x[5])
            repeat = x[6]
            te_class = x[7]
            score = float(x[8])
            strand = x[9]

            overlap_bp = int(x[-1])

            hit_len = bend - bstart

            if overlap_bp <= 0:
                continue

            region_cov = overlap_bp / region_len
            hit_cov = overlap_bp / hit_len

            if overlap_bp < min_overlap_bp:
                continue

            if region_cov < min_region_cov:
                continue

            broad = broad_class(te_class)

            if non_cacta_only and broad == "CACTA":
                continue

            score_per_bp = score / hit_len

            result[region_id]["hits"].append({
                "repeat": repeat,
                "class": te_class,
                "broad_class": broad,
                "score": score,
                "bstart": bstart,
                "bend": bend,
                "hit_len": hit_len,
                "overlap_bp": overlap_bp,
                "region_cov": region_cov,
                "hit_cov": hit_cov,
                "score_per_bp": score_per_bp,
                "strand": strand
            })

    return result


def choose_best(region_dict):

    best = {}

    for region_id, info in region_dict.items():

        hits = info["hits"]

        if len(hits) == 0:
            best[region_id] = None
            continue

        # Select highest SW score among sufficiently overlapping hits
        hit = max(
            hits,
            key=lambda z: z["score"]
        )

        best[region_id] = hit

    return best


parser = argparse.ArgumentParser()

parser.add_argument(
    "--regions",
    required=True
)

parser.add_argument(
    "--cacta-bed",
    required=True
)

parser.add_argument(
    "--full-bed",
    required=True
)

parser.add_argument(
    "--label",
    required=True
)

parser.add_argument(
    "--out",
    required=True
)

parser.add_argument(
    "--min-overlap",
    type=int,
    default=80
)

parser.add_argument(
    "--min-region-cov",
    type=float,
    default=0.5
)

args = parser.parse_args()


# --------------------------------------
# BEDTools intersections
# --------------------------------------

cacta_intersect = run_intersect(
    args.regions,
    args.cacta_bed
)

full_intersect = run_intersect(
    args.regions,
    args.full_bed
)


# --------------------------------------
# Read annotations
# --------------------------------------

cacta_regions = read_intersection(
    cacta_intersect,
    min_overlap_bp=args.min_overlap,
    min_region_cov=args.min_region_cov,
    non_cacta_only=False
)

full_regions = read_intersection(
    full_intersect,
    min_overlap_bp=args.min_overlap,
    min_region_cov=args.min_region_cov,
    non_cacta_only=True
)


cacta_best = choose_best(cacta_regions)
full_best = choose_best(full_regions)


# --------------------------------------
# Construct paired table
# --------------------------------------

rows = []

all_regions = sorted(
    set(cacta_regions.keys()) |
    set(full_regions.keys())
)

for rid in all_regions:

    basic = (
        cacta_regions.get(rid)
        or full_regions.get(rid)
    )

    c = cacta_best.get(rid)
    f = full_best.get(rid)

    row = {
        "region_set": args.label,
        "region_id": rid,
        "chr": basic["chrom"],
        "start": basic["start"],
        "end": basic["end"],
        "region_length": basic["region_len"]
    }

    if c is not None:

        row.update({
            "CACTA_repeat": c["repeat"],
            "CACTA_class": c["class"],
            "CACTA_SW": c["score"],
            "CACTA_overlap_bp": c["overlap_bp"],
            "CACTA_region_cov": c["region_cov"],
            "CACTA_hit_cov": c["hit_cov"],
            "CACTA_score_per_bp": c["score_per_bp"]
        })

    else:

        row.update({
            "CACTA_repeat": None,
            "CACTA_class": None,
            "CACTA_SW": None,
            "CACTA_overlap_bp": None,
            "CACTA_region_cov": None,
            "CACTA_hit_cov": None,
            "CACTA_score_per_bp": None
        })

    if f is not None:

        row.update({
            "Full_repeat": f["repeat"],
            "Full_class": f["class"],
            "Full_broad_class": f["broad_class"],
            "Full_SW": f["score"],
            "Full_overlap_bp": f["overlap_bp"],
            "Full_region_cov": f["region_cov"],
            "Full_hit_cov": f["hit_cov"],
            "Full_score_per_bp": f["score_per_bp"]
        })

    else:

        row.update({
            "Full_repeat": None,
            "Full_class": None,
            "Full_broad_class": None,
            "Full_SW": None,
            "Full_overlap_bp": None,
            "Full_region_cov": None,
            "Full_hit_cov": None,
            "Full_score_per_bp": None
        })

    if c is not None and f is not None:

        row["Delta_SW"] = f["score"] - c["score"]

        row["SW_ratio_full_vs_CACTA"] = (
            f["score"] / c["score"]
            if c["score"] > 0
            else None
        )

        row["log2_SW_ratio"] = (
            math.log2(f["score"] / c["score"])
            if c["score"] > 0 and f["score"] > 0
            else None
        )

        row["Delta_score_per_bp"] = (
            f["score_per_bp"] -
            c["score_per_bp"]
        )

        row["Full_SW_higher"] = (
            "Yes"
            if f["score"] > c["score"]
            else "No"
        )

    else:

        row["Delta_SW"] = None
        row["SW_ratio_full_vs_CACTA"] = None
        row["log2_SW_ratio"] = None
        row["Delta_score_per_bp"] = None
        row["Full_SW_higher"] = None

    rows.append(row)


df = pd.DataFrame(rows)

df.to_csv(
    args.out + ".paired.tsv",
    sep="\t",
    index=False
)


# --------------------------------------
# Statistical summary
# --------------------------------------

paired = df.dropna(
    subset=["CACTA_SW", "Full_SW"]
).copy()


summary = []

summary.append(
    ["region_set", args.label]
)

summary.append(
    ["total_regions", len(df)]
)

summary.append(
    ["paired_regions", len(paired)]
)


if len(paired) > 0:

    summary.append([
        "median_CACTA_SW",
        paired["CACTA_SW"].median()
    ])

    summary.append([
        "median_full_nonCACTA_SW",
        paired["Full_SW"].median()
    ])

    summary.append([
        "median_Delta_SW",
        paired["Delta_SW"].median()
    ])

    summary.append([
        "fraction_full_SW_higher",
        (paired["Delta_SW"] > 0).mean()
    ])

    summary.append([
        "median_log2_SW_ratio",
        paired["log2_SW_ratio"].median()
    ])

    try:

        stat, p = wilcoxon(
            paired["Full_SW"],
            paired["CACTA_SW"],
            alternative="two-sided"
        )

        summary.append([
            "paired_Wilcoxon_p",
            p
        ])

    except Exception:

        summary.append([
            "paired_Wilcoxon_p",
            "NA"
        ])


pd.DataFrame(
    summary,
    columns=["metric", "value"]
).to_csv(
    args.out + ".summary.tsv",
    sep="\t",
    index=False
)

print(
    f"Output: {args.out}.paired.tsv"
)

print(
    f"Output: {args.out}.summary.tsv"
)