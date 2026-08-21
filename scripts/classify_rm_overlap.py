#!/usr/bin/env python3

import argparse
from collections import defaultdict, Counter


def broad_class(te_class):
    if te_class is None:
        return "Unassigned"

    x = te_class.upper()

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

    if x in ("UNKNOWN", "UNCLASSIFIED"):
        return "Unknown"

    return "Other TE"


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True,
                    help="bedtools intersect -wao output")
parser.add_argument("--genome-size", type=int, required=True)
parser.add_argument("--prefix", required=True)

args = parser.parse_args()


# key = (chr, start, end)
hits = defaultdict(list)

with open(args.input) as f:
    for line in f:

        if not line.strip():
            continue

        x = line.rstrip().split("\t")

        # A BED = 3 columns
        achr = x[0]
        astart = int(x[1])
        aend = int(x[2])

        key = (achr, astart, aend)

        # No overlap
        if x[3] == ".":
            hits[key]
            continue

        # B BED starts from column 4
        bstart = int(x[4])
        bend = int(x[5])

        repeat_name = x[6]
        te_class = x[7]
        score = float(x[8])

        s = max(astart, bstart)
        e = min(aend, bend)

        if e > s:
            hits[key].append(
                (s, e, repeat_name, te_class, score)
            )


exact_bp = Counter()
broad_bp = Counter()

total_bp = 0


for (chrom, astart, aend), annotations in hits.items():

    total_bp += aend - astart

    boundaries = {astart, aend}

    for s, e, repeat_name, te_class, score in annotations:
        boundaries.add(s)
        boundaries.add(e)

    boundaries = sorted(boundaries)

    for i in range(len(boundaries) - 1):

        s = boundaries[i]
        e = boundaries[i + 1]

        if e <= s:
            continue

        candidates = []

        for hs, he, repeat_name, te_class, score in annotations:

            if hs <= s and he >= e:
                candidates.append(
                    (score, repeat_name, te_class)
                )

        length = e - s

        if not candidates:

            exact_bp["Unassigned"] += length
            broad_bp["Unassigned"] += length

        else:

            # retain best RepeatMasker SW score
            candidates.sort(reverse=True)

            score, repeat_name, te_class = candidates[0]

            exact_bp[te_class] += length
            broad_bp[broad_class(te_class)] += length


# exact classification
with open(args.prefix + ".exact.tsv", "w") as out:

    out.write(
        "TE_class\tbp\tfraction_of_region\tgenome_percent\n"
    )

    for cls, bp in exact_bp.most_common():

        out.write(
            f"{cls}\t"
            f"{bp}\t"
            f"{bp/total_bp:.6f}\t"
            f"{bp/args.genome_size*100:.6f}\n"
        )


# broad classification
with open(args.prefix + ".broad.tsv", "w") as out:

    out.write(
        "category\tbp\tfraction_of_region\tgenome_percent\n"
    )

    for cls, bp in broad_bp.most_common():

        out.write(
            f"{cls}\t"
            f"{bp}\t"
            f"{bp/total_bp:.6f}\t"
            f"{bp/args.genome_size*100:.6f}\n"
        )


print("Region bp:", total_bp)
print("Genome percent:",
      total_bp / args.genome_size * 100)