#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    sys.exit(f"Usage: python {sys.argv[0]} RepeatMasker.out")

with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        x = line.split()

        # skip header
        if len(x) < 15:
            continue

        try:
            score = float(x[0])
            qstart = int(x[5])
            qend   = int(x[6])
        except ValueError:
            continue

        chrom = x[4]
        strand = x[8]
        repeat = x[9]
        teclass = x[10]
        rm_id = x[14]

        # BED: 0-based, half-open
        bed_start = qstart - 1
        bed_end = qend

        # coordinates on repeat consensus
        if strand == "C":
            repeat_start = int(x[13])
            repeat_end   = int(x[12])
        else:
            repeat_start = int(x[11])
            repeat_end   = int(x[12])

        if repeat_start > repeat_end:
            repeat_start, repeat_end = repeat_end, repeat_start

        print(
            chrom,
            bed_start,
            bed_end,
            repeat,
            teclass,
            score,
            strand,
            repeat_start,
            repeat_end,
            rm_id,
            sep="\t"
        )