#!/usr/bin/env python3
"""Calculate length-weighted average coverage across main chromosomes from samtools coverage output."""

import csv
import glob
import os
import sys

MAIN_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def calc_weighted_avg_coverage(filepath: str) -> float:
    total_len = 0
    weighted_sum = 0.0

    with open(filepath, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # skip comment header(s)

        # Handle possible multiple comment lines
        while header[0].startswith("#"):
            header = next(reader)

        for row in reader:
            rname = row[0]
            if rname not in MAIN_CHROMOSOMES:
                continue

            endpos = int(row[2])
            startpos = int(row[1])
            chrom_len = endpos - startpos + 1
            meandepth = float(row[6])

            weighted_sum += meandepth * chrom_len
            total_len += chrom_len

    return weighted_sum / total_len if total_len else 0.0


def parse_read_count(label: str) -> int:
    """Parse a label like '10m' or '500k' into an integer read count."""
    label = label.lower().strip()
    if label.endswith("m"):
        return int(float(label[:-1]) * 1_000_000)
    elif label.endswith("k"):
        return int(float(label[:-1]) * 1_000)
    else:
        return int(label)


def main():
    coverage_dir = os.path.dirname(os.path.abspath(__file__))
    files = sorted(glob.glob(os.path.join(coverage_dir, "coverage_*.tsv")))

    results = []
    for filepath in files:
        basename = os.path.basename(filepath)
        label = basename.replace("coverage_", "").replace(".tsv", "")
        read_count = parse_read_count(label)
        avg_depth = calc_weighted_avg_coverage(filepath)
        results.append((read_count, label.upper(), avg_depth))

    # Sort by decreasing total read count
    results.sort(key=lambda x: x[0], reverse=True)

    # Write TSV output
    outpath = os.path.join(coverage_dir, "avg_coverage.tsv")
    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Reads", "Weighted_Avg_Depth"])
        for read_count, label, avg_depth in results:
            writer.writerow([label, f"{avg_depth:.4f}"])

    print(f"Results written to {outpath}")
    print()
    with open(outpath) as f:
        print(f.read())


if __name__ == "__main__":
    main()
