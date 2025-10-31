#!/usr/bin/env python3
"""
Discordance analysis for five variant callers (SNP-only lists).
- Defines MISSED: called by all other tools but NOT by this tool.
- Defines UNIQUE: called ONLY by this tool.
- Works for Normal or Abnormal datasets (file suffixes differ).

Input files (produced by venn_pipeline.sh):
  Normal:   <Tool>_Normal_SNPs.txt
  Abnormal: <Tool>_Abnorm_SNPs.txt

Outputs:
  - discordance_summary_<dataset>.csv
  - missed_<Tool>_<dataset>.txt
  - unique_<Tool>_<dataset>.txt
"""

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, Set, List

TOOLS = ["DeepVariant", "GATK", "FreeBayes", "VarScan", "BCFTools"]

def load_snp_set(path: Path) -> Set[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    with path.open() as fh:
        return {line.strip() for line in fh if line.strip()}

def build_file_map(base_dir: Path, dataset: str) -> Dict[str, Path]:
    """
    dataset ∈ {"Normal","Abnorm"}
    Returns mapping: tool -> SNP file Path
    """
    suffix = "Normal_SNPs.txt" if dataset == "Normal" else "Abnorm_SNPs.txt"
    fmap = {tool: base_dir / f"{tool}_{suffix}" for tool in TOOLS}
    return fmap

def compute_discordance(snp_sets: Dict[str, Set[str]]) -> Dict[str, Dict[str, int]]:
    """
    For each tool:
      missed = (intersection of all OTHER sets) - this set
      unique = this set - (union of all OTHER sets)
    Returns dict: tool -> {"missed": int, "unique": int, "total_called": int}
    """
    results: Dict[str, Dict[str, int]] = {}
    # Precompute union and intersections for speed
    for tool in TOOLS:
        others: List[str] = [t for t in TOOLS if t != tool]
        # intersection of others
        inter_others: Set[str] = set.intersection(*(snp_sets[t] for t in others))
        # union of others
        union_others: Set[str] = set.union(*(snp_sets[t] for t in others))

        missed = inter_others - snp_sets[tool]
        unique = snp_sets[tool] - union_others

        results[tool] = {
            "missed": len(missed),
            "unique": len(unique),
            "total_called": len(snp_sets[tool]),
        }
    return results

def save_lists(base_dir: Path, dataset: str, snp_sets: Dict[str, Set[str]]) -> None:
    """
    Save per-tool missed/unique SNP coordinate lists for downstream inspection.
    """
    for tool in TOOLS:
        others = [t for t in TOOLS if t != tool]
        inter_others = set.intersection(*(snp_sets[t] for t in others))
        union_others = set.union(*(snp_sets[t] for t in others))

        missed = sorted(inter_others - snp_sets[tool])
        unique = sorted(snp_sets[tool] - union_others)

        (base_dir / f"missed_{tool}_{dataset}.txt").write_text("\n".join(missed) + ("\n" if missed else ""))
        (base_dir / f"unique_{tool}_{dataset}.txt").write_text("\n".join(unique) + ("\n" if unique else ""))

def save_csv(base_dir: Path, dataset: str, summary: Dict[str, Dict[str, int]]) -> None:
    out_csv = base_dir / f"discordance_summary_{dataset}.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Tool", "Total_Called", "Missed", "Unique"])
        for tool in TOOLS:
            row = summary[tool]
            writer.writerow([tool, row["total_called"], row["missed"], row["unique"]])
    print(f"[+] Saved: {out_csv}")

def main():
    parser = argparse.ArgumentParser(
        description="Compute missed/unique SNP discordance across five callers."
    )
    parser.add_argument(
        "--dataset",
        choices=["Normal", "Abnorm"],
        required=True,
        help="Choose which dataset to analyze."
    )
    parser.add_argument(
        "--base-dir",
        default="/media/fahad/DNA_Work2/concordance/ForvanDiagram",
        help="Base directory containing SNP coordinate files. "
             "For Abnorm, script automatically appends '/Abnormal'."
    )
    args = parser.parse_args()

    base = Path(args.base_dir).expanduser().resolve()
    if args.dataset == "Abnorm":
        base = base / "Abnormal"

    print(f"Dataset: {args.dataset}")
    print(f"Working directory: {base}")

    fmap = build_file_map(base, args.dataset)
    for t, p in fmap.items():
        if not p.exists():
            raise SystemExit(f"ERROR: Required file not found: {p}")

    # Load
    snp_sets = {tool: load_snp_set(path) for tool, path in fmap.items()}

    # Compute
    summary = compute_discordance(snp_sets)

    # Print
    print("\n=== Missed SNPs (missed by one tool, called by ALL others) ===")
    for tool in TOOLS:
        print(f"{tool:10} : {summary[tool]['missed']}")

    print("\n=== Unique SNPs (called ONLY by this tool) ===")
    for tool in TOOLS:
        print(f"{tool:10} : {summary[tool]['unique']}")

    print("\n=== Totals (SNPs called by each tool) ===")
    for tool in TOOLS:
        print(f"{tool:10} : {summary[tool]['total_called']}")

    # Save artifacts
    save_lists(base, args.dataset, snp_sets)
    save_csv(base, args.dataset, summary)

if __name__ == "__main__":
    main()
