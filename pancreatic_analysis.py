#!/usr/bin/env python3
"""Compare selected gene expression with recorded cancer stage in TCGA-PAAD.

This is a research analysis, not a diagnostic tool. It uses primary GDC
RNA-seq TPM values and recorded clinical stage labels when available.
"""

import csv
import json
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
METADATA_FILE = BASE_DIR / "file_metadata.json"
CLINICAL_FILE = BASE_DIR / "clinical.tsv"
OUTPUT_DIR = BASE_DIR / "results"
TARGET_GENES = ["ANGPTL4", "S100P", "KRT19", "COL6A3", "PLAU", "MMP7"]


def load_clinical_stages(path):
    """Load clinical data and choose one usable recorded stage per patient."""
    clinical = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    patient_col = "cases.submitter_id"
    stage_col = "diagnoses.ajcc_pathologic_stage"
    if patient_col not in clinical.columns or stage_col not in clinical.columns:
        raise ValueError(
            "Expected columns were not found in clinical.tsv. "
            "Check the column names in that file."
        )

    stages_by_patient = defaultdict(list)
    for _, row in clinical.iterrows():
        patient = row[patient_col].strip()
        stage = row[stage_col].strip()
        if not patient or stage in ("", "--", "'--", "NA", "NaN"):
            continue
        stages_by_patient[patient].append(stage)

    patient_stages = {}
    conflicting = set()
    for patient, stages in stages_by_patient.items():
        unique_stages = sorted(set(stages))
        if len(unique_stages) == 1:
            patient_stages[patient] = unique_stages[0]
        else:
            conflicting.add(patient)
    return patient_stages, conflicting


def stage_group(stage):
    """Map recorded stages into early or late groups using a regular expression."""
    match = re.search(r"Stage\s+(IV|III|II|I|0)", stage, re.IGNORECASE)
    if not match:
        return None
    roman = match.group(1).upper()
    if roman in ("0", "I"):
        return "Early-stage"
    if roman in ("II", "III", "IV"):
        return "Late-stage"
    return None


def read_gene_tpm(path, target_genes):
    """Read TPM values for the selected genes from one GDC STAR-counts file."""
    values = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[0] == "gene_id":
                continue
            if len(fields) < 7:
                continue
            gene_name = fields[1]
            if gene_name in target_genes:
                try:
                    values[gene_name] = float(fields[6])
                except (ValueError, IndexError):
                    pass
    return values


def main():
    """Load real TCGA-PAAD data, match stages, summarize expression, and plot."""
    if not METADATA_FILE.exists() or not CLINICAL_FILE.exists():
        raise FileNotFoundError(
            "Required files not found. Keep file_metadata.json, clinical.tsv, "
            "and this script in the BINF2111 folder."
        )

    with open(METADATA_FILE, "r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    clinical_stages, conflicting_patients = load_clinical_stages(CLINICAL_FILE)
    OUTPUT_DIR.mkdir(exist_ok=True)

    sample_rows = []
    primary_file_count = 0
    no_patient_count = 0
    no_stage_count = 0
    conflict_count = 0

    for file_record in metadata:
        filepath = BASE_DIR / file_record.get("file_name", "")
        if not filepath.exists():
            candidates = list(BASE_DIR.glob("*/" + file_record.get("file_name", "")))
            if candidates:
                filepath = candidates[0]
            else:
                continue
        cases = file_record.get("cases", [])
        if not cases:
            no_patient_count += 1
            continue

        for case in cases:
            patient_id = case.get("submitter_id", "")
            samples = case.get("samples", [])
            sample_types = [s.get("sample_type", "") for s in samples]
            if "Primary Tumor" not in sample_types:
                continue
            primary_file_count += 1

            if not patient_id:
                no_patient_count += 1
                continue
            if patient_id in conflicting_patients:
                conflict_count += 1
                continue
            stage = clinical_stages.get(patient_id)
            group = stage_group(stage) if stage else None
            if group is None:
                no_stage_count += 1
                continue

            expression = read_gene_tpm(filepath, TARGET_GENES)
            row = {
                "patient_id": patient_id,
                "sample_id": samples[0].get("submitter_id", patient_id),
                "recorded_stage": stage,
                "stage_group": group,
            }
            row.update(expression)
            sample_rows.append(row)

    if not sample_rows:
        print("No samples could be matched to unambiguous recorded stages.")
        print("Check clinical.tsv column names and file paths before running.")
        print("Files in metadata:", len(metadata))
        print("Patients with conflicting stages:", len(conflicting_patients))
        return

    expression_df = pd.DataFrame(sample_rows)
    for gene in TARGET_GENES:
        if gene not in expression_df.columns:
            expression_df[gene] = float("nan")

    expression_df.to_csv(OUTPUT_DIR / "matched_expression_by_stage.csv", index=False)

    summary_rows = []
    for gene in TARGET_GENES:
        for group in ("Early-stage", "Late-stage"):
            values = expression_df.loc[expression_df["stage_group"] == group, gene]
            summary_rows.append({
                "gene": gene,
                "stage_group": group,
                "sample_count": int(values.count()),
                "mean_TPM": float(values.mean()) if len(values) else float("nan"),
                "median_TPM": float(values.median()) if len(values) else float("nan")
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTPUT_DIR / "gene_expression_summary.csv", index=False)

    print("=" * 58)
    print("Pancreatic Cancer Gene Expression Analysis")
    print("Dataset: TCGA-PAAD (real GDC RNA-seq files)")
    print("=" * 58)
    print("Expression files in metadata:", len(metadata))
    print("Primary-tumor file records checked:", primary_file_count)
    print("Unique patients with one unambiguous recorded stage:",
          expression_df["patient_id"].nunique())
    print("Matched samples analyzed:", len(expression_df))
    print("Excluded for missing/unusable stage:", no_stage_count)
    print("Excluded for conflicting clinical stages:", conflict_count)
    print("Excluded because no patient record was available:", no_patient_count)
    print("\nGene expression comparison (mean TPM):")
    print(summary_df.pivot(index="gene", columns="stage_group", values="mean_TPM"))
    print("\nMean and median values saved to:", OUTPUT_DIR / "gene_expression_summary.csv")
    print("Matched sample table saved to:", OUTPUT_DIR / "matched_expression_by_stage.csv")

    plot_df = summary_df.pivot(index="gene", columns="stage_group", values="mean_TPM")
    plot_df.plot(kind="bar", figsize=(11, 6))
    plt.ylabel("Mean TPM (unstranded)")
    plt.xlabel("Selected gene")
    plt.title("TCGA-PAAD selected gene expression by recorded stage group")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plot_path = OUTPUT_DIR / "gene_expression_by_stage.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print("Chart saved to:", plot_path)
    print("\nResearch note: these results show associations in this cohort;")
    print("they are not a validated diagnostic or treatment-response test.")


if __name__ == "__main__":
    main()
