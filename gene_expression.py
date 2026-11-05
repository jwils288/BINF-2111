#!/usr/bin/env python3
"""Explore RNA-seq expression patterns associated with pancreatic cancer stage.

Reads GDC-style STAR-counts files plus file_metadata.json and clinical.tsv
(TCGA-PAAD). Steps:

  1. Match one primary-tumor RNA-seq sample per patient to a recorded stage and
     group patients as Early-stage vs Late-stage.
  2. Remove low-expression genes (TPM filter).
  3. Test every remaining gene for an early-vs-late difference
     (Mann-Whitney U test on log2(TPM + 1)) with Benjamini-Hochberg FDR
     correction, and report where the 7 project target genes rank.
  4. Evaluate a 7-gene model and a genome-wide model with repeated, stratified,
     leakage-aware cross-validation.
  5. Save CSV tables and log2-scale figures to results/.

This is a research/educational script, not a validated diagnostic or treatment
recommendation tool.

Expected input files beside this script:
  file_metadata.json, clinical.tsv, and the RNA-seq files referenced by metadata.

Dependencies: numpy, pandas, scipy, matplotlib, scikit-learn.
"""

import json
import re
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend: works in WSL / servers without a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE_DIR = Path(__file__).resolve().parent
METADATA_FILE = BASE_DIR / "file_metadata.json"
CLINICAL_FILE = BASE_DIR / "clinical.tsv"
OUTPUT_DIR = BASE_DIR / "results"

# ---------------------------------------------------------------- parameters
# Seven markers listed in the project presentation.
TARGET_GENES = [
    "ANGPTL4", "S100P", "KRT19", "COL6A3", "PLAU", "MMP7", "SLC16A3"
]

# Stage grouping. {"0", "I"} -> Stage I vs Stages II-IV.
# Use {"0", "I", "II"} for Stages I-II vs III-IV instead.
EARLY_STAGE_NUMERALS = {"0", "I"}

# Keep only protein-coding genes (uses the gene_type column when present).
PROTEIN_CODING_ONLY = True

# Low-expression filter: keep a gene only if TPM >= MIN_TPM in at least
# MIN_SAMPLE_FRACTION of the matched patients. Target genes are always tested.
MIN_TPM = 1.0
MIN_SAMPLE_FRACTION = 0.20

# Statistics
MIN_PATIENTS_PER_GROUP = 3
FDR_THRESHOLD = 0.05
TOP_N_GENES_TO_REPORT = 100

# Machine learning
MAX_MODEL_FEATURES = 50
CV_FOLDS = 5
CV_REPEATS = 10
RANDOM_STATE = 42

EARLY, LATE = "Early-stage", "Late-stage"


# ------------------------------------------------------------ clinical data
def load_clinical_stages(path):
    """Return one unambiguous recorded stage per patient and conflicting IDs."""
    clinical = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    patient_col = "cases.submitter_id"
    stage_col = "diagnoses.ajcc_pathologic_stage"
    missing = [col for col in (patient_col, stage_col) if col not in clinical.columns]
    if missing:
        raise ValueError(
            f"Missing expected column(s) in clinical.tsv: {missing}. "
            "Check the GDC clinical export column names."
        )

    stages_by_patient = defaultdict(set)
    for _, row in clinical.iterrows():
        patient = row[patient_col].strip()
        stage = row[stage_col].strip()
        if not patient or stage.upper() in {"", "--", "'--", "NA", "NAN", "NOT REPORTED"}:
            continue
        stages_by_patient[patient].add(stage)

    patient_stages = {}
    conflicting = set()
    for patient, stages in stages_by_patient.items():
        if len(stages) == 1:
            patient_stages[patient] = next(iter(stages))
        elif len(stages) > 1:
            conflicting.add(patient)
    return patient_stages, conflicting


def stage_group(stage):
    """Map recorded stage labels to the project's early/late comparison."""
    if not stage:
        return None
    match = re.search(r"\bStage\s+(IV|III|II|I|0)(?:[ABC])?\b", stage, re.IGNORECASE)
    if not match:
        return None
    roman = match.group(1).upper()
    return EARLY if roman in EARLY_STAGE_NUMERALS else LATE


# ------------------------------------------------------------ expression IO
def locate_data_file(file_name):
    """Find a metadata-listed file at the project root or one directory below."""
    if not file_name:
        return None
    direct = BASE_DIR / file_name
    if direct.is_file():
        return direct
    candidates = list(BASE_DIR.glob("*/" + Path(file_name).name))
    return candidates[0] if candidates else None


def read_gene_expression(path, protein_coding_only=PROTEIN_CODING_ONLY):
    """Read gene names and the TPM column from a GDC STAR-counts TSV.

    Fails clearly if no TPM column exists rather than silently reading another
    count column. If several Ensembl IDs share one gene name, their TPM values
    are summed so no gene is silently overwritten.
    """
    values = {}
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        tpm_index = None
        gene_name_index = None
        gene_type_index = None
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "gene_id":
                header = fields
                gene_name_index = header.index("gene_name") if "gene_name" in header else 1
                gene_type_index = header.index("gene_type") if "gene_type" in header else None
                if "tpm_unstranded" in header:
                    tpm_index = header.index("tpm_unstranded")
                elif "TPM" in header:
                    tpm_index = header.index("TPM")
                else:
                    raise ValueError(
                        f"No TPM column found in {path}. Columns: {header}. "
                        "Please confirm this is a GDC STAR-counts file."
                    )
                continue
            if header is None or len(fields) < len(header):
                continue
            if (protein_coding_only and gene_type_index is not None
                    and fields[gene_type_index].strip() != "protein_coding"):
                continue
            gene_name = fields[gene_name_index].strip()
            if not gene_name or gene_name == "N/A":
                continue
            try:
                tpm = float(fields[tpm_index])
            except (ValueError, IndexError):
                continue
            if np.isfinite(tpm) and tpm >= 0:
                values[gene_name] = values.get(gene_name, 0.0) + tpm
    return values


def extract_patient_samples(file_record):
    """Return [(patient_id, sample_id, is_primary_tumor)] for one metadata record.

    Supports two GDC metadata layouts:
      1. nested  : record["cases"][i]["samples"][j]["sample_type"]
      2. cart    : record["associated_entities"][i]["entity_submitter_id"], a TCGA
                   barcode such as TCGA-2J-AAB1-01A-11R-A41B-07. The two digits
                   after the patient ID are the sample-type code (01 = primary tumor).
    """
    results = []
    for case in file_record.get("cases", []) or []:
        patient_id = case.get("submitter_id", "").strip()
        primary = [
            sample for sample in case.get("samples", []) or []
            if sample.get("sample_type", "").strip().lower() == "primary tumor"
        ]
        sample_id = primary[0].get("submitter_id", patient_id) if primary else None
        results.append((patient_id, sample_id, bool(primary)))
    if results:
        return results

    for entity in file_record.get("associated_entities", []) or []:
        barcode = entity.get("entity_submitter_id", "").strip()
        match = re.match(r"^(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4})-(\d{2})", barcode)
        if match:
            results.append((match.group(1), barcode[:16], match.group(2) == "01"))
    return results


def collect_samples(metadata, clinical_stages, conflicting_patients):
    """Match one primary-tumor RNA-seq sample per patient to a stage label."""
    rows = []
    stats = {
        "metadata_records": len(metadata),
        "primary_tumor_files_checked": 0,
        "missing_files": 0,
        "no_patient_record": 0,
        "missing_or_unusable_stage": 0,
        "conflicting_stage": 0,
        "duplicate_patient_samples_skipped": 0,
    }
    seen_patients = set()

    for file_record in metadata:
        filepath = locate_data_file(file_record.get("file_name", ""))
        if filepath is None:
            stats["missing_files"] += 1
            continue
        entries = extract_patient_samples(file_record)
        if not entries:
            stats["no_patient_record"] += 1
            continue

        for patient_id, sample_id, is_primary in entries:
            if not is_primary:
                continue
            stats["primary_tumor_files_checked"] += 1

            if not patient_id:
                stats["no_patient_record"] += 1
                continue
            if patient_id in seen_patients:
                stats["duplicate_patient_samples_skipped"] += 1
                continue
            if patient_id in conflicting_patients:
                stats["conflicting_stage"] += 1
                continue
            stage = clinical_stages.get(patient_id)
            group = stage_group(stage) if stage else None
            if group is None:
                stats["missing_or_unusable_stage"] += 1
                continue

            expression = read_gene_expression(filepath)
            row = {
                "patient_id": patient_id,
                "sample_id": sample_id or patient_id,
                "recorded_stage": stage,
                "stage_group": group,
            }
            row.update(expression)
            rows.append(row)
            seen_patients.add(patient_id)
    return rows, stats


# --------------------------------------------------------------- statistics
def filter_low_expression(tpm):
    """Return genes with TPM >= MIN_TPM in at least MIN_SAMPLE_FRACTION of patients."""
    detected_fraction = (tpm >= MIN_TPM).mean(axis=0)  # NaN compares as False
    return detected_fraction[detected_fraction >= MIN_SAMPLE_FRACTION].index.tolist()


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR-adjusted p-values (NaNs are ignored)."""
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return adjusted
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    result = np.empty(n)
    result[order] = np.clip(ranked, 0, 1)
    adjusted[valid] = result
    return adjusted


def differential_expression(tpm, log_expr, groups, genes, passed_filter):
    """Mann-Whitney U test per gene (late vs early) with BH-adjusted p-values.

    log2_diff_late_minus_early is the difference of mean log2(TPM + 1) values,
    which approximates a log2 fold change. auc_late_vs_early is the probability
    that a random late-stage patient has higher expression than a random
    early-stage patient (0.5 = no separation).
    """
    late_mask = (groups == LATE).to_numpy()
    records = []
    for gene in genes:
        x = log_expr[gene].to_numpy(dtype=float)
        t = tpm[gene].to_numpy(dtype=float)
        late, early = x[late_mask], x[~late_mask]
        late, early = late[~np.isnan(late)], early[~np.isnan(early)]
        late_t, early_t = t[late_mask], t[~late_mask]
        late_t, early_t = late_t[~np.isnan(late_t)], early_t[~np.isnan(early_t)]

        rec = {
            "gene": gene,
            "passed_expression_filter": gene in passed_filter,
            "early_n": len(early), "late_n": len(late),
            "early_mean_TPM": float(early_t.mean()) if len(early_t) else np.nan,
            "late_mean_TPM": float(late_t.mean()) if len(late_t) else np.nan,
            "early_median_TPM": float(np.median(early_t)) if len(early_t) else np.nan,
            "late_median_TPM": float(np.median(late_t)) if len(late_t) else np.nan,
            "log2_diff_late_minus_early": np.nan,
            "auc_late_vs_early": np.nan,
            "p_value": np.nan,
        }
        if len(late) >= MIN_PATIENTS_PER_GROUP and len(early) >= MIN_PATIENTS_PER_GROUP:
            rec["log2_diff_late_minus_early"] = float(late.mean() - early.mean())
            try:
                u_stat, p_val = mannwhitneyu(late, early, alternative="two-sided")
                rec["auc_late_vs_early"] = float(u_stat / (len(late) * len(early)))
                rec["p_value"] = float(p_val)
            except ValueError:  # all values identical
                rec["auc_late_vs_early"] = 0.5
                rec["p_value"] = 1.0
        records.append(rec)

    de = pd.DataFrame(records)
    de["fdr_bh"] = benjamini_hochberg(de["p_value"].to_numpy())
    de["significant_fdr"] = de["fdr_bh"] < FDR_THRESHOLD
    de["direction"] = np.where(
        de["log2_diff_late_minus_early"] > 0, "higher in late",
        np.where(de["log2_diff_late_minus_early"] < 0, "higher in early", "")
    )
    de["abs_log2_diff"] = de["log2_diff_late_minus_early"].abs()
    de = de.sort_values(["p_value", "abs_log2_diff"], ascending=[True, False],
                        na_position="last").reset_index(drop=True)
    de["rank_by_p_value"] = np.where(de["p_value"].notna(), np.arange(1, len(de) + 1), np.nan)
    return de


# ----------------------------------------------------------- classification
def run_cv(log_expr, groups, gene_list, label):
    """Repeated stratified k-fold CV. Gene selection, imputation and scaling are
    all fitted inside each training fold, so test folds never influence them.
    Returns (summary dict or None, fitted-on-all-data pipeline or None, features)."""
    try:
        from sklearn.feature_selection import SelectKBest, f_classif
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import make_scorer, recall_score
        from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                             cross_val_predict, cross_validate)
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print(f"\n[{label}] skipped: install scikit-learn to enable cross-validation.")
        return None, None, [], None, None

    X = log_expr[gene_list].dropna(axis=1, how="all")
    features = X.columns.tolist()
    y = (groups == LATE).astype(int)  # 1 = Late-stage, 0 = Early-stage
    class_counts = y.value_counts()
    if len(features) < 2 or len(class_counts) != 2 or class_counts.min() < 3:
        print(f"\n[{label}] skipped: need >= 2 usable genes and >= 3 patients per stage group "
              f"(counts: {groups.value_counts().to_dict()}).")
        return None, None, features, None, None

    n_splits = int(min(CV_FOLDS, class_counts.min()))
    k = min(MAX_MODEL_FEATURES, len(features))

    def make_pipeline():
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("select_genes", SelectKBest(score_func=f_classif, k=k)),
            ("scale", StandardScaler()),
            ("classifier", LogisticRegression(class_weight="balanced", max_iter=3000,
                                               random_state=RANDOM_STATE)),
        ])

    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=CV_REPEATS,
                                 random_state=RANDOM_STATE)
    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "roc_auc": "roc_auc",
        "late_recall": make_scorer(recall_score, pos_label=1, zero_division=0),
        "early_recall": make_scorer(recall_score, pos_label=0, zero_division=0),
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # constant-feature warnings from f_classif
        scores = cross_validate(make_pipeline(), X, y, cv=cv, scoring=scoring,
                                n_jobs=1, error_score="raise")

    summary = {
        "model": label,
        "features_available": len(features),
        "features_selected_per_fold": k,
        "patients": int(len(y)),
        "early_patients": int((y == 0).sum()),
        "late_patients": int((y == 1).sum()),
        "cv_folds": n_splits,
        "cv_repeats": CV_REPEATS,
        "majority_class_accuracy": float(class_counts.max() / class_counts.sum()),
    }
    for name in scoring:
        vals = scores[f"test_{name}"]
        summary[f"{name}_mean"] = float(np.mean(vals))
        summary[f"{name}_sd"] = float(np.std(vals))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final_pipeline = make_pipeline().fit(X, y)
        # One set of out-of-fold probabilities (each patient predicted by a model
        # that never saw them) for the ROC curve and confusion matrix figures.
        oof_prob = cross_val_predict(
            make_pipeline(), X, y,
            cv=StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE),
            method="predict_proba")[:, 1]
    return summary, final_pipeline, features, pd.Series(oof_prob, index=y.index), y


def save_seven_gene_coefficients(pipeline, features):
    """Standardized logistic-regression coefficients (positive = late-stage)."""
    mask = pipeline.named_steps["select_genes"].get_support()
    genes = [g for g, keep in zip(features, mask) if keep]
    coefs = pipeline.named_steps["classifier"].coef_[0]
    out = pd.DataFrame({"gene": genes, "standardized_coefficient": coefs})
    out["direction"] = np.where(out["standardized_coefficient"] > 0,
                                "pushes toward late", "pushes toward early")
    out = out.reindex(out["standardized_coefficient"].abs().sort_values(ascending=False).index)
    out.to_csv(OUTPUT_DIR / "seven_gene_model_coefficients.csv", index=False)
    return out


# -------------------------------------------------------------------- plots
def significance_label(fdr):
    if pd.isna(fdr):
        return "n/a"
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < FDR_THRESHOLD:
        return "*"
    return "ns"


def plot_target_genes(log_expr, groups, target_de, path):
    genes = [g for g in TARGET_GENES if g in log_expr.columns]
    if not genes:
        return
    early = log_expr.loc[groups == EARLY, genes]
    late = log_expr.loc[groups == LATE, genes]
    x = np.arange(len(genes))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(x - width / 2, early.mean(), width, yerr=early.sem(), capsize=3,
           label=f"Early-stage (n={len(early)})")
    ax.bar(x + width / 2, late.mean(), width, yerr=late.sem(), capsize=3,
           label=f"Late-stage (n={len(late)})")
    fdr = target_de.set_index("gene")["fdr_bh"]
    top = np.maximum(early.mean() + early.sem().fillna(0),
                     late.mean() + late.sem().fillna(0))
    for i, gene in enumerate(genes):
        ax.text(i, top[gene] + 0.1, significance_label(fdr.get(gene, np.nan)),
                ha="center", va="bottom", fontsize=10)
    ax.set_ylim(0, float(np.nanmax(top.to_numpy())) * 1.15 + 0.3)  # headroom for labels
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=35, ha="right")
    ax.set_ylabel("Mean log2(TPM + 1) \u00b1 SEM")
    ax.set_xlabel("Target gene")
    ax.set_title("Target gene expression by stage group (TCGA-PAAD)\n"
                 "stars = BH-adjusted FDR (* <0.05, ** <0.01, *** <0.001, ns = not significant)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_volcano(de, path):
    tested = de.dropna(subset=["p_value", "log2_diff_late_minus_early"]).copy()
    if tested.empty:
        return
    tested["neg_log10_fdr"] = -np.log10(tested["fdr_bh"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(9, 7))
    sig = tested["significant_fdr"]
    ax.scatter(tested.loc[~sig, "log2_diff_late_minus_early"],
               tested.loc[~sig, "neg_log10_fdr"], s=6, alpha=0.4, color="gray",
               label="Not significant")
    ax.scatter(tested.loc[sig, "log2_diff_late_minus_early"],
               tested.loc[sig, "neg_log10_fdr"], s=8, alpha=0.7, color="tab:red",
               label=f"FDR < {FDR_THRESHOLD}")
    ax.axhline(-np.log10(FDR_THRESHOLD), linestyle="--", linewidth=0.8, color="black")
    for _, row in tested[tested["gene"].isin(TARGET_GENES)].iterrows():
        ax.scatter(row["log2_diff_late_minus_early"], row["neg_log10_fdr"],
                   s=40, color="tab:blue", zorder=3)
        ax.annotate(row["gene"], (row["log2_diff_late_minus_early"], row["neg_log10_fdr"]),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Late-stage minus early-stage mean log2(TPM + 1)")
    ax.set_ylabel("-log10(BH-adjusted FDR)")
    ax.set_title("Genome-wide stage differences (blue = project target genes)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_top_genes(de, path, n=20):
    top = de.dropna(subset=["p_value"]).head(n).sort_values("log2_diff_late_minus_early")
    if top.empty:
        return
    colors = np.where(top["log2_diff_late_minus_early"] > 0, "tab:red", "tab:blue")
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(top["gene"], top["log2_diff_late_minus_early"], color=colors)
    ax.set_xlabel("Late-stage minus early-stage mean log2(TPM + 1)")
    ax.set_ylabel("Gene")
    ax.set_title(f"Top {len(top)} genes by p-value (red = higher in late, blue = higher in early)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_pca(log_expr, groups, gene_sets, path):
    """PCA of patients (standardized log2 TPM), one panel per gene set."""
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        return
    panels = [(title, genes) for title, genes in gene_sets if len(genes) >= 2]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(6.5 * len(panels), 5.5), squeeze=False)
    for ax, (title, genes) in zip(axes[0], panels):
        X = log_expr[genes]
        X = X.fillna(X.median())
        X = X.loc[:, X.std() > 0]
        pca = PCA(n_components=2, random_state=RANDOM_STATE)
        coords = pca.fit_transform(StandardScaler().fit_transform(X))
        for label, color in ((EARLY, "tab:blue"), (LATE, "tab:red")):
            mask = (groups == label).to_numpy()
            ax.scatter(coords[mask, 0], coords[mask, 1], s=22, alpha=0.7, color=color,
                       label=f"{label} (n={int(mask.sum())})")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
        ax.set_title(title)
        ax.legend()
    fig.suptitle("Do patients separate by stage group? (PCA of log2 TPM)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_roc_curves(model_results, path):
    """ROC curves from out-of-fold predictions (model_results: label -> (y, prob))."""
    try:
        from sklearn.metrics import roc_auc_score, roc_curve
    except ImportError:
        return
    if not model_results:
        return
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for label, (y, prob) in model_results.items():
        fpr, tpr, _ = roc_curve(y, prob)
        ax.plot(fpr, tpr, linewidth=2, label=f"{label} (AUC = {roc_auc_score(y, prob):.2f})")
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="Chance (AUC = 0.50)")
    ax.set_xlabel("False positive rate (early patients called late)")
    ax.set_ylabel("True positive rate (late patients called late)")
    ax.set_title("Stage classifiers: cross-validated ROC")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_confusion_matrices(model_results, path):
    """Row-normalized confusion matrices from out-of-fold predictions."""
    try:
        from sklearn.metrics import confusion_matrix
    except ImportError:
        return
    if not model_results:
        return
    fig, axes = plt.subplots(1, len(model_results), figsize=(5.2 * len(model_results), 4.8),
                             squeeze=False)
    for ax, (label, (y, prob)) in zip(axes[0], model_results.items()):
        cm = confusion_matrix(y, (prob >= 0.5).astype(int), labels=[0, 1])
        pct = cm / cm.sum(axis=1, keepdims=True)
        ax.imshow(pct, vmin=0, vmax=1, cmap="Blues")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]}\n({pct[i, j]:.0%})", ha="center", va="center",
                        color="white" if pct[i, j] > 0.5 else "black")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Predicted early", "Predicted late"])
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["Actual early", "Actual late"])
        ax.set_title(label, fontsize=9)
    fig.suptitle("Cross-validated confusion matrices (% of each actual group)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_coefficients(coef_df, path):
    """Which of the 7 genes push the model toward late vs early stage."""
    if coef_df is None or coef_df.empty:
        return
    d = coef_df.sort_values("standardized_coefficient")
    colors = np.where(d["standardized_coefficient"] > 0, "tab:red", "tab:blue")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(d["gene"], d["standardized_coefficient"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Standardized logistic-regression coefficient")
    ax.set_title("7-gene model: red = toward late-stage, blue = toward early-stage")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------- main
def main():
    if not METADATA_FILE.exists() or not CLINICAL_FILE.exists():
        raise FileNotFoundError(
            "Required files not found. Keep file_metadata.json, clinical.tsv, "
            "the RNA-seq files, and this script in the same project folder "
            "(or place the RNA-seq files one directory below it)."
        )

    with open(METADATA_FILE, "r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    clinical_stages, conflicting_patients = load_clinical_stages(CLINICAL_FILE)
    OUTPUT_DIR.mkdir(exist_ok=True)

    rows, stats = collect_samples(metadata, clinical_stages, conflicting_patients)
    print("=" * 64)
    print("Pancreatic Cancer Gene Expression and Stage Analysis")
    print("Dataset: TCGA-PAAD / GDC-style RNA-seq files")
    print("=" * 64)
    for key, value in stats.items():
        print(f"{key.replace('_', ' ').capitalize()}: {value}")
    print("Patients with conflicting clinical stages:", len(conflicting_patients))

    if not rows:
        print("No samples could be matched to unambiguous recorded stages.")
        print("Check the clinical column names, metadata paths, and sample labels.")
        return

    # ---- build matrices (patients x genes), indexed by patient ID
    df = pd.DataFrame(rows).set_index("patient_id", drop=False)
    meta_cols = ["patient_id", "sample_id", "recorded_stage", "stage_group"]
    groups = df["stage_group"]
    tpm = df.drop(columns=meta_cols).apply(pd.to_numeric, errors="coerce")
    for gene in TARGET_GENES:
        if gene not in tpm.columns:
            tpm[gene] = np.nan
    log_expr = np.log2(tpm + 1.0)

    df[meta_cols].join(tpm).to_csv(OUTPUT_DIR / "matched_expression_by_stage.csv", index=False)
    df["recorded_stage"].value_counts().rename_axis("recorded_stage").reset_index(
        name="patients").to_csv(OUTPUT_DIR / "stage_counts.csv", index=False)

    n_early = int((groups == EARLY).sum())
    n_late = int((groups == LATE).sum())
    print("\nUnique patients analyzed:", len(df))
    print(f"Early-stage patients: {n_early}   Late-stage patients: {n_late}")
    print("Recorded stage breakdown:", df["recorded_stage"].value_counts().to_dict())
    if min(n_early, n_late) < 10:
        print("WARNING: fewer than 10 patients in one stage group. "
              "Tests and CV results will be very unstable.")

    # ---- low-expression filter
    kept_genes = filter_low_expression(tpm)
    print(f"\nGenes in matched files: {tpm.shape[1]}")
    print(f"Genes kept after low-expression filter (TPM >= {MIN_TPM} in "
          f">= {MIN_SAMPLE_FRACTION:.0%} of patients): {len(kept_genes)}")
    targets_present = [g for g in TARGET_GENES if tpm[g].notna().any()]
    missing_targets = [g for g in TARGET_GENES if g not in targets_present]
    if missing_targets:
        print("Target genes not found in the expression files:", missing_targets)
    dropped_targets = [g for g in targets_present if g not in kept_genes]
    if dropped_targets:
        print("Target genes below the expression filter (still tested, flagged):", dropped_targets)

    # ---- differential expression with FDR correction
    genes_to_test = sorted(set(kept_genes) | set(targets_present))
    de = differential_expression(tpm, log_expr, groups, genes_to_test, set(kept_genes))
    de.to_csv(OUTPUT_DIR / "differential_expression_all_genes.csv", index=False)
    de.head(TOP_N_GENES_TO_REPORT).to_csv(OUTPUT_DIR / "top_gene_stage_differences.csv", index=False)

    target_de = de[de["gene"].isin(TARGET_GENES)].copy()
    target_de["genes_tested"] = int(de["p_value"].notna().sum())
    target_de.to_csv(OUTPUT_DIR / "target_gene_stage_summary.csv", index=False)

    n_sig = int(de["significant_fdr"].sum())
    print(f"\nGenes tested: {int(de['p_value'].notna().sum())}   "
          f"significant at FDR < {FDR_THRESHOLD}: {n_sig}")
    print("\nTarget-gene results (Mann-Whitney U, BH-adjusted):")
    show = target_de[["gene", "early_mean_TPM", "late_mean_TPM", "log2_diff_late_minus_early",
                      "auc_late_vs_early", "p_value", "fdr_bh", "rank_by_p_value"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"(rank is out of {int(de['p_value'].notna().sum())} tested genes; "
          "1 = strongest stage association)")

    # ---- plots (all on log2(TPM + 1) scale)
    plot_target_genes(log_expr, groups, target_de, OUTPUT_DIR / "gene_expression_by_stage.png")
    plot_volcano(de, OUTPUT_DIR / "volcano_stage_differences.png")
    plot_top_genes(de, OUTPUT_DIR / "top_gene_stage_differences.png")
    print("\nFigures saved in:", OUTPUT_DIR)

    # ---- cross-validated classification
    summaries = []
    model_results = {}
    seven_summary, seven_pipe, seven_features, seven_oof, seven_y = run_cv(
        log_expr, groups, targets_present, "7-gene model (LogisticRegression)")
    if seven_summary:
        summaries.append(seven_summary)
        coef_df = save_seven_gene_coefficients(seven_pipe, seven_features)
        plot_coefficients(coef_df, OUTPUT_DIR / "seven_gene_coefficients.png")
        model_results["7-gene model"] = (seven_y, seven_oof)
    genome_summary, genome_pipe, genome_features, genome_oof, genome_y = run_cv(
        log_expr, groups, kept_genes,
        f"Genome-wide model (top {MAX_MODEL_FEATURES} genes per fold, LogisticRegression)")
    if genome_summary:
        summaries.append(genome_summary)
        model_results[f"Genome-wide model (top {MAX_MODEL_FEATURES} genes)"] = (genome_y, genome_oof)
        mask = genome_pipe.named_steps["select_genes"].get_support()
        pd.DataFrame({"selected_gene": [g for g, k in zip(genome_features, mask) if k]}).to_csv(
            OUTPUT_DIR / "genome_model_selected_genes_full_data.csv", index=False)

    plot_roc_curves(model_results, OUTPUT_DIR / "model_roc_curves.png")
    plot_confusion_matrices(model_results, OUTPUT_DIR / "model_confusion_matrices.png")
    plot_pca(log_expr, groups,
             [("7 target genes", targets_present), ("All filtered genes", kept_genes)],
             OUTPUT_DIR / "model_pca_by_stage.png")

    if summaries:
        pd.DataFrame(summaries).to_csv(OUTPUT_DIR / "stage_model_cv_metrics.csv", index=False)
        print("\nCross-validated results "
              f"({CV_REPEATS}x repeated stratified CV; mean +/- SD across folds):")
        for s in summaries:
            print(f"\n  {s['model']}")
            print(f"    patients: {s['patients']} (early {s['early_patients']}, late {s['late_patients']}), "
                  f"{s['cv_folds']}-fold")
            print(f"    majority-class accuracy baseline: {s['majority_class_accuracy']:.3f}")
            for name in ("accuracy", "balanced_accuracy", "roc_auc", "late_recall", "early_recall"):
                print(f"    {name}: {s[name + '_mean']:.3f} +/- {s[name + '_sd']:.3f}")

    with open(OUTPUT_DIR / "run_parameters.json", "w", encoding="utf-8") as handle:
        json.dump({
            "early_stage_numerals": sorted(EARLY_STAGE_NUMERALS),
            "protein_coding_only": PROTEIN_CODING_ONLY,
            "min_tpm": MIN_TPM, "min_sample_fraction": MIN_SAMPLE_FRACTION,
            "fdr_threshold": FDR_THRESHOLD, "cv_folds": CV_FOLDS, "cv_repeats": CV_REPEATS,
            "max_model_features": MAX_MODEL_FEATURES, "random_state": RANDOM_STATE,
            "patients": len(df), "early_patients": n_early, "late_patients": n_late,
            "genes_after_filter": len(kept_genes),
        }, handle, indent=2)

    print("\nImportant limitations:")
    print("- This script analyzes TCGA-PAAD inputs only; it does not automatically load GEO data.")
    print("- It predicts only early vs late stage, not individual stages.")
    print("- It does not predict chemotherapy response because treatment-response labels are not loaded.")
    print("- Stage differences can be confounded by tumor purity, batch, and treatment history.")
    print("- Results are exploratory and need external validation; they are not clinical advice.")


if __name__ == "__main__":
    main()
