import gzip
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Define relative paths pointing to the dataset subfolder and results folder
dataset_dir = "datasets/GSE205154"
matrix_path = os.path.join(dataset_dir, "GSE205154_series_matrix.txt.gz")
tpm_path = os.path.join(dataset_dir, "GSE205154_Gene_Level_TPM_Estimates.txt.gz")

os.makedirs("results", exist_ok=True)

# 1. Parse Series Matrix Metadata for Tumor Types
print("Parsing GEO series matrix metadata...")
tumor_types = []

with gzip.open(matrix_path, "rt") as f:
  for line in f:
    if (
        line.startswith("!Sample_characteristics_ch1\t")
        and "tumor type:" in line
    ):
      parts = line.strip().split("\t")[1:]
      tumor_types = [
          p.strip('"').replace("tumor type: ", "").strip() for p in parts
      ]

# Read sample columns from expression matrix header
with gzip.open(tpm_path, "rt") as f:
  header = f.readline().strip().split("\t")
  sample_columns = header[5:]

# Map sample IDs to standardized tumor type labels
meta_dict = {}
for i, sample_id in enumerate(sample_columns):
  if i < len(tumor_types):
    raw_type = tumor_types[i]
    meta_dict[sample_id] = "Metastatic" if raw_type == "Met" else raw_type
  else:
    meta_dict[sample_id] = "Unknown"

print(f"Successfully mapped metadata for {len(meta_dict)} samples.")

# 2. Target Marker Genes
target_genes = {"ANGPTL4", "S100P", "KRT19", "COL6A3", "PLAU", "MMP7", "SLC16A3"}

# 3. Parse Expression Values
print("Parsing gene expression matrix...")
expression_rows = []
with gzip.open(tpm_path, "rt") as f:
  _ = f.readline()  # Skip header
  for line in f:
    parts = line.strip().split("\t")
    hgnc = parts[2]
    if hgnc in target_genes:
      gene_data = {
          "ensembl_gene_id": parts[1],
          "hgnc_symbol": hgnc,
          "gene_name": parts[4],
      }
      for sample_id, val in zip(sample_columns, parts[5:]):
        gene_data[sample_id] = float(val) if val != "" else 0.0
      expression_rows.append(gene_data)

df_expr = pd.DataFrame(expression_rows)

# Format into long table format for statistics
melted_data = []
for idx, row in df_expr.iterrows():
  gene = row["hgnc_symbol"]
  for sample_id in sample_columns:
    t_type = meta_dict.get(sample_id, "Unknown")
    melted_data.append({
        "gene_symbol": gene,
        "sample_id": sample_id,
        "tumor_type": t_type,
        "tpm": row[sample_id],
    })

df_long = pd.DataFrame(melted_data)

# Export processed data tables directly into results/
df_long.pivot(
    index=["gene_symbol"], columns="sample_id", values="tpm"
).reset_index().to_csv("results/gse205154_gene_expression.csv", index=False)

summary_stats = (
    df_long.groupby(["gene_symbol", "tumor_type"])["tpm"]
    .agg(
        Mean_TPM="mean", Median_TPM="median", Std_TPM="std", Sample_Count="count"
    )
    .reset_index()
)
summary_stats.to_csv("results/gse205154_gene_summary.csv", index=False)

sample_summary = (
    df_long[["sample_id", "tumor_type"]]
    .drop_duplicates()["tumor_type"]
    .value_counts()
    .reset_index()
)
sample_summary.columns = ["tumor_type", "count"]
sample_summary.to_csv("results/gse205154_sample_summary.csv", index=False)

# 4. Generate Visualization Chart
plt.figure(figsize=(12, 6))
genes = sorted(list(target_genes))
x = np.arange(len(genes))
width = 0.35
primary_means, met_means = [], []

for g in genes:
  p_mean = df_long[
      (df_long["gene_symbol"] == g) & (df_long["tumor_type"] == "Primary")
  ]["tpm"].mean()
  m_mean = df_long[
      (df_long["gene_symbol"] == g) & (df_long["tumor_type"] == "Metastatic")
  ]["tpm"].mean()
  primary_means.append(p_mean)
  met_means.append(m_mean)

plt.bar(
    x - width / 2, primary_means, width, label="Primary Tumor", color="#3498db"
)
plt.bar(
    x + width / 2, met_means, width, label="Metastatic Tumor", color="#e74c3c"
)
plt.xlabel("Target Genes", fontweight="bold")
plt.ylabel("Mean TPM Expression", fontweight="bold")
plt.title(
    "GSE205154: Primary vs. Metastatic PDAC Expression for Target Genes",
    fontweight="bold",
)
plt.xticks(x, genes)
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig("results/gse205154_expression_chart.png", dpi=300)

print("Saved all output files successfully to results/")
print("GSE205154 analysis complete!")
