import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Parse Series Matrix Metadata (matching Sample ID to Tumor Type)
titles = []
tumor_types = []

print("Parsing GEO series matrix metadata...")
with gzip.open("GSE205154_series_matrix.txt.gz", "rt") as f:
  for line in f:
    if line.startswith("!Sample_title\t"):
      titles = [s.strip('"') for s in line.strip().split("\t")[1:]]
    elif "tumor type:" in line:
      parts = line.strip().split("\t")[1:]
      candidate_types = [
          p.strip('"').replace("tumor type: ", "")
          for p in parts
          if "tumor type:" in p
      ]
      if len(candidate_types) == len(titles):
        tumor_types = candidate_types

# Fallback line parsing
if len(titles) != len(tumor_types):
  print("Using secondary parsing strategy for metadata alignment...")
  with gzip.open("GSE205154_series_matrix.txt.gz", "rt") as f:
    for line in f:
      if line.startswith("!Sample_characteristics_ch1\t") and (
          "tumor type:" in line
      ):
        parts = line.strip().split("\t")[1:]
        tumor_types = [p.strip('"').replace("tumor type: ", "") for p in parts]

meta_dict = {}
for s_id, t_type in zip(titles, tumor_types):
  meta_dict[s_id] = t_type

print(f"Successfully mapped metadata for {len(meta_dict)} samples.")

# 2. Target genes to extract
target_genes = {"ANGPTL4", "S100P", "KRT19", "COL6A3", "PLAU", "MMP7", "SLC16A3"}

# 3. Parse Gene Expression Matrix
print("Parsing gene expression matrix...")
expression_rows = []
sample_columns = []

with gzip.open("GSE205154_Gene_Level_TPM_Estimates.txt.gz", "rt") as f:
  header = f.readline().strip().split("\t")
  sample_columns = header[5:]  # First 5 columns are gene metadata

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

# Reshape into long format
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

# Save processed tables
df_long.pivot(
    index=["gene_symbol"], columns="sample_id", values="tpm"
).reset_index().to_csv("gse205154_gene_expression.csv", index=False)
print("Saved: gse205154_gene_expression.csv")

summary_stats = (
    df_long.groupby(["gene_symbol", "tumor_type"])["tpm"]
    .agg(
        Mean_TPM="mean",
        Median_TPM="median",
        Std_TPM="std",
        Sample_Count="count",
    )
    .reset_index()
)
summary_stats.to_csv("gse205154_gene_summary.csv", index=False)
print("Saved: gse205154_gene_summary.csv")

sample_summary = (
    df_long[["sample_id", "tumor_type"]]
    .drop_duplicates()["tumor_type"]
    .value_counts()
    .reset_index()
)
sample_summary.columns = ["tumor_type", "count"]
sample_summary.to_csv("gse205154_sample_summary.csv", index=False)
print("Saved: gse205154_sample_summary.csv")

# 4. Generate Comparison Visualization
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
      (df_long["gene_symbol"] == g) & (df_long["tumor_type"] == "Met")
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

plt.savefig("gse205154_expression_chart.png", dpi=300)
print("Saved: gse205154_expression_chart.png")
print("GSE205154 Analysis pipeline complete!")
