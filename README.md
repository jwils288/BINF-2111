# BINF-2111 Midterm Project: Team Genome4!

**Course:** BINF 2111 - 001  
**Project Title:** Can tumor RNA-sequencing data be used to predict early vs. late stage pancreatic cancer and chemotherapy response?

---

## 📌 Background & Research Question
Pancreatic cancer is often diagnosed at an advanced stage and can behave differently from one patient to another. Some tumors progress more aggressively, while others may respond differently to chemotherapy. These differences can be influenced by which genes are active within the tumor.

**Primary Question:**  
*Can genomic alterations combined with gene-expression profiles distinguish early-stage from late-stage pancreatic cancer, and propose potential chemotherapy response indicators?*

This project examines whether tumor RNA-sequencing data can be evaluated in a two-stage analytical approach:
1. **Stage Classification:** Analyze gene-expression patterns to estimate tumor aggressiveness, progression risk, and stage (Early vs. Late).
2. **Chemotherapy Profile Context:** Evaluate gene-expression signatures associated with chemotherapy sensitivity or resistance to better understand treatment profiles.

---

## 👥 Team Members & Roles
* **Lean Saavedra Lopez:** Project Concept, Lead Coding (`pancreatic_analysis.py`), GitHub Management
* **Sesina:** Presentation Development
* **Jalen:** GitHub Support
* **Melissa:** Specific Genomic Variant Search (NCBI) & RNA-seq Matrix Retrieval, Research

---

## 🧬 Target Marker Genes
Our analysis evaluates gene expression variations for key diagnostic markers that are frequently overexpressed in advanced pancreatic adenocarcinoma and squamous cell carcinomas:
* **ANGPTL4** (Angiopoietin 4) — [NCBI NM_139314.3](https://www.ncbi.nlm.nih.gov/nuccore/NM_139314.3/)
* **S100P** — [NCBI NM_005980.3](https://www.ncbi.nlm.nih.gov/nuccore/NM_005980.3/)
* **KRT19** — [NCBI NM_002276.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_002276.5/)
* **COL6A3** (Collagen Type VI Alpha 3) — [NCBI NM_004369.4](https://www.ncbi.nlm.nih.gov/nuccore/NM_004369.4)
* **PLAU** — [NCBI NM_002658.6](https://www.ncbi.nlm.nih.gov/nuccore/NM_002658.6/)
* **MMP7** — [NCBI NM_002423.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_002423.5/)
* **SLC16A3** — [NCBI NM_004207.4](https://www.ncbi.nlm.nih.gov/nuccore/NM_004207.4/)

---

## 📊 Datasets Used
* **TCGA-PAAD:** Pancreatic Adenocarcinoma RNA-seq expression and clinical stage data ([GDC Portal](https://portal.gdc.cancer.gov/projects/TCGA-PAAD)).
* **GSE252327 / GSE252347:** Pancreatic Squamous Cell Carcinoma datasets evaluating tumor sub-types and treatment contexts.

---

## 🛠️ Tools & Approach
* **Language:** Python 3 (standard libraries including `sys`, `csv`, and `math` / `matplotlib` for visualizations).
* **Approach:** Filter expression matrices for target diagnostic markers, process patient stage classifications, calculate stage-specific mean expression levels, and output summary statistics.

### Pseudocode Workflow
BEGIN
    Load RNA-seq expression data and clinical patient data
    Load the target marker genes
    Match expression samples to patient clinical records

    FOR each matched patient sample
        Determine whether the sample is early-stage or late-stage
        Extract expression values for the target marker genes
    END FOR

    FOR each target gene
        Calculate mean expression in the early-stage group
        Calculate mean expression in the late-stage group
        Compare the two mean expression values
    END FOR

    Save summary statistics to CSV files
    Generate gene-expression comparison charts
END
---

## 🚀 How to Run the Analysis

1. Clone the repository:
   ```bash
   git clone [https://github.com/jwils288/BINF-2111.git](https://github.com/jwils288/BINF-2111.git)
   cd BINF-2111
