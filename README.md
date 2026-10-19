# Input data setup

The analysis script expects public TCGA-PAAD data retrieved from the [NCI Genomic Data Commons project page](https://portal.gdc.cancer.gov/projects/TCGA-PAAD).

Place files in these locations:

```text
data/
├── README.md
├── clinical.tsv
├── file_metadata.json
└── raw/
    ├── <expression-file-name-from-metadata>
    └── ...
```

## Required files

1. **`clinical.tsv`** — clinical TSV export with columns named `cases.submitter_id` and `diagnoses.ajcc_pathologic_stage`.
2. **`file_metadata.json`** — metadata JSON export with each expression file's `file_name`, associated `cases`, patient `submitter_id`, and sample metadata including `sample_type`.
3. **Expression files** — the STAR-counts TSV files named in `file_metadata.json`, placed in `data/raw/`.

When retrieving files, choose RNA-seq gene-expression files compatible with the parser in `pancreatic_analysis.py`, which expects a GDC STAR-counts table with a gene name in the second column and unstranded TPM in the seventh column (index 6). Check a downloaded file's header before running the script; GDC file formats can vary.

The script also supports `clinical.tsv` and `file_metadata.json` placed beside `pancreatic_analysis.py` for backward compatibility. Do not commit large raw data files to GitHub. If your team shares small metadata or clinical exports, follow your instructor's requirements and the data source's terms. Record the download date, filters, and source in your project notes so another team member can reproduce the selection.
