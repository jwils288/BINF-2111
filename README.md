# BINF-2111
Midterm Project for BINF 2111 - 001
Can genomic alterations combined with gene-expression profiles, distinguish early-stage from late-stage pancreatic cancer?

Background: 
Pancreatic cancer is often difficult to detect and is frequently diagnosed after it has already progressed, making treatment more challenging. Understanding how the molecular characteristics of a tumor change as the cancer advances could help researchers identify patterns associated with disease progression.

This project examines whether combining genomic alterations with gene-expression profiles can distinguish early-stage from late-stage pancreatic cancer. Identifying these molecular differences could improve our understanding of how pancreatic cancer develops and support future efforts to create better biomarkers and treatment strategies.

--------------------------------------
Team [Name] — Midterm Progress Report
**Course:** BINF 2111 - 001

## Team Members & Roles - Roles: (Research, Handling github, coding, overall style/keeping group on track)
* **Sesina:** [To Be Finalized]
* **Lean Saavedra Lopez:** Coding, GitHub, & Research
* **Jalen:** [To Be Finalized]
* **Melissa:** [To Be Finalized] ---- Searching for Specific Genomic Variants

--------------------------------------
**Our Pseudocode**

**START**
  1. Open and read the file containing patient gene mutations
  2. Open and read the file containing patient gene expression data (how active each gene is)
  3. Clean up the gene numbers so all patient data can be compared fairly.
  4. Check if the patient has any major cancer-linked gene mutations.
  5. Calculate a total score based on the patient's gene activity and mutations.
  6. IF the total score is higher than our cutoff limit:
       Mark the patient as "Late-Stage Pancreatic Cancer"
     ELSE:
       Mark the patient as "Early-Stage Pancreatic Cancer"
  7. Print out the final diagnosis and report how accurate our script was.
END
