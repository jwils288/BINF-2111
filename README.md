# BINF-2111
Midterm Project for BINF 2111 - 001

Can tumor RNA-sequencing data be used to predict early vs late stage pancreatic cancer and chemotherapy response?

Background: 
Pancreatic cancer is often diagnosed at an advanced stage and can behave differently from one patient to another. Some tumors progress more aggressively, while others may respond differently to chemotherapy. These differences can be influenced by which genes are active within the tumor.

This project examines whether tumor RNA-sequencing data can be used in a two-stage machine learning approach. The first stage will analyze gene-expression patterns to estimate tumor aggressiveness, progression risk, and metastatic potential. The second stage will evaluate gene-expression signatures associated with chemotherapy sensitivity or resistance in order to predict how the tumor may respond to different treatments.

The overall goal is to produce a tailored tumor risk assessment and compare potential chemotherapy responses based on the molecular profile of the tumor. This could help improve our understanding of how gene-expression patterns relate to pancreatic cancer progression and treatment response, while supporting future research into more personalized cancer treatment strategies.

--------------------------------------
Team [Name] — Midterm Progress Report
**Course:** BINF 2111 - 001

## Team Members & Roles - Roles: (Research, Handling github, coding, overall style/keeping group on track)
* **Sesina:** [To Be Finalized]
* **Lean Saavedra Lopez:** Coding, GitHub, & Research
* **Jalen:** [To Be Finalized]
* **Melissa:** Searching for Specific Genomic Variants

--------------------------------------

**Proposed Solution**
Can genomic alterations combined with gene-expression profiles distinguish early-stage from late-stage pancreatic cancer then propose the best potential chemotherapy response?

### Tool & Coding Language
* **Language:** Python 3 (using standard libraries such as `sys` and `csv`)
* **Approach:** Filter variant lists and gene expression levels for key diagnostic markers, score each patient profile based on combined genomic features, and output stage classifications (Early vs. Late Stage).

--------------------------------------
**Approach/Pseudocode**

START  
1. Open and read the file containing tumor RNA sequencing data.
2. Open and read the patient clinical data, including cancer stage.
3. Clean and organize the gene-expression data so all patient samples can be compared fairly.
4. Label each patient as either Early-Stage or Late-Stage Pancreatic Cancer based on their known clinical stage. 
5. Look for gene-expression patterns that are different between early-stage and late-stage patients. Also look for patterns from past patients linked to chemotherapy resistance.
6. Use those gene-expression patterns to build a model that predicts cancer stage and best chemotherapy options.
7. Give the model a new patient’s RNA data.
8. IF the model predicts the patient matches late-stage patterns:
   Mark the patient as “Late-Stage Pancreatic Cancer.”
   ELSE:
   Mark the patient as “Early-Stage Pancreatic Cancer.”
9. Print the predicted stage.
10. Compare the prediction to the real patient stage and report how accurate the model is.
END
--------------------------------------
--------------------------------------

