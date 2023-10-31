# Code for the PS4 paper
## Background
The American College of Medical Genetics and Genomics and Association for Molecular Pathology joint guidelines have been widely adapted for variant classification in genetic testing of Mendelian disorders. Meanwhile, rules in the guidelines are continuously developed to improve pathogenicity interpretation for accurate and consistent clinical application. One specific criterion in the ACMG/AMP guidelines is PS4, designated as 'Strong' level evidence. This pathogenicity criterion refers to a significantly increased prevalence of a variant in affected individuals compared to ancestry-matched controls. However, quantitative support for PS4 thresholds from real-world Mendelian case-control cohorts is lacking.

To address this gap, we evaluated and defined the PS4 thresholds using a large cohort dataset including 13,845 cases with disabling hearing loss (HL) and 6,570 ancestry-matched controls from the Chinese Deafness Genetics Consortium (CDGC) project. As shown in Figure 1, positive likelihood ratio and local positive likelihood ratio values were calculated to determine the thresholds corresponding to each strength of evidence across three variant subsets: subset 1, consisting of variants present in both cases and controls with an allele frequency (AF) in cases ≥ 0.0005; subset 2, which encompassed variants present in both cases and controls with a case AF < 0.0005; and subset 3, comprising variants found only in cases and absent from controls.

![Evaluation workflow](https://github.com/liusihan/PS4_paper_2023/blob/main/Figure1.png)
<p align="center"> Figure 1. Study design. </p>




## Repository structure
This repository contains the code relevant to the PS4 paper. All data relevant to the PS4 paper are included in the article or uploaded as Supplemental Material. Additional data are available from the corresponding author upon request. The bootstrap result files for the truth subset 2 are hosted in the `result` directory. 

This `code` directory contains the actual implementation of the algorithm to calculate posterior probabilities and local posterior probabilities. Additionally, the code used to make the plots in the paper was stored in the `Figure_plot.R`.

1. **`quantitative_thresholds.R`**: R Script for calculating the **`posterior probability`** and the **`local posterior probability`**. We defined the function of `VUS_classify` to classify VUS into 6 levels. Also, the function of `LR` was utilized for calculating Prevalence, Accuracy, PPV, NPV, F1 score, Sensitivity, Specificity, and posterior probability in truth subset 1 and truth subset 3. The local posterior probability was first developed by Pejaver, Vikas et al. to evaluate and define the thresholds for PP3 and BP4. We described the `local_lr` and `local_bootstrapped_lr` functions to calculate the local posterior probability and the one-sided 95% confidence bound for each estimated lr+. Briefly, the posterior probability was calculated for each OR value within the interval, considering a minimum of 100 pathogenic and non-pathogenic variants. Additionally, the one-sided 95% confidence bound for each estimated lr+  was determined through 10,000 bootstrapping iterations, enabling the assessment of evidence strength.  


2. **`Figure_plot.R`**: R Script for generating the main and supplementary figures in the PS4 paper.
