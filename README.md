# GWAS-program
A script that reads DNA data in input VCF files, compares it to a physical trait, adjusts for age/sex, and outputs a list of relevant genetic markers.

An example command line execution will be: 
!python3/gwas.py
-g chr4.1000G.bcf 
-p phenotype.txt
-c covariates.txt
-o chr4.gwas.tsv.gz

Required Command Line Arguments
1. **Genotype Source File**

**Input:** A VCF or BCF format file.

Description: This file serves as the primary dataset containing the multi-sample genotypes. The script must be able to parse this file to extract the dosage or state of each genetic variant for every individual included in the study.

2. **Quantitative Phenotype Data**

**Input:** A tab-delimited text file.

Description: This file provides the dependent variable for the regression analysis. It must contain a column for the Individual Identifier (IID) to facilitate sample matching, followed by a column representing the continuous phenotypic values (such as HDL levels or other biomarkers).

3. **Model Covariates**

**Input:** A multi-column tab-delimited text file.

Description: This file contains the independent variables to be included in the multivariate regression model. After the IID column, it should include metadata, such as biological sex, age, or principal components (PCs) required to adjust for population stratification and other confounding factors.

4. **Association Statistics Output**

**Input:** A filename for the tab-delimited results.

Description: The script must generate a summary table where each row represents a tested variant. The output must explicitly include:

Variant Metadata: Genomic coordinates (CHROM and POSITION) and the specific alleles (REF and ALT).

Sample Size: The total number of individuals (N) passing quality control for that specific test.

Statistical Metrics: The regression coefficient (EFFECT), representing the slope of the association, and the corresponding P-VALUE, indicating the statistical significance of the variant-phenotype correlation.
