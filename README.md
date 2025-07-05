# GWAS_AssocTest
This repository contains code and notes for running a GWAS (See GWAS_Phenotype and GWAS_QC for steps prior)

## Association tests

- The GWAS was run using **PLINK 2.0** on SwissArmyknife because the PLINK GWAS applet on DNANEXUS only accepts BGEN files (which limits QC options) or BIM/BED/FAM formats, making it unsuitable for this dataset.
- UKB partition used for this was `mem1_hdd1_v2_x16`

### Command used:

```bash
plink2 --pfile 6.plinkmaf_c22_snps --pheno ukb_phenotype_data.pheno --covar ukb_covariates.cov --glm --no-input-missing-phenotype --out final_gwas_results
```

- The `--no-input-missing-phenotype` flag was included to prevent PLINK from treating -9 values in the phenotype file as missing data.
- This outputs:
  1) `final_gwas_results.log`
  2) `final_gwas_results.MetaboHealth_Score.glm.linear`

## Visualisation tests
# QQ plot 
- run `qq_plot.R`

## Manhattan plot
FOR BONFERRONI ADJUSTED
- run `manhattan_plot.R`

FOR Benjamini-Hochberg/FDR ADJUSTED 
- run `manhattan_plot_FDR.r`

### If Using Benjamini-Hochberg/FDR
- go to LDmatrix
- input the `.txt` file created by `manhattan_plot_FDR.r`
- set it to `GRCH38 High Coverage`
- get D and r2 stats - saved here as `d_ldsnps.txt` and `r2_inlds.txt`
