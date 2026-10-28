# Input/output manifest

Maps each script in `Analyses_codes/` to its required inputs, outputs, and the manuscript figures and tables it supports.
Draft written against commit `65fa129`.

The repository contains analysis code, not a self-contained end-to-end workflow. Individual-level data are controlled-access, so each script starts from objects or files that must be prepared beforehand. The "Prepared externally" rows below say what those are.

## Workflow overview

```
cohort-level XWAS (BOLT-LMM / Regenie / PLINK2; not in repo)
        |   per-cohort summary statistics
        v
01A  rXCI, Female, Male meta-analysis (GWAMA, fixed effects) --+
        |  Female + Male results                                |
        v                                                       |
01B  eXCI meta-analysis (GWAMA, random effects)                 |
        |                                                       |
        +--------------------+----------------------------------+
                             v
                  01C  Circos plot
                             |
   XWAS + trait/xQTL summary statistics, harmonized into `merged` (not in repo)
                             v
        02A  coloc.abf     02B  coloc.susie
                             |
   prioritized genes + AMP-AD DEG results (not in repo)
                             v
                    03  bulk-brain DEG enrichment
```

## 01A_GWAMA_master_non-eXCI.R

Fixed-effects meta-analysis of cohort summary statistics for the rXCI, female-only and male-only XWAS.

**Prepared externally (must exist before running)**

| Object | Description |
|---|---|
| `adgc`, `adsp`, `ukb`, `fg`, `mvp1`, `mvp2` | Cohort summary-statistic tables (data.tables), one per cohort, with the 14 columns listed below. |
| `filepath` | Character vector of paths to the per-cohort summary-statistic files that GWAMA reads. |

**Cohort input file format (one row per variant, chromosome X coded as 23)**

| Column | Type | Description |
|---|---|---|
| `CHR` | integer | Chromosome (23 = X). |
| `POS` | integer | Base-pair position. |
| `MARKER` | character | Variant ID in the form `23:POS:NEA:EA`, for example `23:2781604:G:T`. |
| `SNP` | character | rsID. |
| `NEA` | character | Non-effect allele. |
| `EA` | character | Effect allele. |
| `N` | integer | Sample size for the variant. |
| `OR` | numeric | Odds ratio for the effect allele. |
| `SE` | numeric | Standard error. |
| `OR_95L`, `OR_95U` | numeric | 95% confidence limits of the odds ratio. |
| `Z_stat` | numeric | Z statistic. |
| `P` | numeric | P-value. |
| `EAF` | numeric | Effect-allele frequency in that cohort. |

`EMP` is a placeholder column that is not in the cohort files. `01A` adds it when each cohort file is read (`fread(filepath, header = T) %>% mutate(EMP = 0)`) because the GWAMA map file it builds has the columns `CHR`, `MARKER`, `EMP`, `POS`. `01B` builds its map file without it (`CHR`, `MARKER`, `POS`).

**Outputs**

| File | Contents |
|---|---|
| `gwama.map.txt` | Tab-delimited map file: `CHR`, `MARKER`, `EMP`, `POS` (`EMP` is the placeholder 0 described above). No header. |
| `gwama.key.txt` | Key file: `MARKER`, `SNP`. Used to restore rsIDs after GWAMA. |
| `gwama.merge_list.txt` | One input-file path per line. |
| GWAMA output `.out` | Raw GWAMA result. |
| `stats.gen090` | Filtered meta-analysis, written to `out_dir` under this fixed file name, columns: `CHR`, `BP`, `MARKER`, `ALLELE1`, `ALLELE0`, `A1FREQ`, `BETA`, `SE`, `BETA_95L`, `BETA_95U`, `Z`, `P`, `P_neg_log10`, `Q_V`, `Q_P`, `i2`, `N_studies`, `N_incl`, `META_DIR`, `SNP`. |

The output file name is fixed, so each run (rXCI, female, male, and each APOE4 stratum) needs its own `out_dir` or its own file name. Otherwise the files overwrite each other. The female and male versions of this file are the inputs to `01B`.

**Filters:** drops variants whose first `META_DIR` character is `?` and keeps variants with `N_incl` above 90% of the maximum `N_incl`.

**Supports:** rXCI, female and male XWAS results (Fig 2, Table 1, Table S4, Fig S1).

## 01B_GWAMA_master_eXCI.R

eXCI model. Random-effects meta-analysis of the female and male XWAS, after rescaling male estimates to the 0/1 encoding.

**Inputs**

| File | Contents |
|---|---|
| Female summary statistics (`stats.gen090`) | Output of 01A for the female XWAS. Needs `MARKER`, `SNP`, `CHR`, `POS`, `BETA`, `SE`. |
| Male summary statistics (`stats.gen090`) | Output of 01A for the male XWAS (0/2 encoding). |

Male `BETA` and `SE` are multiplied by 2 to convert from 0/2 to 0/1 encoding. If males were already coded 0/1, skip that step. The 01A outputs carry `MARKER`, which 01B reuses as the GWAMA marker name.

**Outputs:** `XWAS_female.txt`, `XWAS_male.txt` (male rescaled), `gwama.map.txt` (`CHR`, `MARKER`, `POS`), `gwama.key.txt`, `gwama_list.txt`, the GWAMA `.out` file, and the filtered eXCI results (same columns and 90% filter as 01A).

**Supports:** eXCI XWAS (Fig 2, Table 1, Table S4).

## 01C_XWAS_circos_plot.R

**Inputs:** four summary-statistic files (rXCI, eXCI, female, male; hg38) with at least `SNP`, `BP`, `P`.
**Output:** `plot_name.jpg` in the working directory (11 x 11 in, 320 dpi), one circular Manhattan plot per APOE4 stratum.
**Supports:** Fig 2B.

## 02A_XWAS_coloc_abf.R and 02B_XWAS_coloc_susie.R

Run for both trait pleiotropy and xQTL colocalization. One call = one locus and one paired dataset.

**Prepared externally**

| Object | Description |
|---|---|
| `merged` | AD XWAS and paired trait or QTL data, harmonized to the same alleles over one locus (+/-1 Mb). Columns: `SNP`, `BP`, `AD_P`, `AD_N`, `s`, `AD_BETA`, `varbeta_ad`, `AD_MAF`, `QTL_P`, `QTL_N`, `QTL_BETA`, `varbeta_t2`, `QTL_MAF`. `s` is the case fraction of the AD dataset (number of cases / total sample size) and is repeated on every row. `02A` passes its first value (`merged$s[1]`) to coloc, which expects a single number. |
| `snps`, `LD` (02B only) | Variant IDs and the LD matrix for the locus, with `SNP` as row and column names. Computed from an in-house TOPMed-imputed European LD panel. |

`AD_MAF` is the minor-allele frequency from the AD XWAS. When an external dataset is merged with the AD XWAS, its alleles are aligned to the AD MAF, so `QTL_MAF` is on the same basis. This happens in the step that builds `merged`, upstream of these scripts.

**Procedure (02B):** `runsusie` is run on both datasets at coverage 0.95, then 0.85, 0.75, ... down to 0.05 (steps of 0.1) until both datasets give at least one credible set. `coloc.susie` is then run, and the row with the highest PP4 is kept.

**Outputs**

| Output | Contents |
|---|---|
| `coloc_results_default` (02A) | `coloc.abf` result object: PP.H0 to PP.H4. |
| `out$summary`, `coloc_PP0` to `coloc_PP4`, `coloc_hit1`, `coloc_hit2` (02B) | `coloc.susie` summary and the best-PP4 row. |
| LocusCompare plot | AD against QTL or trait, anchored on the lead AD variant (hg38). |

The scripts print results; they do not write them to disk. Results were compiled outside the scripts into one row per locus and paired dataset, with these columns:

| Column | Description |
|---|---|
| `rsID`, `BP` | Lead variant and its base-pair position. |
| `t1_dset` | AD XWAS dataset (model and APOE4 stratum). |
| `t2_dset` | Paired trait or QTL dataset. |
| `t2_group` | Group label of the paired dataset (for example `Pathology`). |
| `gene_name` | Gene tested (xQTL analyses). |
| `PP0` to `PP4` | Posterior probabilities for hypotheses H0 to H4 (`PP.H0.abf` to `PP.H4.abf`). |

These rows are the source of the colocalization supplementary tables. Both the ABF (`02A`) and SuSiE (`02B`) results are compiled with these same columns.

**Supports:** pleiotropy results (Fig 3B, Table S8, Table S11, Fig 4B), xQTL colocalization (Fig 5A, Fig S4, Table S13), and the LocusCompare plots on Zenodo (10.5281/zenodo.18807726). The Fig 6 and Fig S5 to S9 panels draw on the same results.

## 03_bulk_DEG.R

Tests whether the prioritized X-linked genes show larger absolute differential-expression t-statistics than random X-chromosome genes.

**Inputs**

| File | Columns used |
|---|---|
| `gene_reference.csv` | `gene_name`, `ens_id`. This is the list of prioritized genes. |
| `DEG_results.csv` | `chromosome_name`, `hgnc_symbol`, `ensembl_gene_id`, `Comparison`, `t`. AMP-AD DEG results (Synapse syn26720676). |

The script keeps X-chromosome rows and one `Comparison` contrast. The contrast written in the script, `diagnosis_tissue*age_death.AD_CBE - diagnosis_tissue*age_death.CT_CBE`, is an example from the cerebellum (CBE). The same script was applied to each of the eight brain regions by changing the `Comparison` string to that region's AD-versus-control contrast.

**Procedure:** 1,000 random draws (seed 123) of as many X genes as were prioritized, then a Z-test of the prioritized mean against the simulated mean.
**Outputs:** printed `z` and `p`; histogram of the null and violin plot, combined with patchwork (not saved to file by the script).
**Supports:** Fig 5B, which has one panel per brain region (eight in total); only the ACC enrichment is significant (P = 0.03).

## External data and software releases

| Item | Version or release |
|---|---|
| GWAMA | v1.2.6 |
| coloc | v4.2.1 |
| CMplot | v4.5.1 |
| locuscompare | v1.0.0 |
| FinnGen | R12 |
| VEP | release 115 |
| BOLT-LMM | v2.5 |
| GWAS/XWAS traits | See Table S5 (566 traits) |
| xQTL resources | See Table S12 |
| Genome build | hg38 |
