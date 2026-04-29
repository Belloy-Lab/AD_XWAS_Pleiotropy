**An APOE*4-Informed Genomic Atlas of the X Chromosome in Alzheimer’s Disease**

## Background
To elucidate sex differences in Alzheimer’s disease (AD), we present the most comprehensive analysis to date of X chromosome genetic variation in Alzheimer’s disease (AD), integrating:
- X-chromosome-wide association studies (XWAS)
- Sex-stratified and APOE*4-stratified analyses
- Modeling of escape from X-chromosome inactivation (eXCI)
- Multi-trait pleiotropy with brain imaging, lipid, and hormone traits
- Functional genomics (xQTL colocalization and differential expression)
This work provides a resource of X-linked AD risk loci, candidate genes, and biological insights to guide future research.

---

## Repository Structure
```
├── Analysis_codes/
│ ├── 01A_GWAMA_master_non-eXCI.R # XWAS Meta-analyses: rXCI, Female, and Male Modeling (GWAMA)
│ ├── 01B_GWAMA_master_eXCI.R # XWAS Meta-analyses: eXCI Modeling (GWAMA)
│ ├── 01C_XWAS_circos_plot.R # XWAS: Circos Plot Creation
│ ├── 02A_XWAS_coloc_abf.R # Trait Pleiotropy and xQTL COLOC example script (ABF)
│ ├── 02B_XWAS_coloc_susie.R # Trait Pleiotropy and xQTL COLOC example script (SuSiE)
│ └── 03_bulk_DEG.R # Bulk brain tissue DEG analyses
└── README.md
└── LICENSE
```

---

## Data availability

### XWAS Analyses
Meta-analyzed XWAS summary statistics will be made available upon publication in GWAS catalog and Zenodo. Individual-level genetic data used in the X-chromosome-wide association analyses are not publicly hosted in this repository due to controlled-access restrictions. Data are available upon application through the following repositories:

| Datasource | Link                                          |
|------------|-----------------------------------------------|
| dbGaP      | https://www.ncbi.nlm.nih.gov/gap/             |
| NIAGADS    | https://www.niagads.org/                      |
| LONI       | https://ida.loni.usc.edu/                     |
| AMP-AD Knowledge Portal (Synapse) | https://www.synapse.org/ |
| Rush ADRC  | https://www.radc.rush.edu/                    |
| NACC       | https://naccdata.org/                         |
| UK Biobank | https://www.ukbiobank.ac.uk/                  |
| FinnGen    | https://www.finngen.fi/en                     |
| Million Veteran Program (MVP)    | https://www.mvp.va.gov/ |

---

### Multi-Trait Pleiotropy Analyses
All external datasets used for pleiotropy analyses are publicly available:

| Trait      | Link                                          |
|------------|-----------------------------------------------|
| Age at menopause      | https://www.reprogen.org/data_download.html               |
| Age at menarche    | https://www.repository.cam.ac.uk/items/8c5f7afb-5fa2-45ea-b52d-4e643bc2a5b7    |
| Hormone exposure traits (UK Biobank)   | https://www.ukbiobank.ac.uk/                     |
| Global Lipids Genetics Consortium (GLGC) | https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/chrx_summary_stats/  |
| Brain MRI XWAS (Jiang et al.)  | https://zenodo.org/records/12676622                    |

---

### xQTL Colocalization Analyses
Publicly available xQTL resources include:
 
| Datasource | Link                                          |
|------------|-----------------------------------------------|
| Wingo et al. DLPFC pQTL/eQTL                         | https://www.synapse.org/Synapse:syn51150434/wiki/621280 |
| Fujita et al. single-cell brain eQTL data            | https://www.synapse.org/Synapse:syn52335807  |
| eQTL Catalogue (GTEx v8, CommonMind, et.c eQTL/sQTL) | https://www.ebi.ac.uk/eqtl/Data_access/      |

---

### Differential Expression Analyses
Bulk brain RNA-seq datasets used for differential expression analyses are available upon application via:

- AMP-AD Knowledge Portal:  
  https://www.synapse.org/Synapse:syn26720676

---

## Citation 
If you use these scripts, please cite our paper:\
[insert pre-print/citation once uploaded]

## Contact
Noah Cook\
noahc@wustl.edu

## License (MIT)
Copyright (c) 2026 Michael Belloy
