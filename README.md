# 2025-nestedness-project
Project on nestedness in phonological inventories

## Pre-registration

The pre-registration can be accessed on OSF [here](https://osf.io/sfc6p/).

## Reproduction

### Folders

- **Figures**: contains all figures used in the pre-registration or shared on Workflowy
- **Results Datasets**: contains all datasets obtained after running the analysis scripts (old NODF, new NODF, temperature, consonant/vowel inventories)
- **Analysis**: contains the scripts for the different analyses

### Data

- **Phoible** data were used for all main analyses. The dataset can be downloaded [here](https://zenodo.org/records/2677911)
- **Grambank** and **Glottolog** data were used for preliminary analyses. Grambank dataset can be dowloaded [here](https://github.com/annagraff/crossling-curated). This curated version originates from the article [“Curating global datasets of structural linguistic features for independence”](https://www.nature.com/articles/s41597-024-04319-4) by Anna Graff et al. It is specifically designed to minimize logical and statistical dependencies among Grambank features. Glottolog dataset can be downloaded [here](https://glottolog.org/meta/downloads).

### Analysis

<b><u>Phoible</u></b>

- **analysis.R**: performs the nestedness analysis of phonological inventories grouped by language family. It includes a section for matrix visualization and a preliminary exploration of borrowings using SegBo.
- **analysis_all_lang.R**: conducts the same nestedness analysis, this time including all languages in Phoible without grouping by family.
- **analysis_cons_vow.R**: repeats the nestedness analysis by language family, but focuses separately on either vowel inventories or consonant inventories only.


<b><u>Grambank</u></b>

- **analysis_grambank.R**: visualizes the distribution of Grambank data through various plots, and merges it with Glottolog data to retrieve language family names.
