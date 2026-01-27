# CD20-Negativity in Relapsed B-Cell Lymphomas

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code for a nationwide population-based study investigating **CD20 expression at relapse in patients with B-cell lymphomas** treated with rituximab-containing first-line therapy.

The study examines the prognostic impact of CD20-negativity at relapse in three major B-cell lymphoma subtypes:
- **DLBCL** (Diffuse Large B-Cell Lymphoma)
- **FL** (Follicular Lymphoma)
- **MCL** (Mantle Cell Lymphoma)

## Study Design

### Cohort Selection
Patients were identified from the **Danish National Lymphoma Registry (RKKP_LYFO)** with the following inclusion criteria:
- Diagnosis of DLBCL, FL, or MCL
- Received rituximab as part of first-line immunochemotherapy
- Relapsed and received second-line treatment

### CD20 Status Determination
CD20 expression status was determined through:
1. **Text mining** of pathology reports from the Danish Pathology Data Bank (SDS_pato)
2. **Manual review** of pathology reports within a defined time window (60 days before to 30 days after second-line treatment initiation)

### Exclusions
- Primary CNS lymphoma (PCNSL)
- Plasmablastic lymphoma (PBL)
- CNS-directed first-line therapy

## Repository Structure

```
CD20-negativity_in_relapse_B-cell_lymphomas/
├── 1_CD20_textmining_SDSpato.R         # Text mining of pathology reports
├── 2_defining_studypopulation_CONSORT_CD20.R  # Study population definition
├── 3_Cleaning_cohort.R                 # Cohort cleaning and exclusions
├── 4_1st_line_categorization.R         # First-line treatment categorization
├── 5_ASCT_polyrx.R                     # ASCT and polypharmacy calculations
├── 6_CD20_cohort_baseline_charac.R     # Baseline characteristics compilation
├── 7_Table1_Table2.R                   # Descriptive tables generation
├── 8_all_survivalanalysis_and_sensitivity.R  # Survival analysis (full)
├── 8_survivalanalysis_and_sensitivity.R      # Survival analysis (alternative)
├── 9_2L_categorization_table3.R        # Second-line treatment categorization
├── 9_masterscript_.R                   # Master script to run full pipeline
└── README.md                           # This documentation
```

## Data Sources

This study utilizes data from the **DALYCARE** research database hosted at the **Danish National Genome Center (NGC)**. The following data sources are integrated:

| Data Source | Description |
|-------------|-------------|
| **RKKP_LYFO** | Danish National Lymphoma Registry |
| **SDS_pato** | Danish Pathology Data Bank |
| **SDS_t_mikro_ny** | Pathology report free-text |
| **SDS_epikur/ekokur** | National Prescription Registry |
| **patient** | Patient demographics |
| **IPI** | International Prognostic Index scores |

## Analysis Pipeline

### 1. Text Mining (`1_CD20_textmining_SDSpato.R`)
- Loads pathology reports for the rituximab-treated lymphoma cohort
- Filters sentences containing "CD20" mentions
- Prepares data for manual review

### 2. Study Population Definition (`2_defining_studypopulation_CONSORT_CD20.R`)
- Applies inclusion/exclusion criteria
- Links pathology reports to treatment dates
- Filters for reports within the relapse time window

### 3. Cohort Cleaning (`3_Cleaning_cohort.R`)
- Excludes PBL and PCNSL patients
- Categorizes CD20 status (Positive, Reduced, Negative)
- Handles ASCT patients separately

### 4. First-Line Treatment (`4_1st_line_categorization.R`)
Categorizes first-line regimens into:
- CHOP-like (including CHOEP, EPOCH-R, Bio-CHIC)
- Bendamustine/BAC
- R-maxi-CHOP/TRIANGLE (for MCL)
- R-monotherapy
- CVP
- Other

### 5. Polypharmacy Assessment (`5_ASCT_polyrx.R`)
- Calculates number of unique ATC codes in the year before diagnosis
- Defines polypharmacy as ≥5 medications

### 6. Baseline Characteristics (`6_CD20_cohort_baseline_charac.R`)
Compiles clinical variables:
- Age at diagnosis and second-line
- Sex, ECOG performance status
- Ann Arbor stage, IPI score
- LDH levels, extranodal sites
- Response to first-line therapy
- Time to second-line treatment (TT2L)

### 7. Tables (`7_Table1_Table2.R`)
Generates:
- **Table 1**: Baseline characteristics by subtype and CD20 status
- **Table 2**: Response rates and time to relapse

### 8. Survival Analysis (`8_all_survivalanalysis_and_sensitivity.R`)
- **Primary analysis**: Cox regression adjusted for age, sex, and IPI
- **Kaplan-Meier curves**: Overall survival from second-line treatment
- **Sensitivity analyses**:
  - Adjusting for response to first-line therapy
  - Adjusting for early vs. late relapse (POD12/POD24)
  - Excluding secondary CNS relapse
  - Three-category CD20 analysis (Positive/Reduced/Negative)

### 9. Second-Line Treatment (`9_2L_categorization_table3.R`)
Categorizes second-line regimens:
- Salvage with/without HDT (high-dose therapy)
- CHOP, R-CVP, BR
- Targeted therapies (±anti-CD20)
- Bispecific antibodies (BsAbs)
- CNS-directed therapy
- Palliative

## Key Variables

### CD20 Status Categories
| Code | Description |
|------|-------------|
| E | Expressed (Positive) |
| L | Low/Reduced expression |
| N | Negative |

### Time to Second-Line Groups
| Subtype | Early Relapse | Late Relapse |
|---------|---------------|--------------|
| DLBCL | POD12 (≤12 months) | >12 months |
| FL | POD24 (≤24 months) | >24 months |
| MCL | POD12 (≤12 months) | >12 months |

## Dependencies

### R Packages
```r
# Data manipulation
library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(arrow)
library(lubridate)

# Survival analysis
library(survival)
library(survminer)
library(ggsurvfit)

# Tables and visualization
library(table1)
library(tableone)
library(Publish)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(broom)
```

### Custom DALYCARE Functions
The analysis relies on the DALYCARE package:
```r
source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
```

Key functions:
- `load_dataset()` - Load data tables
- `clean_RKKP_LYFO()` - Clean lymphoma registry data
- `clean_SDS_t_mikro()` - Clean pathology text
- `filter_sentence()` - Text mining helper
- `diff_days()` / `diff_years()` - Date calculations
- `nrow_npatients()` - Count rows and unique patients
- `KM_plot()` - Kaplan-Meier visualization

## Running the Analysis

### Prerequisites
1. Access to the NGC secure computing environment
2. DALYCARE database credentials
3. R installation with required packages

### Execution
Run the master script to execute the full pipeline:
```r
source('9_masterscript_.R')
```

Or run individual scripts in numerical order:
```r
source('1_CD20_textmining_SDSpato.R')
source('2_defining_studypopulation_CONSORT_CD20.R')
# ... continue through script 9
```

### Output
The analysis generates:
- Intermediate data files in `/ngc/projects2/dalyca_r/sangho_r/data/`
- Figures in `/ngc/projects2/dalyca_r/sangho_r/figures/`
- Tables exported via `fwrite()` and `write_utable()`

## Statistical Methods

### Primary Analysis
- **Multivariable Cox proportional hazards regression**
  - Outcome: Overall survival from second-line treatment
  - Exposure: CD20 status (Positive vs. Negative)
  - Covariates: Age, Sex, IPI score

### Sensitivity Analyses
1. Additional adjustment for response to first-line therapy
2. Additional adjustment for time to relapse (early vs. late)
3. Exclusion of patients with secondary CNS relapse
4. Three-category CD20 analysis with pairwise log-rank tests

### Visualization
- Forest plots with hazard ratios (custom `ggforestfixed()` function)
- Kaplan-Meier survival curves with risk tables

## Authors

- **Sanaz M. Gholy** - Primary analyst
- **Jojo Biel-Nielsen Dietz** - Co-analyst
- **Christian Brieghel** - Senior author and methodology

## Affiliation

**RH-CLL-LAB**  
Department of Hematology  
Rigshospitalet, Copenhagen University Hospital  
Denmark

## Data Availability

This study uses Danish national health registry data, which cannot be shared publicly due to data protection regulations. Access to similar data can be requested through the Danish Health Data Authority and Statistics Denmark.

## Ethics

The study was conducted in accordance with Danish data protection regulations and approved by relevant authorities. Individual patient consent was not required for this registry-based study.

## Citation

If you use this code, please cite the associated publication (details to be added upon publication).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- DALYCARE Research Database team
- Danish National Lymphoma Registry (LYFO/RKKP)
- Nationalt Genom Center for computational resources

