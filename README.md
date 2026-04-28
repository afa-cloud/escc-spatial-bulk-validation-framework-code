# ESCC spatial-to-bulk validation framework

This repository contains the analysis code and reproducibility materials for the manuscript:

**A tiered spatial-to-bulk reproducibility framework supports an association-level CAF/ECM stromal-remodeling phenotype in esophageal squamous cell carcinoma**

The study is a public-data bioinformatics workflow. It applies a tiered validation framework that connects:

1. spatially nominated hypotheses,
2. bulk transcriptomic validation, and
3. source-table checks from published supplementary quantitative tables.

The repository is intended to support code availability and reproducibility for the public-data analyses. It does not contain controlled-access data, newly generated human-subject data, manuscript drafts, administrative submission files, full-text PDFs, or temporary local files.

## Repository contents

- `scripts/`: Python scripts used for the public-data workflow, independent patient/source-table checks, supplemental transferability analysis, and S1 workbook assembly from generated TSV outputs.
- `project_config.yaml`: project-level configuration and declared study scope.
- `requirements.txt`: Python dependencies used by the submitted code package.
- `manifests/`: public-data source manifest.
- `supporting_information/`: submitted processed tables copied from the manuscript package.
- `reproducibility_check.tsv`: local check summary.
- `S2_Code_manifest.json`: manifest for the submitted code package.
- `.zenodo.json` and `CITATION.cff`: metadata used for DOI archiving.

## Data sources

The workflow uses public resources only, including TCGA/GDC/UCSC Xena, GEO cohorts GSE47404 and GSE53625, HRA003627 source data, HRA008846 supplementary tables, and GDSC2. Dataset accessions, URLs, dates, checksums where available, and processed-output paths are listed in the manifest files and in the manuscript supporting information.

## Reproducibility scope

This repository is designed to reproduce the reported public-data analyses and submitted tables from the available public data layers. Some source tables are redistributed as processed supporting information in the manuscript package; raw third-party datasets remain available from their original repositories. The workflow does not perform causal inference, wet-lab validation, prospective clinical validation, or expression-drug response modeling.

## Quick start

Use Python 3.12 or a compatible Python 3 environment.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Run scripts from the repository root after reviewing `project_config.yaml` and the source manifest. Output paths are configured inside the scripts and may need to be adjusted if the repository is moved.

## Citation

Cite the clean-text public code release for the manuscript-submission version:

https://github.com/afa-cloud/escc-spatial-bulk-validation-framework-code

Release tag: `v1.0.10-submission-consistency`.

The GitHub release page records the exact release commit.

Zenodo concept DOI for this repository: https://doi.org/10.5281/zenodo.19826728.

The code-availability statement for peer review should cite the archived release tag, its Zenodo version DOI, and the submitted S2 Code package.

## License

Code in this repository is released under the MIT License. Third-party datasets remain governed by the terms of their original repositories and publications.
