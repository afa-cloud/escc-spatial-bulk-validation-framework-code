# S2 Code

This archive contains the reproducible scripts and configuration files used for the public-data ESCC spatial-to-bulk validation workflow.

## Environment

- Python 3.12.13 was used for the public-data analysis scripts in the local environment.
- The submitted S2 Code package is Python-based. R was not used for the reported analyses in this code package.

## Main scripts

- `scripts/run_spatial_axis_deep_validation.py`
- `scripts/run_independent_patient_and_spatial_quant.py`
- `scripts/run_transferability_supplement.py`
- `scripts/assemble_s1_table.py`
- `scripts/public_data_helpers.py`

## Data

All primary datasets are public: TCGA/GDC/UCSC Xena, GEO GSE47404, GEO GSE53625, HRA003627, HRA008846 and GDSC2. This code archive does not contain large source datasets.

This archive contains analysis code and configuration files only. Administrative submission files and local utility scripts are not included. The S1 assembly helper can combine generated TSV outputs into the submitted S1 workbook layout when the generated source tables are present locally.

## Archive

This repository uses Zenodo GitHub integration for release archiving. The repository concept DOI is https://doi.org/10.5281/zenodo.19826728; the release-specific version DOI is shown on the Zenodo record for the archived GitHub release.
