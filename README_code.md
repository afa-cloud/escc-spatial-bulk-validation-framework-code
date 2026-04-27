# Reproducibility Code

These scripts reproduce the public-data analysis tables used by the manuscript. The public code package intentionally excludes manuscript drafting, DOCX generation, cover-letter generation, and internal audit-rendering scripts.

- run_spatial_axis_deep_validation.py
- run_independent_patient_and_spatial_quant.py
- run_transferability_supplement.py
- public_data_helpers.py
- project_config.yaml
- requirements.txt

Run order: deep validation, independent patient/spatial quantification, then transferability supplement.
The workflow uses public data only and records executor/reviewer gates in output tables.
