# Reproducibility Code

These scripts reproduce the public-data analysis tables used by the manuscript. The public code package contains analysis code only; administrative submission files and local quality-control utilities are not included.

- run_spatial_axis_deep_validation.py
- run_independent_patient_and_spatial_quant.py
- run_transferability_supplement.py
- assemble_s1_table.py
- public_data_helpers.py
- project_config.yaml
- requirements.txt

Run order: deep validation, independent patient/spatial quantification, transferability supplement, then optional S1 workbook assembly when generated TSV outputs are available.
The workflow uses public data only and records quality-control checks in output tables.
