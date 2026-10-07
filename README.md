# PsychAD snMultiome — manuscript analysis code

## Code organization

| Directory | Intended contents |
|---|---|
| `01_preprocessing/` | QC, metadata harmonization, paired/single-modality integration, cell annotation |
| `02_disease_scores/` | Dx_donor construction; Dx_cell features, training, cross-validation, saved-model use |
| `03_cellular_responses/` | Dx_dev computation; pathology/covariate adjustment; composition/cognition models |
| `04_state_resolved_analysis/` | Within-sample cellular contrasts (Dx_contrast), state assignment, pseudobulk and differential models |
| `05_regulatory_programs/` | Chromatin and RNA regulatory analysis, TF/motif associations and eRegulons |
| `06_donor_factors/` | Donor-associated axes, factor loadings, pathway-level interpretation |
| `07_external_validation/` | Analyses in independent datasets and projection of learned scores or axes |

## Data, and frozen models

| Resource | Contents | Availability |
|---|---|---|
| [PsychAD study on Synapse (syn52160016)](https://www.synapse.org/Synapse:syn52160016) | Raw and processed multiome data, cell- and sample-level metadata, and intended complete DEG/DAC analysis outputs | **Submission in progress.** This link currently identifies the broader MSSM_PsychAD study folder; file-level accessions and permissions for this dataset should be confirmed after deposition. |
| [UCSC Cell Browser: PsychAD snMultiome](https://cells.ucsc.edu/?ds=psychad-snmultiome) | Interactive visualization of processed single-nucleus data, annotations, and available score fields | Browser link supplied by the study authors. Confirm public visibility and displayed content before release. |
| [Zenodo: frozen Dx_cell models](https://zenodo.org/records/23190790) | Fixed model artifacts for reproducing Dx_cell predictions | Record link supplied by the study authors. **Verify record availability, artifact checksums, and correspondence to the manuscript models** before release. |

## citation
Please cite the following:
```
[TBD]
```
## contact
For questions regarding the analysis code, processed data, or computational methods, please contact:
Pengfei Dong
Center for Disease Neurogenomics
Icahn School of Medicine at Mount Sinai
Email: pengfei.dong@mssm.edu

