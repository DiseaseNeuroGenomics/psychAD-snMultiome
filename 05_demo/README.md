# SEA-AD MTG demonstration: Dx_cell projection, Dx_dev and neuronal composition

This example applies the **frozen PsychAD-snMultiome Dx_cell models** to the independent **SEA-AD middle temporal gyrus (MTG) snRNA-seq** dataset. It then derives **SEA-AD-specific Dx_dev** and tests associations between neuronal subtype composition, neuropathology (CPS) and glia-derived Dx_dev. **The pretrained models are not retrained.**

**Scope:** This is a **full-dataset, HPC-scale reproduction example**, not a small example validated on an ordinary desktop. The SEA-AD MTG input file is large (34Gb). The recorded ~1-hour / 329.6-GB peak benchmark predates the addition of the QC export and R analysis. 

### 1. Directory layout

```text
PsychAD-snMultiome/
├── python.yml                  # Python analysis environment
└── 05_demo/
    ├── 01_Dx_cell_projection.py
    ├── 02_Dx_dev_crumblr.R
    ├── input/                 # Downloaded SEA-AD h5ad (not tracked by Git)
    ├── models/                # Extracted Zenodo model archive
    └── output/                # Predicted scores, QC and association results
```

### 2. Software and environment

The repository-level [`python.yml`](../python.yml) specifies the original Python environment, including Python 3.8.15, Scanpy 1.9.3, AnnData 0.9.2, pandas 2.0.3, NumPy 1.22.4 and LightGBM 3.3.5. To install from the `05_demo/` directory:

```bash
conda env create -f ../python.yml -n psychad-demo
conda activate psychad-demo
python -c 'import scanpy, pandas, numpy, lightgbm; print("Python dependencies loaded")'
```

The R analysis requires `readxl`, `data.table`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `ggrepel`, `crumblr` and `variancePartition` (with their dependencies). The reported R library versions are listed below. 

Tested R package versions used in the analysis: `R 4.4.1`, `readxl 1.4.3`, `data.table 1.16.2`, `variancePartition 1.36.3`, `BiocParallel 1.40.0`, `limma 3.62.2`, `rhdf5 2.50.0`, `ggrepel 0.9.6`, `crumblr 0.99.19`, `purrr 1.0.2`, `tidyr 1.3.1`, `dplyr 1.1.4`, and `ggplot2 4.0.3`. 

### 3. Download SEA-AD data and frozen models

The workflow uses the **archived February 13, 2024 SEA-AD MTG release**.

```bash
mkdir -p input models output

wget -O input/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/previous_objects/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad

wget -O model.zip \
  https://zenodo.org/records/23190790/files/model.zip
unzip -q model.zip -d models

# Published SEA-AD MTG scCODA result table (needed for Step 2)
wget -O input/SEAAD_Supplementary_Table_4.xlsx \
  "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/Supplementary+Information/Supplementary+Table+4.xlsx"
```

The Zenodo ZIP includes both the pretrained models and each class-specific `ref_avg_expression.csv`. 


### 4. Step 1 — project Dx_cell and export nucleus-level QC

```bash
python 01_Dx_cell_projection.py \
  --input input/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  --model_dir models \
  --output_path output/SEAAD_MTG_Dx_cell.csv.gz \
  --obs_output output/SEAAD_MTG_obs_QC.csv.gz
```


**Outputs:**

- `output/SEAAD_MTG_Dx_cell.csv.gz`: original metadata plus mapped `class`, predicted `Dx_cat_cell` (`control`, `early`, `late`) and continuous `Dx_cell`, for six major cell classes
- `output/SEAAD_MTG_obs_QC.csv.gz`: all original nucleus metadata plus Scanpy QC metrics such as `n_genes_by_counts`, `pct_counts_mt` and `pct_counts_ribo`; this file covers the original AnnData observations

Both CSV files retain nucleus identifiers as the first column. **Do not change those identifiers**: the R analysis joins score rows back to the QC table by row name.

### 5. Step 2 — estimate SEA-AD Dx_dev and test neuronal composition

From `05_demo/`, run the following after Step 1 has completed:

```bash
cd output
Rscript ../02_Dx_dev_crumblr.R \
  SEAAD_MTG_Dx_cell.csv.gz \
  SEAAD_MTG_obs_QC.csv.gz \
  ../input/SEAAD_Supplementary_Table_4.xlsx
cd ..
```

The R script expects **three positional command-line arguments**, in this order: the Dx_cell prediction CSV, the nucleus-level observation/QC CSV, and the downloaded SEA-AD scCODA Excel workbook (`MTG_RNAseq Neu CPS` sheet). Its outputs are written to the **current working directory** (`output/` in this example).

**What the R script computes:**

1. Aggregates per-nucleus `Dx_cell` to the **donor × major cell-class** mean for MTG.
2. Computes a class-specific Dx_dev residual from a fixed-effects model of mean Dx_cell on **SEA-AD CPS**, age at death, sex, PMI, and the mean RNA QC covariates (`n_genes_by_counts`, `pct_counts_mt`, `pct_counts_ribo`).
3. Standardizes the four glial-class Dx_dev estimates (`Astro`, `Oligo`, `OPC`, `Micro_PVM`) and takes their mean for each donor.
4. Constructs neuronal subtype counts from nuclei annotated as `EN` or `IN`, using the source `Supertype` labels. It uses **crumblr** and **dream** to test neuronal composition associations with CPS and with glia-derived Dx_dev after adjustment for CPS, age, sex and PMI.
5. Loads the published MTG neuronal CPS scCODA results, averages the thresholded `Final Parameter` and posterior inclusion probability across reference-cell models for each neuronal supertype, and compares the resulting exploratory scCODA CPS effect with the crumblr CPS effect.
6. Plots subclass-level estimates from the two crumblr analyses and the scCODA–crumblr CPS comparison. The variable labeled `Dx_donor` internally in the R plotting data actually represents the **SEA-AD CPS association coefficient**, not the PsychAD Dx_donor score.

The analysis is limited to **neuronal subtype composition** because the SEA-AD FACS sampling strategy influences broad neuron/glia proportions. The use of **glial-derived Dx_dev** also reduces direct overlap between the molecular predictor and the neuronal composition outcome.


**Outputs (written under `output/`):**

| File | Description |
| --- | --- |
| `SEAAD-MTG-Dx_dev.csv` | per sample Dx_dev for each major cell class |
| `SEAAD_CPS_subtype_crumblr.csv` | Neuronal-subtype association statistics for CPS |
| `SEAAD_Dx_dev_subtype_crumblr.csv` | Neuronal-subtype association statistics for glia-derived Dx_dev, adjusted for CPS |
| `SEAAD_scCODA_CPS_summary.csv` | Reference-averaged scCODA CPS `Final Parameter`, SEM across reference models, and mean inclusion probability |
| `SEAAD_MTG_neuron_cell_composition_change_Dx_dev_vs_CPS.pdf` | Four panels: two CPS-versus-Dx_dev crumblr plots and two scCODA-versus-crumblr CPS plots |


### 6. Recorded full-data compute resources

Previous Dx_cell projection of the full SEA-AD MTG input was run on an LSF-managed HPC cluster with **8 allocated CPU threads**, completing in **approximately 1 hour**. The scheduler reported **329.6 GB peak memory** and **194 GB average memory**; approximately **400 GB** was requested.

### 7. Reproduction notes

The published SEA-AD scCODA output contains CPS coefficients fitted using rotating reference cell populations. The demo reduces this reference-by-cell-type table to one entry per cell type by **averaging the thresholded `Final Parameter` values across available reference models**, along with the posterior inclusion probabilities. `mean_pip > 0.8` is used as a reference-robustness indicator. The reported SEM quantifies variation **across dependent reference fits**, not independent donor-level sampling uncertainty. Do not interpret this average as an unshrunk or unbiased CLR regression coefficient.

The crumblr comparison uses the **CPS coefficient** from a centered log-ratio analysis (`scale(CPS)`). The scCODA coefficient might use a different CPS scaling and covariate set; effect magnitudes should not be compared directly without harmonizing both. 

### References

1. **SEA-AD dataset**

   Gabitto, M. I., Travaglini, K. J., Rachleff, V. M., et al. (2024). Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience*, 27, 2366–2383. https://doi.org/10.1038/s41593-024-01774-5

2. **PsychAD multiome analytical framework**

   Roussos, P., Dong, P., Cao, X., et al. (2026). Coordinated cellular responses track cognitive outcomes in Alzheimer's disease. *Research Square*, preprint, version 1. https://www.researchsquare.com/article/rs-11296978/v1
