# UV-C probabilistic bounds (Coffee/Cacao) — Code + Data (Zenodo)

[![DOI (Code)](https://zenodo.org/badge/DOI/10.5281/zenodo.18133165.svg)](https://doi.org/10.5281/zenodo.18133165)
[![DOI (Dataset)](https://zenodo.org/badge/DOI/10.5281/zenodo.18133112.svg)](https://doi.org/10.5281/zenodo.18133112)

Reproducible analysis pipeline for:
1) **probabilistic UV-C inactivation bounds** (exact Clopper–Pearson 95% CI),
2) **physics-informed translation** (σ\_LB, t90/t99/t99.9, energy-per-kill vs. density),
3) a **shape-only morphometric fingerprint** (Gradient Boosting with plate-grouped CV).

Raw Excel inputs **do not need to live inside this repository**. Provide them via paths at runtime.

---

## Archives (Zenodo)

### Code (Software)
- **Concept DOI (all versions):** https://doi.org/10.5281/zenodo.18133165  
- **Version DOI (v1.0.0):** https://doi.org/10.5281/zenodo.18133166  

### Coffee dataset (Dataset)
- **Concept DOI (all versions):** https://doi.org/10.5281/zenodo.18133112  
- **Version DOI (v1.0.0):** https://doi.org/10.5281/zenodo.18133113  

---

## Data sources (Coffee vs. Cacao panels)

### Coffee panel (generated for the Food Control study)
- The Coffee UV-C cohort was generated for the *Food Control* manuscript “From minutes to bounds …” (accepted, in press).
- Dataset archive (Zenodo): https://doi.org/10.5281/zenodo.18133112

### Cacao panel (published reference dataset)
- The Cacao panel originates from:
  Baek et al. (2025) *Scientific Reports* 15:36256. https://doi.org/10.1038/s41598-025-20277-2
- Download the **Supplementary material 2** from the article page and save/rename locally as:
  - `Cacao.xlsx`

---

## Inputs
Place these files in a folder (example: `C:\projects\UV_data\`):
- `Coffee.xlsx`  *(from Zenodo Coffee dataset, above)*
- `Cacao.xlsx`   *(from the Scientific Reports supplementary material, above)*

> Tip: If you prefer not to rename files, pass explicit paths via `--coffee` and `--cacao`.

---

## Setup (Anaconda Prompt)

Create a conda environment:
```bat
conda env create -f environment.yml
conda activate uvc-analysis


## Run
Execute the analysis pipeline by pointing to your data files:

### A) Run with default directory
Place `Coffee.xlsx` and `Cacao.xlsx` in one folder and run:
```bat
python uvc_run_all.py --data-dir "C:/path/to/your/data" --out-dir "outputs"


python uvc_run_all.py ^
  --coffee "C:/data/Coffee.xlsx" ^
  --cacao  "C:/data/Cacao.xlsx" ^
  --out-dir "outputs"



Citation

If you use this code or dataset, please cite:

Ahn et al. (2026). From minutes to bounds: A probabilistic UV-C control and a shape-only morphological fingerprint for postharvest Colletotrichum inactivation in cacao and coffee processing. Food Control. (Accepted, In Press)

Baek et al. (2025). Scientific Reports 15:36256. https://doi.org/10.1038/s41598-025-20277-2
