# GDFlex

## Overview
This repository accompanies the KDD 2025 paper

> **Generalized Discords for Time‑Series Anomaly Detection with Flexible Subsequence Lengths**  
> DOI: <https://doi.org/10.1145/3711896.3736977>

It contains all MATLAB source code, example scripts, and documentation needed to reproduce the results reported in the paper.

---

## Installation

No external libraries or compilation are required.  
All MATLAB scripts are self‑contained and have been verified to run on **MATLAB R2023b** (R2021a+ should also work).

```bash
git clone https://github.com/<your‑fork>/gdflex.git
cd gdflex
```

---

## Quick Start

```matlab
cd src

% Example 1 – run on a single dataset (ID = 66) with the final GDFlex mode
Pub_GDFlex_execute(66,"interNoise");

% Example 2 – run on every dataset (~250) with the final GDFlex mode
Pub_GDFlex_execute("all","interNoise");
```

`Pub_GDFlex_execute(dataID, mode)`

* **dataID** — an ML/UCR dataset ID (integer) or `"all"`  
* **mode**   — `baseline | znormBias | intra | locMis | inter | interNoise`  
  (`interNoise` = full GDFlex)

### Typical Runtime  
*(Intel i5‑1335U @ 1.30 GHz, 16 GB RAM)*

| Target & Mode | Typical Time | Notes |
|---------------|--------------------------------------------------------|-------|
| **ML/UCR (250)** `interNoise` | ≈ 7 h | Full GDFlex pipeline |
| **ML/UCR (250)** `baseline`   | ≥ 24 h | Slowest (no target-set restriction) |
| **ML/UCR (250)** `znormBias`, `intra`, `locMis`, `inter` | *Faster than* `interNoise` | Ablation modes |

Generated artifacts (per dataset)

* `result_all.csv`     — summary spreadsheet  
* `InfoDemo_<dataID>.mat` — intermediate results  
* Five explanatory figures

For a step‑by‑step guide see **`doc/0_HowToRunGDFlex.pdf`**.

---

## Documentation (`doc/`)

| File                               | Contents / Figures Reproduced |
|------------------------------------|--------------------------------|
| **0_HowToRunGDFlex.pdf**           | Running instructions, output description |
| **1_KeyFigures.pdf**               | Fig. 1, Fig. 7 (A–C) |
| **2_Observations.pdf**             | Fig. 8–10, Fig. 3, Fig. 7 (E), Fig. 15 |
| **3_SpikeletDecomposition.pdf**    | Fig. 6, Fig. 7 |
| **4_Algorithm.pdf**                | Fig. 12–14 |
| **5_Evaluations.pdf**              | Table 2, Fig. 16, Supplementary Figure A, Supplementary Table A |

---

## Citation

If you use GDFlex or the processed datasets in your work, please cite:

Imamura, M., and Nakamura, T. 2025. Generalized Discords for Time Series Anomaly Detection
with Flexible Subsequence Lengths. ACM SIGKDD. <https://doi.org/10.1145/3711896.3736977>

For the original datasets:

Dau, A. H., Bagnall, A., Kamgar, K., Yeh, C.-C. M., Zhu, Y., Gharghabi, S., Ratanamahatana, C. A., and Keogh, E. 2018.  The UCR Time Series Archive. Retrieved from <https://www.cs.ucr.edu/~eamonn/time_series_data_2018>

GDFlex makes use of the distance profile computation programs published in the following paper or provided on its supporting website:

Zhu, Y., Imamura, M., Nikovski, D., and Keogh, E. 2017. Matrix Profile VII: Time Series Chains—A New Primitive for Time Series Data Mining.* In *Proceedings of IEEE ICDM, pp. 695–704.  
<https://doi.org/10.1109/ICDM.2017.79>
---

## Acknowledgments
* The **ML/UCR Time‑Series Archive** for providing benchmark datasets.  
* Additional supporters are listed in the paper’s acknowledgments.
