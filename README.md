# A New Look at Avoidance Phenomena in Maltese Roots

This repository contains the **data and analysis code** accompanying [*A New Look at Avoidance Phenomena in Maltese Roots*](https://doi.org/10.5281/zenodo.23450593) by Carlo Y. Meloni and Kellen Parker van Dam (2026).

## About the study

The paper investigates restrictions on the co-occurrence of similar consonants in Maltese roots, comparing **Semitic-origin** and **loanword-derived** roots. It also considers corresponding patterns in **Arabic** and **Italian**. Bayesian hurdle and offset models are used to examine whether consonant pairs occur and how frequently they occur relative to expectation.

## Repository contents

| Directory | Contents |
| --- | --- |
| `raw_data/` | Original Maltese, Arabic, and Italian lexical data. |
| `processed_data/` | Root lists and consonant-pair datasets used in the analyses. |
| `data_processing/` | R scripts for preparing roots and a Python script (`build_pairs.py`) for constructing pair datasets. |
| `models/` | R scripts fitting the Maltese and Arabic Bayesian hurdle and offset models. |
| `visualizations/` | R scripts for plotting model predictions, contrasts, and diagnostics. |

## Reproducing the analyses

The processed datasets are already included, so the models can be fitted without rerunning data preparation. Install **R** with the packages used by the scripts (notably `tidyverse`, `brms`, `posterior`, `marginaleffects`, and `ggdist`), together with a working **Stan** backend.

From the repository root, run:

```bash
mkdir -p results
Rscript models/hurdlemodel_maltese.R
Rscript models/hurdlemodel_arabic.R
Rscript models/offset_maltese.R
Rscript models/offset_arabic.R
```

The resulting `.rds` model files are saved in `results/`. Once they have been generated, the scripts in `visualizations/` can be used to examine and plot the results.

To regenerate the processed data, use the R scripts in `data_processing/` followed by `build_pairs.py`. This additionally requires **Python 3** with `numpy`, `pandas`, `pyclts`, `soundvectors`, and compatible CLTS data. Note that the Maltese preprocessing script's output filenames differ slightly from the supplied processed-data filenames; check these before rebuilding datasets.

**Note:** Fitted model objects are not included in this repository. The model scripts provided here cover Maltese and Arabic; pooled cross-language model-fitting scripts are not included.

## Citation

If you use this repository, please cite the paper:

Meloni, Carlo Y., and Kellen Parker van Dam. 2026. *A New Look at Avoidance Phenomena in Maltese Roots*. Authors' accepted manuscript, forthcoming in the *Journal of Maltese Studies*. https://doi.org/10.5281/zenodo.23450593

A BibTeX entry for the manuscript is:

```text
@misc{meloni_vandam_2026_maltese_roots,
  author       = {Meloni, Carlo Y. and van Dam, Kellen Parker},
  title        = {A New Look at Avoidance Phenomena in Maltese Roots},
  year         = {2026},
  doi          = {10.5281/zenodo.23450593},
  url          = {https://doi.org/10.5281/zenodo.23450593},
  note         = {Authors' accepted manuscript; forthcoming in the Journal of Maltese Studies}
}
```