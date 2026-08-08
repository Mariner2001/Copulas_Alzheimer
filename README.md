# A tail-aware vine-copula model for biomarker dependence in Alzheimer's disease

Code accompanying the manuscript *"A tail-aware vine-copula model for biomarker dependence in
Alzheimer's disease"*.

The analysis fits a truncated regular vine copula to ADNI-GO/ADNI-2 baseline data
(n = 630, p = 14), annotates each edge with a co-severity coefficient κ = λ_L, and tests whether
the lower-tail dependence found in the pooled cohort survives removal of clinical-stage pooling.
It does not: this repository reproduces that negative result.

---

## ⚠️ Data are not included

ADNI's Data Use Agreement prohibits redistribution. **No ADNI data are in this repository, and none
should ever be committed to it.** `.gitignore` is configured to prevent this; please do not override
it. Obtain the data yourself from <https://ida.loni.usc.edu> and place the files in `data/` as
described in [`data/README.md`](data/README.md).

---

## Installation

```bash
git clone https://github.com/Mariner2001/Copulas_Alzheimer
cd Copulas_Alzheimer
pip install -r requirements.txt        # or: conda env create -f environment.yml
```

R is needed for two steps only: reading the `ADNIMERGE2` package (step 00) and the Gaussian
copula benchmark. R 4.3.3 with `BDgraph` 2.74; the script installs `BDgraph` if missing.

---

## Reproducing the analysis

Run from the repository root, in order. Total runtime ≈ 80 minutes on one CPU core, almost all of
it step 02.

| # | Command | Produces | Time |
|---|---------|----------|------|
| 00 | `Rscript src/00_extract_adnimerge2.R` | `data/adsl.csv`, `arm.csv`, `uwnp.csv` | ~10 s |
| 01 | `python src/01_merge.py` | `data/adni_merged.csv`, `output/vine_edges.csv`, `output/flow.csv` | ~30 s |
| 02 | `python src/02_bootstrap.py 1000` | `output/boot_raw.pkl` | **~75 min** |
| 03 | `python src/03_aggregate.py` | `output/bootstrap_agg.csv` | ~2 s |
| 04 | `python src/04_tail_evidence.py` | `output/tail_evidence.csv`, `lambda_curves.npz` | ~20 s |
| 05 | `python src/05_stratified.py` | `output/stratified.csv` | ~60 s |
| 06 | `python src/06_stratified_evidence.py` | `output/strat_evidence.csv` | ~30 s |
| 07 | `python src/07_jitter_only.py 110` | `output/jit_raw.pkl`, `jitter_agg.csv` | ~9 min |
| 08 | `python src/08_figures.py` | `figures/fig1–fig4.pdf` | ~40 s |
| 09 | `Rscript gcgm_benchmark.R` | `output/gcgm_results.csv`, `gcgm_summary.txt` | ~4 min |

**Step 02 is resumable.** It appends to `output/boot_raw.pkl` and saves every 10 resamples, so it
can be run in chunks and interrupted safely:

```bash
python src/02_bootstrap.py 200     # run four more times to reach 1000
```

Check progress at any point:

```bash
python -c "import pickle; print(len(pickle.load(open('output/boot_raw.pkl','rb'))))"
```

All random seeds are fixed, so results are reproducible bit-for-bit given the same ADNI release.
Different ADNI releases will shift the numbers slightly.

---

## What each script does

- **`src/vinelib.py`** — shared library. Variable list and orientation (biomarkers oriented so that
  smaller = more pathology; amyloid reflected), jittering of discrete margins, vine fitting,
  extraction of each edge's conditioned pair and conditioning set, and **analytic** tail-dependence
  coefficients per copula family and rotation (not numerical approximations).
- **`00_extract_adnimerge2.R`** — pulls `ADSL`, `ARM` and `UWNPSYCHSUM` out of the `ADNIMERGE2` R
  package. `ARM` is what preserves the EMCI/LMCI distinction that derived tables collapse.
- **`01_merge.py`** — the merge. Filters `ORIGPROT ∈ {ADNIGO, ADNI2}`, joins the five sources,
  averages left/right hemispheres, scales volumes by ICV, applies `OVERALLQC == "Pass"`, and
  reduces to complete cases. Prints the attrition flow.
- **`02_bootstrap.py`** — case bootstrap, B resamples, fresh jitter at each. Records for every
  monitored pair whether it is retained, its tree, family, τ and λ_L, and its conditioning set.
- **`03_aggregate.py`** — inclusion frequencies with Monte Carlo standard errors, percentile
  intervals pooled only over resamples sharing the same conditioning set, and family distributions.
- **`04_tail_evidence.py`** — ΔBIC and Schwarz-corrected Vuong tests against the best radially
  symmetric alternative; observed corner shares against a Gaussian copula at the same τ; λ̂_L(u).
- **`05_stratified.py` / `06_stratified_evidence.py`** — the three ways of breaking the
  stage mixture: within-stratum ranks, diagnosis as an ordinal node, and stratified refits.
- **`07_jitter_only.py`** — sample held fixed, jitter redrawn. Separates identification ambiguity
  (Proposition 4) from sampling variability.
- **`08_figures.py`** — the four manuscript figures.
- **`gcgm_benchmark.R`** — Gaussian copula graphical model on the same data. This is the controlled
  comparison in the paper. Run step 01 first so `output/vine_edges.csv` exists; the script then
  also reports the vine/GCGM edge overlap.

---

## Expected results

Step 01 should print exactly:

```
>>> CASOS COMPLETOS (p=14)   n=630
CN 156 · EMCI 249 · LMCI 128 · AD 97
```

Headline numbers after the full pipeline:

| | Memory–FDG | Memory–Amyloid | PCC–Precuneus |
|---|---|---|---|
| Family (full sample) | Clayton, κ = 0.47 | Clayton, κ = 0.39 | Gaussian, κ = 0 |
| Bootstrap EIP | 1.00 | 1.00 | 0.95 |
| % resamples with λ_L > 0 | 100 % | 57 % | 50 % |
| κ median [5–95 %] | 0.48 [0.42, 0.53] | 0.36 [0.00, 0.44] | 0.00 [0.00, 0.39] |
| ΔBIC vs symmetric | 58.1 | 7.4 | — |
| Vuong | z = 3.72, p < 0.001 | z = 0.55, p = 0.58 | — |
| Lower/upper corner ratio | 2.00 | 1.88 | 1.24 |
| **Within-stratum ranks** | **Student-t, ratio 1.00** | Frank, κ = 0 | Gaussian, κ = 0 |

GCGM benchmark (step 09): 33 edges at posterior EIP > 0.5, of which 20 are shared with the vine,
5 are vine-only and 13 GCGM-only. The age–memory edge has EIP = 0.01 and partial correlation
2×10⁻⁵ — essentially zero, matching the vine, which does not retain it either.

Two calibration figures worth noting: the Gaussian-fitted ExecFun–Memory reference pair still
selects a lower-tail family in **36 %** of resamples, and PCC–Precuneus's corner ratio (1.24) is
*lower* than that same Gaussian-fitted pair's (1.27).

---

Aggregate outputs are included in [`results_published/`](results_published/) so the figures in the
paper can be checked without applying for ADNI access.

## Citation

```bibtex
@article{TBD,
  title   = {A tail-aware vine-copula model for biomarker dependence in Alzheimer's
             disease},
  author  = {TBD},
  journal = {TBD},
  year    = {2026}
}
```

Data collection and sharing for this project was funded by the Alzheimer's Disease Neuroimaging
Initiative (ADNI) (National Institutes of Health Grant U01 AG024904). ADNI investigators
contributed to the design and implementation of ADNI and provided data but did not participate in
the analysis or writing of this report. A complete listing is available at
<http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf>.

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). Data: not included; governed by ADNI's Data Use Agreement.
