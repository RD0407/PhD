# Efficiency_FitML : Smooth & Compare Efficiency Curves (ROOT/TMVA)

This repository contains a single ROOT macro, **`Efficiency_FitML.c`**, that reads detector efficiency curves (as `TEfficiency` objects) and compares four tail–smoothing/fit strategies:

- **Fit1: Erf-plateau** (sigmoid-like rise with plateau)
- **Fit2: Logistic+tilt** (logistic with a small linear term)
- **ML: TMVA BDTG** + post-smoothing (moving average + isotonic monotonicity + clamping to [0,1])
- **S: Constrained smoother** (moving average + isotonic + clamp; *no TMVA required*)

The macro produces side-by-side plots of **original points vs. all models**, tail **ratios** (model/original), **reliability metrics** (weighted RMSE and χ² with heuristic NDF), and writes **hybrid efficiency histograms** that use **data below the tail cut** and **model above**.

> ✔️ Out-of-the-box batch usage:
>
> ```bash
> root -l -b -q 'Efficiency_FitML.c(1, "TPC_efficiencies.root", "TPC_eff_")'
> ```
> where the first argument selects the particle/charge (see below).

---

## What the macro does

For each centrality bin (defaults: `0–5, 5–10, …, 80–100`), the macro:

1. Loads a `TEfficiency` from your ROOT file using a **name prefix** + `_cmin_cmax` (e.g. `TPC_eff_0_5`). A few common fallbacks are tried automatically.
2. Converts it to a `TGraphAsymmErrors` and extracts the **tail** region starting at `kPTcut` (default **1.0**).
3. Builds four tail models:
   - **Fit1**: erf-plateau fit (weighted)
   - **Fit2**: logistic with linear tilt (weighted)
   - **ML**: TMVA `BDTG` regression trained on the tail and evaluated on a dense grid; then **moving average → isotonic non-decreasing → clamp to [0,1]**
   - **S**: the same post-smoothing chain applied directly to the tail data (no TMVA)
4. Computes **weighted RMSE** and **χ²** on the tail points only, shows them in a stats box, and highlights the **best** method.
5. Draws a two-pad canvas per centrality:
   - **Top**: original points + all models
   - **Bottom**: **ratios** (model/original) over a common x-range (default **0–5**)
6. Writes **hybrid histograms** for each method: data for `x < kPTcut`, model for `x ≥ kPTcut`:
   - `*_eff_fit1`, `*_eff_fit2`, `*_eff_ml`, `*_eff_mlconstr`

---

## Inputs

- A ROOT file with `TEfficiency` objects named as:
  - `"<prefix><cmin>_<cmax>"` (e.g. `TPC_eff_0_5`)
  - Fallbacks tested automatically include: `effTPC_cmin_cmax`, `TPCeff_cmin_cmax`, `TPC_eff_cmin_cmax`, `TOF_eff_cmin_cmax`.
- The macro infers the **detector tag** (`TPC`/`TOF`) **from the prefix** to label outputs.

### Particle / charge selection (first argument)

`Efficiency_FitML(int option, const char* finname, const char* prefix)`

| option | Particle | Charge | Labels |
|-------:|----------|--------|--------|
| 1 | Pion | + | `#pi^{+}` |
| 2 | Pion | − | `#pi^{-}` |
| 3 | Kaon | + | `K^{+}` |
| 4 | Kaon | − | `K^{-}` |
| 5 | Proton | + | `p` |
| 6 | Proton | − | `#bar{p}` |

---

## Outputs

For each `option` run the macro creates a **particle/charge-specific folder** and puts detector-tagged results inside it.


- **PNG plots**: one per centrality bin (`<cmin>_<cmax>_compare.png`), e.g. `0_5_compare.png`.
- **ROOT file**: `<DET>_smoothCompare.root` contains:
  - Original graphs: `"<cmin>_<cmax>_orig"`
  - Hybrid histograms using data (< `kPTcut`) + model (≥ `kPTcut`):
    - `"<cmin>_<cmax>_eff_fit1"`
    - `"<cmin>_<cmax>_eff_fit2"`
    - `"<cmin>_<cmax>_eff_ml"`
    - `"<cmin>_<cmax>_eff_mlconstr"`

At the end you’ll also see a console summary like:
✅ Done. Plots in smooth_out_<DET> and graphs in <DET>_smoothCompare.root


---

## Requirements

- **ROOT 6** (tested with TMVA enabled). Headers used include:
  - `TFile`, `TEfficiency`, `TGraph(Asymm)Errors`, `TF1`, `TH1D`, `TCanvas`, `TLatex`, `TLegend`, `TStyle`, `TPad`, etc.
  - **TMVA** (optional but recommended) for BDTG:
    - `TMVA/Tools.h`, `TMVA/Factory.h`, `TMVA/DataLoader.h`, `TMVA/Reader.h`
  - If TMVA is **not** available or an exception occurs, **ML** is skipped and the macro still produces **S** (constrained smoother) and the two parametric fits.
- A UNIX-like shell for creating output directories (the macro uses `gSystem->Exec("mkdir -p …")`).

---

## Usage

Run in batch mode from your project directory:

```bash
# General form
root -l -b -q 'Efficiency_FitML.c(<option>, "<in.root>", "<prefix>")'

# Examples
root -l -b -q 'Efficiency_FitML.c(1, "TPC_efficiencies.root", "TPC_eff_")'   # π+
root -l -b -q 'Efficiency_FitML.c(4, "TOF_efficiencies.root", "TOF_eff_")'   # K−
```
The macro loops over centrality bins:
0–5, 5–10, 10–20, 20–30, 30–40, 40–50, 50–60, 60–70, 70–80, 80–100.


Tunable knobs (edit inside the macro)
static const double kPTcut   = 1.0;  // Tail starts here (modeling/metrics)
static const double kXminAll = 0.0;  // Common x-range for the canvas
static const double kXmaxAll = 5.0;
static const double kStep    = 0.02; // Step for drawing smooth curves

    kPTcut controls where the tail begins. Ratios and metrics are computed only for x ≥ kPTcut. Hybrid histograms use data below, model above.

    The BDTG is configured with sensible defaults (NTrees=400, Shrinkage=0.08, MaxDepth=2, etc.) and can be tweaked inside the SmoothML helper.

Methods in a nutshell

    Fit1 (Erf-plateau)
    f(x) = p0 + p1 * 0.5*(1 + erf((x - p2)/p3)) with bounded parameters and weighted fit.

    Fit2 (Logistic+tilt)
    f(x) = p0 + p1 / (1 + exp(-(x - p2)/p3)) + p4 * x (weighted; small tilt allowed).

    ML (BDTG + smoothing)
    Train BDTG on tail (x=pt, y=eff); predict on a dense grid; then:

        Moving average (window 2–3)

        Isotonic regression (enforce non-decreasing)

        Clamp to [0,1]

    S (Constrained smoother)
    Applies the same three post-steps directly to the tail points.

Metrics: on tail points only:

    Weighted RMSE (rmse_w) using inverse-variance weights from asymmetric y-errors condensed to a single σ.

    χ² (with heuristic NDF = max(1, N−4)) reported as “reduced χ²” in the stats box.

Troubleshooting

    “ML skipped: weights not found …”
    TMVA training writes weights under dl/weights/F_BDTG.weights.xml in a temporary output dir. If the file is missing or TMVA isn’t enabled, the macro will skip ML and continue.

    “Too few tail points”
    If there are < 4–5 tail points after kPTcut, parametric fits or ML may be skipped; you will still get the original graph and (usually) the constrained smoother.

    Names don’t match
    Verify your object names match "<prefix><cmin>_<cmax>" or adjust GetEffByCent in the macro.

Development notes

    Output paths are automatically namespaced by particle/charge and detector inferred from prefix.

    Styling (SetupNiceStyle, pad layout, legends, latex labels) is handled inside the macro for publication-quality PNGs.

    The macro writes/updates the ROOT output in UPDATE mode to avoid losing previous results
