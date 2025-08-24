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

