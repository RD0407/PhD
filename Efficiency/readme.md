# Efficiency_allCentralities

This program is a ROOT-based C++ analysis macro for computing and plotting reconstruction efficiencies as a function of transverse momentum (pT) across different centrality classes in heavy-ion physics analyses.

It is designed to take pre-filled multi-dimensional histograms (THnSparse) that contain numerator and denominator counts, project them onto the pT axis, and calculate efficiencies using TEfficiency. The results are drawn as styled efficiency curves, with optional zoomed inset panels for different detectors (e.g., TPC, TOF).


### 🔍 What the Code Does (Step by Step)

#### Input Handling

The macro expects an input ROOT file containing THnSparse histograms:

Numerator → events that passed selection.

Denominator → all candidate events.

Centrality and pT are stored as different axes of the THnSparse.

The macro defines axis indices (kPtAxis, kCentAxis) so it knows how to project the data.

Projection onto pT

For a given centrality range (cmin, cmax), the function ProjectPt() sets a range on the centrality axis of the THnSparse, and then projects onto the pT axis.

This produces a 1D histogram of yields vs pT.

#### Efficiency Calculation

The function MakeEff() takes the numerator and denominator THnSparse objects.

For each centrality interval, it:

Projects numerator and denominator onto pT.

Constructs a TEfficiency object to calculate efficiency bin-by-bin with proper statistical treatment (Clopper-Pearson confidence intervals).

The efficiency is then styled with a unique color and marker.

#### Plotting

Efficiencies for different centrality bins are drawn on the same canvas.

Each bin gets a different color/marker style (arrays colors[] and markers[] are predefined).

A legend is added to distinguish centrality classes.

Axis ranges (kMainYmin, kMainYmax) can be fixed or set to automatic.

Inset Panels (Zoom Regions)

The code defines inset pads (TPad) that zoom into specific pT ranges.

#### Example:

TPC inset → zoom on efficiency between 0–5 GeV/c with efficiency range 0.35–0.7.

TOF inset → zoom on efficiency 0–5 GeV/c with efficiency range 0.10–0.30.

These are useful to highlight low-efficiency detector regions.

#### Styling & Output

The macro sets consistent styles: marker types, colors, axis labels, legend placement.

The final figure is drawn to a canvas (TCanvas) and can be saved to file.



## ⚙️ Code Structure (Functions)

ResetRanges(THnSparseD* h)
Clears all axis ranges on a THnSparse (so projections are not biased by old selections).

ProjectPt(THnSparseD* h, double cmin, double cmax, const char* name)
Selects a centrality interval and projects the THnSparse onto the pT axis, returning a 1D histogram.

MakeEff(THnSparseD* num, THnSparseD* den, double cmin, double cmax, Color_t col, int mstyle)
Creates a TEfficiency object for a given centrality range, styled with a chosen color and marker.

main() (or macro execution block)

Opens input ROOT file.

Defines centrality bins.

Loops through bins → calls MakeEff() for each.

Draws all efficiency curves on a shared canvas with insets.
