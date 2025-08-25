// Efficiency_AutoCent.C
#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TSystem.h>
#include <vector>
#include <utility>
#include <iostream>
#include <TString.h>
#include <TLatex.h>


using std::cout;
using std::endl;

// Axes
constexpr int kPtAxis   = 0;  // pT
constexpr int kCentAxis = 1;  // centrality

// ---------- Styling / Options ----------
int colors[]  = {kRed+1, kBlue+1, kGreen+2, kOrange+7, kMagenta+1, kCyan+2, kAzure+1, kPink+1, kTeal+3, kViolet+1};
int markers[] = {20, 21, 22, 23, 33, 34, 25, 26, 27, 28};

// Inset options (enable + position + y range)
const bool   kUseInset    = true;
// Inset in top-left
const double kInsetX1 = 0.15;  // left edge
const double kInsetY1 = 0.70;  // bottom edge
const double kInsetX2 = 0.45;  // right edge
const double kInsetY2 = 0.88;  // top edge

const bool   kInsetZoomY  = true;   // if true, apply custom y-range
// ---- TPC inset zoom ----
const bool   kInsetZoomTPC  = true;
const double kInsetYminTPC  = 0.35;
const double kInsetYmaxTPC  = 0.7;
const double kInsetXminTPC  = 0.;
const double kInsetXmaxTPC  = 5.;

// ---- TOF inset zoom ----
const bool   kInsetZoomTOF  = true;
const double kInsetYminTOF  = 0.10;
const double kInsetYmaxTOF  = 0.30;
const double kInsetXminTOF  = 0.;
const double kInsetXmaxTOF  = 5.;


// Main pad y-range (set to <0 to auto)
const double kMainYmin = -1;      // e.g. 0.0 to force lower bound; <0 means auto
const double kMainYmax = -1;      // e.g. 1.05 to force upper bound; <0 means auto

// ---------- Helpers ----------
static void ResetRanges(THnSparseD* h) {
  if (!h) return;
  for (int ia=0; ia<h->GetNdimensions(); ++ia) h->GetAxis(ia)->SetRange(0,0);
}

static TH1D* ProjectPt(THnSparseD* h, double cmin, double cmax, const char* name) {
  if (!h) return nullptr;
  TAxis* axCent = h->GetAxis(kCentAxis);
  if (!axCent) return nullptr;
  axCent->SetRangeUser(cmin, cmax);
  TH1D* hp = (TH1D*)h->Projection(kPtAxis);
  if (!hp) return nullptr;
  hp->SetDirectory(nullptr);
  hp->SetName(name);
  hp->Sumw2();
  return hp;
}

static TEfficiency* MakeEff(THnSparseD* num, THnSparseD* den,
                            double cmin, double cmax,
                            Color_t col, int mstyle) {
  ResetRanges(den);
  ResetRanges(num);
  TH1D* hD = ProjectPt(den, cmin, cmax, Form("den_%.0f_%.0f", cmin, cmax));
  TH1D* hN = ProjectPt(num, cmin, cmax, Form("num_%.0f_%.0f", cmin, cmax));
  if (!hD || !hN) { delete hD; delete hN; return nullptr; }
  if (hD->Integral() <= 0) { delete hD; delete hN; return nullptr; }

  auto eff = new TEfficiency(*hN, *hD);
  eff->SetStatisticOption(TEfficiency::kFCP);
  eff->SetLineColor(col);
  eff->SetMarkerColor(col);
  eff->SetMarkerStyle(mstyle);
  eff->SetMarkerSize(1.);
  eff->SetLineWidth(2);
  // optional: vary line style to separate curves more
  // eff->SetLineStyle(1 + (mstyle % 5));

  delete hD; delete hN;
  return eff;
}

// ---- Your custom centrality bins ----
static std::vector<std::pair<double,double>> BuildCentGroupsCustom() {
  return {
    {0,5}, {5,10}, {10,20}, {20,30}, {30,40},
    {40,50}, {50,60}, {60,70}, {70,80}
  };
}

// Draw a full panel (main) and optional inset with the same curves but a custom y-zoom.
static void DrawEfficiencyPanel(THnSparseD* num, THnSparseD* den,
                                const std::vector<std::pair<double,double>>& centGroups,
                                const char* ctitle, const char* ytitle,
                                const char* pngOut,
                                bool doZoom, double yMinZoom, double yMaxZoom, double xMinZoom, double xMaxZoom,
                                const TString& pLabel,
                                const char* pdfOut=nullptr)

{
  TCanvas* c = new TCanvas(Form("c_%s", pngOut), ctitle, 1000, 750);
 // c->SetGrid();
 c->SetLeftMargin(0.1); 
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);

  auto* leg = new TLegend(0.55,0.60,0.88,0.88);
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);

  bool first = true;
  int didx = 0;
  for (const auto& g : centGroups) {
    double cmin = g.first;
    double cmax = g.second;

    TEfficiency* eff = MakeEff(num, den, cmin, cmax,
                               colors[didx % (int)(sizeof(colors)/sizeof(int))],
                               markers[didx % (int)(sizeof(markers)/sizeof(int))]);
    if (!eff) { ++didx; continue; }
// Main pad y-range (set to <0 to auto)
const double kMainYmin = 0.0;   // put your desired lower bound here (e.g. 0.40)
const double kMainYmax = 1.;   // and upper bound (e.g. 0.62)

    if (first) {
      eff->Draw("AP");
      gPad->Update();
if (auto gr = eff->GetPaintedGraph()) {
  if (kMainYmin >= 0 || kMainYmax >= 0) {
    double ymin = (kMainYmin >= 0 ? kMainYmin : gr->GetYaxis()->GetXmin());
    double ymax = (kMainYmax >= 0 ? kMainYmax : gr->GetYaxis()->GetXmax());
    gr->GetYaxis()->SetRangeUser(ymin, ymax);
  }
  gr->GetYaxis()->SetTitleOffset(1.2);   // <<< bring closer (default is ~1.4–1.6)
  gr->GetYaxis()->SetTitleSize(0.05);    // adjust size
  gr->GetYaxis()->SetTitle("Efficiency");
}

      first = false;
    } else {
      eff->Draw("P SAME");
    }
    leg->AddEntry(eff, Form("%.0f-%.0f%%", cmin, cmax), "lep");
    ++didx;
  }
  leg->Draw();

  // ---- Inset zoom (redraw same curves with different y range) ----
 if (kUseInset) {
    TPad* inset = new TPad("inset","inset", kInsetX1, kInsetY1, kInsetX2, kInsetY2);
    inset->SetFillStyle(0);
    inset->SetLineColor(kGray+2);
    inset->SetGrid();
    inset->Draw();
    inset->cd();

    bool firstIn = true;
    int didxIn = 0;
    for (const auto& g : centGroups) {
      double cmin = g.first;
      double cmax = g.second;

      TEfficiency* eff = MakeEff(num, den, cmin, cmax,
                                 colors[didxIn % (int)(sizeof(colors)/sizeof(int))],
                                 markers[didxIn % (int)(sizeof(markers)/sizeof(int))]);
      if (!eff) { ++didxIn; continue; }

      if (firstIn) {
        eff->Draw("AP");
        gPad->Update();
        if (auto gr = eff->GetPaintedGraph()) {
          if (doZoom) gr->GetYaxis()->SetRangeUser(yMinZoom, yMaxZoom);
          if (doZoom) gr->GetXaxis()->SetRangeUser(xMinZoom, xMaxZoom);
          gr->GetXaxis()->SetTitle("");
          gr->GetYaxis()->SetTitle("Efficiency");
          gr->GetXaxis()->SetLabelSize(0.06);
          gr->GetYaxis()->SetLabelSize(0.06);
        }
        firstIn = false;
      } else {
        eff->Draw("P SAME");
      }
      ++didxIn;
    }
    c->cd();
  }
//INfo
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.08);
latex.SetTextFont(42);
latex.DrawLatex(0.85, 0.8, pLabel.Data());

  c->SaveAs(pngOut);
  if (pdfOut) c->SaveAs(pdfOut);
}
// ===================
// Main entry
// ===================
void Efficiency_allCentralities(int option) {
  gStyle->SetOptStat(0);
  
TString particle_name, charge, particle, p;
  
  switch (option) {
  case 1:
    particle_name = "Pion";
    charge        = "pos";
    particle      = "pi";
    p             = "#pi^{+}";
    break;

  case 2:
    particle_name = "Pion";
    charge        = "neg";
    particle      = "pi";
    p             = "#pi^{-}";
    break;
    
  case 3:
    particle_name = "Kaon";
    charge        = "pos";
    particle      = "ka";
    p             = "K^{+}";
    break;
    
  case 4:
    particle_name = "Kaon";
    charge        = "neg";
    particle      = "ka";
    p             = "K^{-}";
    break;

  case 5:
    particle_name = "Proton";
    charge        = "pos";
    particle      = "pr";
    p             = "p";
    break;

  case 6:
    particle_name = "Proton";
    charge        = "neg";
    particle      = "pr";
    p             = "#bar{p}";
    break;

  default:
    std::cout << "Invalid option chosen" << std::endl;
    return;
}

// Inputs
TString kInputFile   = "AnalysisResults_LHC25f3.root";
TString kDenPath     = Form("tof-spectra/MC/%s/%s/prm/pt/den",     particle.Data(), charge.Data());
TString kNumPath     = Form("tof-spectra/MC/withPID/%s/%s/prm/pt/num",     particle.Data(), charge.Data());
TString kNumTOFPath  = Form("tof-spectra/MC/withPID/%s/%s/prm/pt/numtof", particle.Data(), charge.Data());


  TFile* fIn = TFile::Open(kInputFile,"READ");
  if (!fIn || fIn->IsZombie()) { cout << "Cannot open " << kInputFile << endl; return; }

  auto denHist    = dynamic_cast<THnSparseD*>(fIn->Get(kDenPath));
  auto numHist    = dynamic_cast<THnSparseD*>(fIn->Get(kNumPath));
  auto numTOFHist = dynamic_cast<THnSparseD*>(fIn->Get(kNumTOFPath)); // optional

  if (!denHist || !numHist) {
    cout << "Missing THnSparse:\n"
         << "  den: "    << (denHist? "OK":"MISSING")   << " @ " << kDenPath << "\n"
         << "  num: "    << (numHist? "OK":"MISSING")   << " @ " << kNumPath << "\n"
         << "  numTOF: " << (numTOFHist? "OK":"MISSING")<< " @ " << kNumTOFPath << "\n";
    return;
  }

  TAxis* axCent = denHist->GetAxis(kCentAxis);
  if (!axCent) { cout << "No centrality axis at index " << kCentAxis << endl; return; }

  auto centGroups = BuildCentGroupsCustom();
  cout << "Custom centrality groups: " << centGroups.size() << endl;

  gSystem->Exec(Form("mkdir -p %s/%s", particle_name.Data(), charge.Data()));

// -------- TPC efficiency --------
DrawEfficiencyPanel(numHist, denHist, centGroups,
                    "TPC Efficiency (num/den)",
                    "TPC efficiency",
                    Form("%s/%s/TPC_efficiency_AllCents.png", particle_name.Data(), charge.Data()),
                    kInsetZoomTPC, kInsetYminTPC, kInsetYmaxTPC,
                    kInsetXminTPC, kInsetXmaxTPC,
                    p);

// -------- TOF efficiency --------
if (numTOFHist) {
  DrawEfficiencyPanel(numTOFHist, denHist, centGroups,
                      "TOF Efficiency (numTOF/den)",
                      "TOF efficiency",
                      Form("%s/%s/TOF_efficiency_AllCents.png", particle_name.Data(), charge.Data()),
                      kInsetZoomTOF, kInsetYminTOF, kInsetYmaxTOF,
                      kInsetXminTOF, kInsetXmaxTOF,
                    p);
} else {
    cout << "numTOF not found; skipping TOF panel.\n";
  }
  
// -------- SAVE TEfficiencies to ROOT files --------

// TPC
TFile* fOutTPC = new TFile(Form("%s/%s/TPC_efficiencies.root", particle_name.Data(), charge.Data()),"RECREATE");
int idxTPC = 0;
for (const auto& g : centGroups) {
  double cmin = g.first;
  double cmax = g.second;
  TEfficiency* eff = MakeEff(numHist, denHist, cmin, cmax,
                             colors[idxTPC % (int)(sizeof(colors)/sizeof(int))],
                             markers[idxTPC % (int)(sizeof(markers)/sizeof(int))]);
  if (!eff) continue;
  eff->SetName(Form("TPC_eff_%.0f_%.0f", cmin, cmax));
  eff->Write();
  ++idxTPC;
}
fOutTPC->Close();

// TOF
if (numTOFHist) {
  TFile* fOutTOF = new TFile(Form("%s/%s/TOF_efficiencies.root", particle_name.Data(), charge.Data()),"RECREATE");
  int idxTOF = 0;
  for (const auto& g : centGroups) {
    double cmin = g.first;
    double cmax = g.second;
    TEfficiency* eff = MakeEff(numTOFHist, denHist, cmin, cmax,
                               colors[idxTOF % (int)(sizeof(colors)/sizeof(int))],
                               markers[idxTOF % (int)(sizeof(markers)/sizeof(int))]);
    if (!eff) continue;
    eff->SetName(Form("TOF_eff_%.0f_%.0f", cmin, cmax));
    eff->Write();
    ++idxTOF;
  }
  fOutTOF->Close();
}

  

  cout << "Done.\n";
}

