// Efficiency_FitML.c  — Beautified & with reliability metrics
// - Top pad: Original (points) vs Fit1 (erf-plateau) vs Fit2 (logistic+tilt) vs ML (BDTG+smooth)
// - Bottom pad: Ratios (model/original), same x-range (0–5). Ratios exist only from kPTcut.
// - Metrics: weighted RMSE and reduced chi2 shown in a stats box; best method highlighted.
// - Outputs: PNGs in smooth_out/ and graphs in TPC_smoothCompare.root
// Run it as root -l -b -q 'Efficiency_FitML.c("TOF_efficiencies.root","TOF_eff_")' or TPC according to what is needed


#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <TString.h>
#include <TLatex.h>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>

// TMVA (optional ML smoother)
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Reader.h>

using std::cout;
using std::endl;

// ---------------- User knobs ----------------
static const double kPTcut   =1.0;    // start of "tail" region to smooth/compare
static const double kXminAll = 0.0;    // common x-range top+bottom
static const double kXmaxAll = 5.0;
static const double kStep    = 0.02;   // step for drawing smooth curves
static const char*  kOutDir  = "smooth_out";
static const char*  kOutRoot = "TPC_smoothCompare.root";
// ---- Global outputs set by Efficiency_FitML() ----

static TString gDet    = "DET";
static TString gOutDir = ".";
static TString gOutRoot= "out.root";


// -------------------------------------------

static std::vector<std::pair<int,int>> kCentBins = {
  {0,5},{5,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,70},{70,80},{80,100}
};

// ---------- style ----------
static void SetupNiceStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextFont(42);
  gStyle->SetLegendFont(42);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.12);
}

// ---------- helpers ----------

// Detect detector label from prefix ("TPC_eff_", "TOF_eff_", ...)
static TString DetectorFromPrefix(const char* prefix) {
  TString p(prefix);
  if (p.BeginsWith("TOF") || p.Contains("TOF")) return "TOF";
  if (p.BeginsWith("TPC") || p.Contains("TPC")) return "TPC";
  return "DET";
}

static TEfficiency* GetEffByCent(TFile* fin, const char* prefix, int cmin, int cmax) {
  if (!fin) return nullptr;
  const TString candidates[] = {
    Form("%s%d_%d", prefix, cmin, cmax),
    Form("effTPC_%d_%d", cmin, cmax),
    Form("TPCeff_%d_%d", cmin, cmax),
    Form("TPC_eff_%d_%d", cmin, cmax),
    Form("TOF_eff_%d_%d", cmin, cmax) // if prefix mismatch
  };
  for (auto& nm : candidates) { TEfficiency* e=nullptr; fin->GetObject(nm, e); if (e) return e; }
  // fallback: scan keys and pick TEfficiency whose name ends with _cmin_cmax
  TIter it(fin->GetListOfKeys()); TKey* k; TString suf=Form("_%d_%d",cmin,cmax);
  while ((k=(TKey*)it())) {
    if (!TString(k->GetClassName()).Contains("TEfficiency")) continue;
    TString nm = k->GetName();
    if (nm.EndsWith(suf)) { TEfficiency* e=nullptr; fin->GetObject(nm, e); if (e) return e; }
  }
  return nullptr;
}

static std::unique_ptr<TGraphAsymmErrors> EffToGraph(TEfficiency* eff) {
  if (!eff) return nullptr;
  return std::unique_ptr<TGraphAsymmErrors>(eff->CreateGraph());
}

// Tail-only symmetric, with weights from y-errors
static std::unique_ptr<TGraphErrors> TailWithWeights(const TGraphAsymmErrors* g,
                                                     double& xmin, double& xmax) {
  auto out = std::make_unique<TGraphErrors>(); xmin=1e9; xmax=-1e9;
  for (int i=0;i<g->GetN(); ++i) {
    double x,y; g->GetPoint(i,x,y);
    if (!(x==x && y==y)) continue;
    if (x < kPTcut) continue;
    double ey = 0.5*(g->GetErrorYlow(i)+g->GetErrorYhigh(i)); if (ey<=0) ey = 1e-3;
    int p = out->GetN(); out->SetPoint(p, x, y); out->SetPointError(p, 0.0, ey);
    xmin = std::min(xmin,x); xmax = std::max(xmax,x);
  }
  return out;
}

static double Median(std::vector<double> v) {
  if (v.empty()) return NAN; std::sort(v.begin(), v.end()); size_t n=v.size();
  return (n%2? v[n/2] : 0.5*(v[n/2-1]+v[n/2]));
}

static TGraph* GraphFromTF1(TF1* f, double xmin, double xmax) {
  auto* gg = new TGraph();
  for (double x=xmin; x<=xmax; x+=kStep) {
    double y=f->Eval(x); if (y<0) y=0; if (y>1) y=1;
    gg->SetPoint(gg->GetN(), x, y);
  }
  return gg;
}

static TGraph* RatioToOriginal(const TGraph* curve, const TGraphAsymmErrors* orig) {
  if (!curve || !orig) return nullptr;
  auto* r = new TGraph();
  for (int i=0;i<orig->GetN();++i){ double x,y; orig->GetPoint(i,x,y);
    if (x < kPTcut || y<=0) continue; double yc=curve->Eval(x); r->SetPoint(r->GetN(), x, yc/y); }
  return r;
}

// -------- Reliability metrics (tail only) --------
struct Metrics { double rmse_w=INFINITY; double chi2=INFINITY; int ndf=0; };
static Metrics EvalMetrics(const TGraph* model, const TGraphAsymmErrors* orig) {
  Metrics m; double sse_w=0, sumw=0, chi2=0; int n=0;
  for (int i=0;i<orig->GetN();++i) {
    double x,y; orig->GetPoint(i,x,y);
    if (x < kPTcut) continue;
    double ey = 0.5*(orig->GetErrorYlow(i)+orig->GetErrorYhigh(i));
    if (ey<=0) continue;
    double yp = model->Eval(x);
    double e  = y-yp;
    double w  = 1.0/(ey*ey);
    sse_w += w*e*e;
    sumw  += w;
    chi2  += (e*e)/(ey*ey);
    n++;
  }
  if (n>0) {
    m.rmse_w = std::sqrt(sse_w / sumw);
    m.chi2   = chi2;
    // No free dof count (nonlinear fit); use n-#params as a heuristic if you wish.
    m.ndf    = std::max(1, n-4);
  }
  return m;
}

// -------- Two robust fit functions --------
static TF1* Fit_ErfPlateau(TGraphErrors* gt, const char* name) {
  if (!gt || gt->GetN()<4) return nullptr;
  std::vector<double> ys; ys.reserve(gt->GetN());
  double xmin=1e9,xmax=-1e9; for (int i=0;i<gt->GetN();++i){ double x,y; gt->GetPoint(i,x,y); ys.push_back(y); xmin=std::min(xmin,x); xmax=std::max(xmax,x); }
  double ymed = Median(ys);
  TF1* f = new TF1(name, "[0] + [1]*0.5*(1+TMath::Erf((x-[2])/[3]))", xmin, xmax);
  f->SetParameters(ymed*0.8, 0.2, std::min(xmin+0.3, xmax-0.2), 0.6);
  f->SetParLimits(0, 0.0, 1.0);
  f->SetParLimits(1,-0.5, 0.5);
  f->SetParLimits(2, xmin-1.0, xmax+1.0);
  f->SetParLimits(3, 0.1, 5.0);
  gt->Fit(f, "Q0W"); // weighted
  return f;
}

static TF1* Fit_LogisticTilt(TGraphErrors* gt, const char* name) {
  if (!gt || gt->GetN()<4) return nullptr;
  std::vector<double> ys; ys.reserve(gt->GetN());
  double xmin=1e9,xmax=-1e9; for (int i=0;i<gt->GetN();++i){ double x,y; gt->GetPoint(i,x,y); ys.push_back(y); xmin=std::min(xmin,x); xmax=std::max(xmax,x); }
  double ymed = Median(ys);
  TF1* f = new TF1(name, "[0] + [1]/(1+exp(-(x-[2])/[3])) + [4]*x", xmin, xmax);
  f->SetParameters(ymed*0.7, 0.3, std::min(xmin+0.4, xmax-0.2), 0.7, 0.01);
  f->SetParLimits(0, 0.0, 1.0);
  f->SetParLimits(1,-1.0, 1.0);
  f->SetParLimits(2, xmin-1.0, xmax+1.0);
  f->SetParLimits(3, 0.1, 5.0);
  f->SetParLimits(4,-0.1, 0.1);
  gt->Fit(f, "Q0W");
  return f;
}

// -------- ML smoother + post-smoothing --------
static void moving_average(std::vector<double>& y, int win=2) {
  if ((int)y.size()<=2*win) return;
  std::vector<double> s(y.size());
  for (int i=0;i<(int)y.size();++i){
    int a=std::max(0,i-win), b=std::min((int)y.size()-1,i+win);
    double sum=0; int n=0; for (int j=a;j<=b;++j){ sum+=y[j]; ++n; }
    s[i]=sum/n;
  }
  y.swap(s);
}

static void isotonic_non_decreasing(std::vector<double>& y) {
  if (y.empty()) return;
  std::vector<double> vals; std::vector<int> lens;
  vals.reserve(y.size()); lens.reserve(y.size());
  for (double v : y) {
    vals.push_back(v); lens.push_back(1);
    while (vals.size()>=2 && vals[vals.size()-2] > vals.back()) {
      int m=lens.back(); lens.pop_back();
      int n=lens.back(); lens.pop_back();
      double a=vals.back(); vals.pop_back();
      double b=vals.back(); vals.pop_back();
      double merged=(a*m+b*n)/(m+n);
      vals.push_back(merged); lens.push_back(m+n);
    }
  }
  std::vector<double> out; out.reserve(y.size());
  for (size_t i=0;i<vals.size();++i) for (int k=0;k<lens[i];++k) out.push_back(vals[i]);
  y.swap(out);
}

static void clamp01(std::vector<double>& y){ for (double& v: y){ if (v<0) v=0; if (v>1) v=1; } }

// -------- ML smoother + post-smoothing (robust against tiny samples) --------
// -------- ML smoother + post-smoothing (robust for small samples) --------
static TGraph* SmoothML(const TGraphAsymmErrors* g, const char* tag) {
  // collect tail points (x >= kPTcut)
  std::vector<double> xs, ys;
  xs.reserve(g->GetN()); ys.reserve(g->GetN());
  for (int i=0; i<g->GetN(); ++i) {
    double x,y; g->GetPoint(i,x,y);
    if (!(x==x && y==y)) continue;
    if (x < kPTcut) continue;
    xs.push_back(x); ys.push_back(y);
  }

  // Guard: skip ML when not enough points (avoids TMVA edge cases)
  const int kMinTailPoints = 20;   // try 20 if you want to allow smaller tails
  if ((int)xs.size() < kMinTailPoints) {
    std::cout << "[" << tag << "] ML skipped: too few tail points (" 
              << xs.size() << " < " << kMinTailPoints << ")\n";
    return nullptr;
  }

  // Prepare TMVA
  TMVA::Tools::Instance();
  TString outdir = Form("%s/tmva_%s", kOutDir, tag);
  gSystem->mkdir(outdir, kTRUE);
  TString cwd = gSystem->WorkingDirectory();

  TGraph* gml = nullptr;  // default if anything goes wrong

  try {
    gSystem->ChangeDirectory(outdir);

    // training tree from tail points
    auto* tTrain = new TTree("Train","Train");
    float pt, eff;
    tTrain->Branch("pt",&pt,"pt/F");
    tTrain->Branch("eff",&eff,"eff/F");
    for (size_t i=0; i<xs.size(); ++i) { pt = xs[i]; eff = ys[i]; tTrain->Fill(); }

    TMVA::DataLoader loader("dl");
    loader.AddVariable("pt","pT","GeV/c",'F');
    loader.AddTarget("eff");
    loader.AddRegressionTree(tTrain, 1.0); // all goes to "Train" (we won't split)

    // Only train; do not Test/Evaluate on tiny sets
    TFile* fout = TFile::Open("tmva.root","RECREATE");
    TMVA::Factory factory("F", fout, "!V:!Silent:Color:AnalysisType=Regression");

    // IMPORTANT: disable bagging on small samples to avoid empty bag crashes
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
                       "!H:!V:"
                       "NTrees=400:"
                       "Shrinkage=0.08:"
                       "MaxDepth=2:"
                       "BoostType=Grad:"
                       "nCuts=20:"
                       "UseBaggedBoost=False");

    factory.TrainAllMethods();
    // factory.TestAllMethods();      // intentionally skipped
    // factory.EvaluateAllMethods();  // intentionally skipped
    fout->Close();

    // Absolute path to weights
    TString wAbs = gSystem->ConcatFileName(outdir, "dl/weights/F_BDTG.weights.xml");

    gSystem->ChangeDirectory(cwd);  // restore CWD

    if (gSystem->AccessPathName(wAbs)) {
      std::cout << "[" << tag << "] ML skipped: weights not found at " << wAbs << "\n";
      return nullptr;
    }

    // Reader & prediction on a dense tail grid
    TMVA::Reader reader("!Color:Silent");
    float pt_in; reader.AddVariable("pt",&pt_in);
    reader.BookMVA("BDTG", wAbs);

    double xmin = *std::min_element(xs.begin(), xs.end());
    double xmax = *std::max_element(xs.begin(), xs.end());
    std::vector<double> gx, gy;
    for (double x=xmin; x<=xmax; x+=kStep) {
      pt_in = x;
      double yhat = reader.EvaluateRegression("BDTG")[0];
      gx.push_back(x); gy.push_back(yhat);
    }

    // Post-smoothing: moving average + isotonic (monotone) + clamp to [0,1]
    moving_average(gy, 2);           // increase to 3 if you want even smoother
    isotonic_non_decreasing(gy);     // removes random rises/dips
    clamp01(gy);

    gml = new TGraph();
    for (size_t i=0; i<gx.size(); ++i) gml->SetPoint(gml->GetN(), gx[i], gy[i]);
    gml->SetLineColor(kGreen+2); gml->SetLineStyle(7); gml->SetLineWidth(4);
    return gml;
  }
  catch (const std::exception& e) {
    gSystem->ChangeDirectory(cwd);
    std::cerr << "[" << tag << "] ML skipped due to TMVA exception: " << e.what() << "\n";
    if (gml) { delete gml; gml=nullptr; }
    return nullptr;
  }
  catch (...) {
    gSystem->ChangeDirectory(cwd);
    std::cerr << "[" << tag << "] ML skipped due to unknown TMVA exception\n";
    if (gml) { delete gml; gml=nullptr; }
    return nullptr;
  }
}

static TGraph* Smooth(const TGraphAsymmErrors* g, const char* tag) {
  // collect tail points (x >= kPTcut)
  std::vector<double> xs, ys;
  for (int i=0; i<g->GetN(); ++i) {
    double x,y; g->GetPoint(i,x,y);
    if (!(x==x && y==y)) continue;
    if (x < kPTcut) continue;
    xs.push_back(x); ys.push_back(y);
  }

  if (xs.size() < 5) { // too few points → skip
    std::cout << "[" << tag << "] ML skipped: too few points (" 
              << xs.size() << ")\n";
    return nullptr;
  }

  // sort by x just in case
  std::vector<std::pair<double,double>> pts;
  for (size_t i=0;i<xs.size();++i) pts.push_back({xs[i],ys[i]});
  std::sort(pts.begin(), pts.end());
  for (size_t i=0;i<pts.size();++i){ xs[i]=pts[i].first; ys[i]=pts[i].second; }

  // smoothing: moving average
  moving_average(ys,2);

  // enforce monotonic (non-decreasing)
  isotonic_non_decreasing(ys);

  // clamp to [0,1]
  clamp01(ys);

  // --- tail fix: blend to plateau after 3 GeV ---
  const double FIX_PULL_START = 3.0;
  const double FIX_BLEND_LEN  = 1.0;

  // plateau = median of data y for x >= 3
  double plateau = 0.9;
  {
    std::vector<double> pv;
    for (int i=0;i<g->GetN();++i){
      double xx,yy; g->GetPoint(i,xx,yy);
      if (xx>=FIX_PULL_START) pv.push_back(yy);
    }
    if (!pv.empty()) {
      std::sort(pv.begin(), pv.end());
      plateau = pv[pv.size()/2];
      plateau = std::min(1.0,std::max(0.0,plateau));
    }
  }

  for (size_t i=0;i<xs.size();++i) {
    if (xs[i] >= FIX_PULL_START) {
      double w = std::min(1.0,(xs[i]-FIX_PULL_START)/FIX_BLEND_LEN);
      ys[i] = (1-w)*ys[i] + w*plateau;
      if (ys[i] > plateau) ys[i] = plateau; // prevent overshoot
    }
  }

  // build graph
  TGraph* gml = new TGraph();
  for (size_t i=0;i<xs.size();++i) gml->SetPoint(gml->GetN(), xs[i], ys[i]);
  gml->SetLineColor(kGreen+2); gml->SetLineStyle(7); gml->SetLineWidth(3);
  return gml;
}
// =============================================================
// Hybrid histogram: data bins below kPTcut, model bins above
// =============================================================
static TH1D* MakeHybridHist(const TEfficiency* eff,
                            const char* hname,
                            TGraph* model) {
  if (!eff) return nullptr;
  auto* htot = eff->GetTotalHistogram();
  auto* hout = (TH1D*)htot->Clone(hname);
  hout->Reset();

  for (int ib=1; ib<=htot->GetNbinsX(); ++ib) {
    double x = htot->GetXaxis()->GetBinCenter(ib);
    double y = 0.0;

    if (x < kPTcut || !model) {
      // use data directly
      y = eff->GetEfficiency(ib);
    } else {
      // evaluate model
      y = model->Eval(x);
    }
    if (y<0) y=0; if (y>1) y=1;
    hout->SetBinContent(ib, y);
  }
  return hout;
}

// ---------- Draw/save (beautified) ----------
static void DrawAndSave(const char* tag,
                        TGraphAsymmErrors* gOrig,
                        TGraph* gFit1, TGraph* gFit2, TGraph* gML, TGraph* gS,
                        double xminTail, double xmaxTail)
{
  gSystem->mkdir(kOutDir, kTRUE);
  TCanvas* c = new TCanvas(Form("c_%s", tag), tag, 1100, 850);
  c->cd();
  TPad* pTop = new TPad("pTop","",0,0.30,1,1); pTop->SetBottomMargin(0.03);/* pTop->SetGrid();*/ pTop->Draw();
  TPad* pBot = new TPad("pBot","",0,0.00,1,0.30); pBot->SetTopMargin(0.05); pBot->SetBottomMargin(0.35); pBot->SetGrid(); pBot->Draw();

  // --- top pad: 0–5, large fonts
  pTop->cd();
  gOrig->SetMarkerStyle(20); gOrig->SetMarkerSize(1.2); gOrig->SetMarkerColor(kBlack);
  gOrig->SetLineColor(kBlack);
  gOrig->GetXaxis()->SetLimits(kXminAll, kXmaxAll);
  gOrig->GetYaxis()->SetRangeUser(0.0, 1.05);
  gOrig->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gOrig->GetYaxis()->SetTitle("Efficiency");
  gOrig->GetXaxis()->SetTitleSize(0.055);
  gOrig->GetYaxis()->SetTitleSize(0.055);
  gOrig->GetXaxis()->SetLabelSize(0.045);
  gOrig->GetYaxis()->SetLabelSize(0.045);
  gOrig->Draw("AP");
  // --- Reliability metrics and "best" model (lowest weighted RMSE)
  Metrics m1, m2, mM, mS;
  if (gFit1) m1 = EvalMetrics(gFit1, gOrig);
  if (gFit2) m2 = EvalMetrics(gFit2, gOrig);
  if (gML)   mM = EvalMetrics(gML,  gOrig);
  if (gS)   mS = EvalMetrics(gS,  gOrig);

  // Find best (by weighted RMSE)
  TString best = "N/A"; double bestRMSE = 1e9;
  if (gFit1 && m1.rmse_w < bestRMSE) { bestRMSE=m1.rmse_w; best="Fit1 (erf)"; }
  if (gFit2 && m2.rmse_w < bestRMSE) { bestRMSE=m2.rmse_w; best="Fit2 (logistic+tilt)"; }
  if (gML  && mM.rmse_w < bestRMSE) { bestRMSE=mM.rmse_w; best="ML (BDTG)"; }
  if (gS  && mS.rmse_w < bestRMSE) { bestRMSE=mS.rmse_w; best="Smooth Constrain"; }

TPaveText* box = new TPaveText(0.47,0.58,0.88,0.88, "NDC"); // <-- top band
box->SetFillColor(0);
box->SetFillStyle(0);
box->SetBorderSize(0);
box->SetTextFont(42);
box->SetTextSize(0.030);  // small enough to avoid hiding data
box->AddText(Form("#bf{Tail metrics (x #geq %.1f GeV/c)}", kPTcut));
if (gFit1) box->AddText(Form("Fit1: RMSE_{w}=%.4f   #chi^{2}/ndf=%.2f", m1.rmse_w, m1.chi2/std::max(1,m1.ndf)));
if (gFit2) box->AddText(Form("Fit2: RMSE_{w}=%.4f   #chi^{2}/ndf=%.2f", m2.rmse_w, m2.chi2/std::max(1,m2.ndf)));
if (gML)   box->AddText(Form("ML  : RMSE_{w}=%.4f   #chi^{2}/ndf=%.2f", mM.rmse_w, mM.chi2/std::max(1,mM.ndf)));
if (gS)   box->AddText(Form("Smooth Constrain  : RMSE_{w}=%.4f   #chi^{2}/ndf=%.2f", mS.rmse_w, mS.chi2/std::max(1,mS.ndf)));
box->AddText(Form("#bf{Best: %s}", best.Data()));
box->Draw();
  if (gFit1){ gFit1->SetLineColor(kRed+1); gFit1->SetLineWidth(4); gFit1->Draw("L SAME"); }
  if (gFit2){ gFit2->SetLineColor(kBlue);   gFit2->SetLineWidth(4); gFit2->SetLineStyle(2); gFit2->Draw("L SAME"); }
  if (gML)  {                                 /* already styled */      gML->Draw("L SAME"); }
  if (gS)  { gS->SetLineColor(kCyan-7);;   gS->SetLineWidth(4); gS->SetLineStyle(4); gS->Draw("L SAME"); }

  // Legend with full functional forms
 auto* leg = new TLegend(0.15, 0.65, 0.45, 0.88);


  leg->SetBorderSize(0); leg->SetFillStyle(0);
leg->SetHeader(Form("%s", gDet.Data()), "C");

  leg->AddEntry(gOrig, "Original efficiency (data points)","lep");
  if (gFit1) leg->AddEntry(gFit1, "Fit1: p0 + p1 * 0.5 * (1+erf((x-p2)/p3))","l");
  if (gFit2) leg->AddEntry(gFit2, "Fit2: a + b/(1+exp(-(x-c)/d)) + e*x","l");
  if (gML)   leg->AddEntry(gML,   "ML smoother (TMVA BDTG + smoothing)","l");
  if (gS)   leg->AddEntry(gS,   "Smooth Constrain","l");
  leg->Draw();

  // vertical line at ptcut
 // TLine lcut(kPTcut, 0.0, kPTcut, 1.05); lcut.SetLineStyle(3); lcut.SetLineColor(kGray+2); lcut.Draw("SAME");

  // --- bottom pad: same 0–5, thick ratio lines, reference=1
  pBot->cd();
  TH1D* frame = new TH1D("frame","",100, kXminAll, kXmaxAll);
  frame->GetYaxis()->SetTitle("model / data");
  frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetRangeUser(0.8, 1.2);
  frame->GetXaxis()->SetTitleSize(0.12);
  frame->GetYaxis()->SetTitleSize(0.11);
  frame->GetXaxis()->SetLabelSize(0.10);
  frame->GetYaxis()->SetLabelSize(0.10);
  frame->Draw("AXIS");

  TLine l1(kXminAll, 1.0, kXmaxAll, 1.0); l1.SetLineStyle(7); l1.SetLineColor(kGray+2); l1.Draw("SAME");
  //TLine lcut2(kPTcut, 0.8, kPTcut, 1.2); lcut2.SetLineStyle(3); lcut2.SetLineColor(kGray+2); lcut2.Draw("SAME");

  auto* r1 = gFit1 ? RatioToOriginal(gFit1, gOrig) : nullptr;
  auto* r2 = gFit2 ? RatioToOriginal(gFit2, gOrig) : nullptr;
  auto* rM = gML   ? RatioToOriginal(gML,  gOrig)  : nullptr;
  auto* rS = gS   ? RatioToOriginal(gS,  gOrig)  : nullptr;

  if (r1){ r1->SetLineColor(kRed+1); r1->SetLineWidth(4); r1->Draw("L SAME"); }
  if (r2){ r2->SetLineColor(kBlue);   r2->SetLineStyle(2); r2->SetLineWidth(4); r2->Draw("L SAME"); }
  if (rM){ rM->SetLineColor(kGreen+2);   rM->SetLineStyle(7); rM->SetLineWidth(4); rM->Draw("L SAME"); }
  if (rS){ rS->SetLineColor(kCyan-7);   rS->SetLineStyle(7); rS->SetLineWidth(4); rS->Draw("L SAME"); }

  /*
  // Make the "best" ratio line slightly thicker
  if (best=="Fit1 (erf)" && r1) r1->SetLineWidth(6);
  if (best=="Fit2 (logistic+tilt)" && r2) r2->SetLineWidth(6);
  if (best=="ML (BDTG)" && rM) rM->SetLineWidth(6);

  // Stats box
  auto* box = new TPaveText(0.12, 0.18, 0.55, 0.89, "NDC"); // left side
  box->SetFillColor(0); box->SetTextFont(42); box->SetTextSize(0.08);
  box->AddText(Form("#bf{Tail metrics (x #geq %.1f GeV/c)}", kPTcut));
  if (gFit1) box->AddText(Form("Fit1: RMSE_{w}=%.4f  #chi^{2}/ndf=%.2f", m1.rmse_w, m1.chi2/std::max(1,m1.ndf)));
  if (gFit2) box->AddText(Form("Fit2: RMSE_{w}=%.4f  #chi^{2}/ndf=%.2f", m2.rmse_w, m2.chi2/std::max(1,m2.ndf)));
  if (gML)   box->AddText(Form("ML  : RMSE_{w}=%.4f  #chi^{2}/ndf=%.2f", mM.rmse_w, mM.chi2/std::max(1,mM.ndf)));
  box->AddText(Form("#bf{Best: %s}", best.Data()));
  box->Draw("same");
*/
c->SaveAs(Form("%s/%s_compare.png", gOutDir.Data(), tag));
TFile* fout = TFile::Open(gOutRoot, "UPDATE");



  fout->cd();
  gOrig->Write(Form("%s_orig",tag), TObject::kOverwrite);
 
  fout->Close();
 delete frame; delete c;
}

// ---------------- Main ----------------
// option: 1=pi+, 2=pi-, 3=K+, 4=K-, 5=p, 6=anti-p
void Efficiency_FitML(int option = 1,
                      const char* finname = "TPC_efficiencies.root",
                      const char* prefix  = "TPC_eff_")

{

TString particle_name, charge, particle, p;
switch (option) {
  case 1: particle_name="Pion";   charge="pos"; particle="pi"; p="#pi^{+}";  break;
  case 2: particle_name="Pion";   charge="neg"; particle="pi"; p="#pi^{-}";  break;
  case 3: particle_name="Kaon";   charge="pos"; particle="ka"; p="K^{+}";    break;
  case 4: particle_name="Kaon";   charge="neg"; particle="ka"; p="K^{-}";    break;
  case 5: particle_name="Proton"; charge="pos"; particle="pr"; p="p";        break;
  case 6: particle_name="Proton"; charge="neg"; particle="pr"; p="#bar{p}";  break;
  default: std::cout<<"Invalid option chosen\n"; return;
}
  TString det     = DetectorFromPrefix(prefix); 
// -------- detector from prefix --------
//TString det = (TString(prefix).BeginsWith("TOF") || TString(prefix).Contains("TOF")) ? "TOF" :
          //    (TString(prefix).BeginsWith("TPC") || TString(prefix).Contains("TPC")) ? "TPC" : "DET";

// -------- base output directories (particle/charge) --------
TString baseDir = Form("%s/%s", particle_name.Data(), charge.Data());
gSystem->Exec(Form("mkdir -p %s", baseDir.Data()));

// -------- global outputs for this run (used by DrawAndSave) --------
gDet     = det;
gOutDir  = Form("%s/smooth_out_%s", baseDir.Data(), det.Data()); 
gOutRoot = Form("%s/%s_smoothCompare.root", baseDir.Data(), det.Data());
gSystem->Exec(Form("mkdir -p %s", gOutDir.Data()));


  SetupNiceStyle();

// Ensure output directory exists
gSystem->Exec(Form("mkdir -p %s", gOutDir.Data()));

  // Detector-specific outputs
//  TString det     = DetectorFromPrefix(prefix);            // "TPC" or "TOF"
  TString outDir  = Form("smooth_out_%s", det.Data());     // e.g. smooth_out_TPC
  TString outRoot = Form("%s_smoothCompare.root", det.Data()); // e.g. TPC_smoothCompare.root

  gSystem->Exec(Form("mkdir -p %s", outDir.Data()));


TString inFilePath = Form("%s/%s/%s", particle_name.Data(), charge.Data(), finname);

TFile* fin = TFile::Open(inFilePath, "READ");
if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open " << inFilePath << "\n";
    return;
}

  for (auto &cb : kCentBins) {
    int cmin = cb.first, cmax = cb.second;

    TEfficiency* eff = GetEffByCent(fin, prefix, cmin, cmax);
    if (!eff) { cout << "Missing " << cmin << "_" << cmax << endl; continue; }

    auto gOrig = EffToGraph(eff);
    if (!gOrig || gOrig->GetN()==0) { cout<<"Empty graph "<<cmin<<"_"<<cmax<<endl; continue; }

    double xminTail, xmaxTail;
    auto gTail = TailWithWeights(gOrig.get(), xminTail, xmaxTail);
    if (gTail->GetN() < 4) {
      cout << "["<<cmin<<"_"<<cmax<<"] Too few tail points; plotting only original.\n";
      DrawAndSave(Form("%d_%d", cmin, cmax), gOrig.get(), nullptr, nullptr, nullptr, nullptr, kPTcut, kPTcut+1.0);
      continue;
    }

    TF1* f1 = Fit_ErfPlateau(gTail.get(), Form("f1_%d_%d", cmin, cmax));
    TF1* f2 = Fit_LogisticTilt(gTail.get(),Form("f2_%d_%d", cmin, cmax));
    TGraph* gFit1 = f1 ? GraphFromTF1(f1, xminTail, xmaxTail) : nullptr;
    TGraph* gFit2 = f2 ? GraphFromTF1(f2, xminTail, xmaxTail) : nullptr;

    TGraph* gML = SmoothML(gOrig.get(), Form("%d_%d", cmin, cmax));
    TGraph* gS = Smooth(gOrig.get(), Form("%d_%d", cmin, cmax));
    
    // Make hybrid efficiency histograms (data <2 GeV, model >=2 GeV)
TH1D* hFit1 = MakeHybridHist(eff, Form("%d_%d_eff_fit1",cmin,cmax), gFit1);
TH1D* hFit2 = MakeHybridHist(eff, Form("%d_%d_eff_fit2",cmin,cmax), gFit2);
TH1D* hML   = MakeHybridHist(eff, Form("%d_%d_eff_ml",  cmin,cmax), gML);
TH1D* hMLc  = MakeHybridHist(eff, Form("%d_%d_eff_mlconstr",cmin,cmax), gS);

TFile* fout = TFile::Open(gOutRoot,"UPDATE");
if (fout && !fout->IsZombie()) {
  if (gOrig) gOrig->Write(Form("%d_%d_orig",cmin,cmax), TObject::kOverwrite);
  if (hFit1) hFit1->Write("",TObject::kOverwrite);
  if (hFit2) hFit2->Write("",TObject::kOverwrite);
  if (hML)   hML->Write("",TObject::kOverwrite);
  if (hMLc)  hMLc->Write("",TObject::kOverwrite);
  fout->Close();
}


    DrawAndSave(Form("%d_%d", cmin, cmax), gOrig.get(), gFit1, gFit2, gML, gS, xminTail, xmaxTail);

    delete f1; delete f2; delete gFit1; delete gFit2; delete gML; delete gS;
  }

  fin->Close();
  cout << "✅ Done. Plots in " << kOutDir << " and graphs in " << kOutRoot << endl;
}

