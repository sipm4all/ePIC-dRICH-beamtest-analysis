//  --- --- ---
//  GENERAL UTILITY
//  author: nicola.rubini@bo.infn.it
//  --- --- ---
#pragma once
//  Include
#include "Database.h"
#include "../macros/recowriter.C"
//  --- --- ---
//
inline float get_beta(float mass, float momentum)
{
  return 1. / sqrt(1 + (mass * mass) / (momentum * momentum));
}

inline float get_ctheta(float beta, float rif_index)
{
  return TMath::ACos(1. / (beta * rif_index));
}

inline float get_radius(float mass, float momentum, float rif_index, float arm_length)
{
  return TMath::Tan(get_ctheta(get_beta(mass, momentum), rif_index)) * arm_length;
}

inline float measure_ctheta(float radius, float arm_length)
{
  return TMath::ATan(radius / arm_length);
}

inline float measure_beta(float rif_index, float radius, float arm_length)
{
  return 1. / (rif_index * TMath::Cos(measure_ctheta(radius, arm_length)));
}

inline float get_mass_hypothesis(float momentum, float rif_index, float radius, float arm_length)
{
  return TMath::Sqrt((momentum * momentum) / (measure_beta(rif_index, radius, arm_length) *
                                              measure_beta(rif_index, radius, arm_length)) -
                     (momentum * momentum));
}

std::string data_dir = "./Data/";
std::string sim_dir = "./simulation/data/";

//  X, Y, R and errors
using circle_fit_results = std::array<std::array<float, 2>, 3>;

struct recodata
{
  unsigned short n = 0;
  float x[1024];
  float y[1024];
  float t[1024];
};

struct TObjects
{
  std::map<std::string, TH1F *> _TH1F;
  std::map<std::string, TH1D *> _TH1D;
  std::map<std::string, TH2F *> _TH2F;
  std::map<std::string, TH2I *> _TH2I;
  std::map<std::string, TH3F *> _TH3F;
  std::map<std::string, TProfile *> _TProfile;
  std::map<std::string, TGraph *> _TGraph;
  std::map<std::string, TGraphErrors *> _TGraphErrors;
};

void write_all(TObjects target)
{
  for (auto [label, object] : target._TProfile)
    object->Write(label.c_str());
  for (auto [label, object] : target._TH1F)
    object->Write(label.c_str());
  for (auto [label, object] : target._TH1D)
    object->Write(label.c_str());
  for (auto [label, object] : target._TH2F)
    object->Write(label.c_str());
  for (auto [label, object] : target._TH2I)
    object->Write(label.c_str());
  for (auto [label, object] : target._TH3F)
    object->Write(label.c_str());
  for (auto [label, object] : target._TGraph)
    object->Write(label.c_str());
  for (auto [label, object] : target._TGraphErrors)
    object->Write(label.c_str());
}

void load_production_step(TObjects &target, TFile *input_file, std::string run_tag = "")
{
  target._TH3F["hDeltaRadius"] = (TH3F *)(input_file->Get("hDeltaRadius"));
  target._TH3F["hFitRadius"] = (TH3F *)(input_file->Get("hFitRadius"));
  target._TH2F["hFitXcoordinate"] = (TH2F *)(input_file->Get("hFitXcoordinate"));
  target._TH2F["hFitYcoordinate"] = (TH2F *)(input_file->Get("hFitYcoordinate"));
  target._TH2F["hFitRadiusVal"] = (TH2F *)(input_file->Get("hFitRadiusVal"));
  target._TH2F["hFitMassPrediction"] = (TH2F *)(input_file->Get("hFitMassPrediction"));
  target._TH2F["hPersistance2D"] = (TH2F *)(input_file->Get("hPersistance2D"));
  target._TH3F["hPersistance3D"] = (TH3F *)(input_file->Get("hPersistance3D"));
  target._TH2F["hAveragePhotonEv"] = (TH2F *)(input_file->Get("hAveragePhotonEv"));
  target._TH2F["hAveragePhotonEvDisc"] = (TH2F *)(input_file->Get("hAveragePhotonEvDisc"));
  target._TH2F["hEfficiencyMap"] = (TH2F *)(input_file->Get("hEfficiencyMap"));
  if (!run_tag.empty())
  {
    target._TH3F["hDeltaRadius"]->SetName(Form("hDeltaRadius_%s", run_tag.c_str()));
    target._TH3F["hFitRadius"]->SetName(Form("hFitRadius_%s", run_tag.c_str()));
    target._TH2F["hFitXcoordinate"]->SetName(Form("hFitXcoordinate_%s", run_tag.c_str()));
    target._TH2F["hFitYcoordinate"]->SetName(Form("hFitYcoordinate_%s", run_tag.c_str()));
    target._TH2F["hFitRadiusVal"]->SetName(Form("hFitRadiusVal_%s", run_tag.c_str()));
    target._TH2F["hFitMassPrediction"]->SetName(Form("hFitMassPrediction_%s", run_tag.c_str()));
    target._TH2F["hPersistance2D"]->SetName(Form("hPersistance2D_%s", run_tag.c_str()));
    target._TH3F["hPersistance3D"]->SetName(Form("hPersistance3D_%s", run_tag.c_str()));
    target._TH2F["hAveragePhotonEv"]->SetName(Form("hAveragePhotonEv_%s", run_tag.c_str()));
    target._TH2F["hAveragePhotonEvDisc"]->SetName(Form("hAveragePhotonEvDisc_%s", run_tag.c_str()));
    target._TH2F["hEfficiencyMap"]->SetName(Form("hEfficiencyMap_%s", run_tag.c_str()));
  }
}

void check_recodata(std::string run_tag)
{
  TFile *input_recodatafile = new TFile(Form("%s/%s/recodata.root", data_dir.c_str(), run_tag.c_str()));
  if (input_recodatafile->IsZombie())
    recowriter((data_dir + "/" + run_tag + "/lightdata.root").c_str(), (data_dir + "/" + run_tag + "/recodata.root").c_str());
  input_recodatafile->Close();
}

void load_data(TTree *reco_data_tree, recodata &target_data_struct)
{
  reco_data_tree->SetBranchAddress("n", &target_data_struct.n);
  reco_data_tree->SetBranchAddress("x", &target_data_struct.x);
  reco_data_tree->SetBranchAddress("y", &target_data_struct.y);
  reco_data_tree->SetBranchAddress("t", &target_data_struct.t);
}

TTree *load_data(std::string filename, recodata &target_data_struct)
{
  auto input_file = TFile::Open(filename.c_str());
  auto input_tree = (TTree *)input_file->Get("recodata");
  load_data(input_tree, target_data_struct);
  return input_tree;
}

struct ringdata
{
  unsigned short n = 0;
  float x0[1024];
  float y0[1024];
  float r0[1024];
};

void load_data(TTree *ring_data_tree, ringdata &target_data_struct)
{
  ring_data_tree->SetBranchAddress("N", &target_data_struct.n);
  ring_data_tree->SetBranchAddress("X0", &target_data_struct.x0);
  ring_data_tree->SetBranchAddress("Y0", &target_data_struct.y0);
  ring_data_tree->SetBranchAddress("R", &target_data_struct.r0);
}

TTree *load_data(std::string filename, ringdata &target_data_struct)
{
  auto input_file = TFile::Open(filename.c_str());
  auto input_tree = (TTree *)input_file->Get("ringdata");
  load_data(input_tree, target_data_struct);
  return input_tree;
}

template <bool return_selected = true>
recodata select_points(recodata reco_data, std::array<float, 3> ring_data, float radius_cut)
{
  recodata selected_reco_data;
  recodata discarded_reco_data;
  selected_reco_data.n = 0;
  discarded_reco_data.n = 0;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    auto x_shift = reco_data.x[iPnt] - ring_data[0];
    auto y_shift = reco_data.y[iPnt] - ring_data[1];
    auto current_radius = sqrt(x_shift * x_shift + y_shift * y_shift);
    if (fabs(reco_data.t[0]) > 20)
    {
      discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
      discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
      discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
      discarded_reco_data.n++;
      continue;
    }
    if (fabs(current_radius - ring_data[2]) > radius_cut)
    {
      discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
      discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
      discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
      discarded_reco_data.n++;
      continue;
    }
    selected_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
    selected_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
    selected_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
    selected_reco_data.n++;
  }
  return return_selected ? selected_reco_data : discarded_reco_data;
}

template <bool return_selected = true>
recodata select_points(recodata reco_data, ringdata ring_data, float radius_cut)
{
  recodata selected_reco_data;
  recodata discarded_reco_data;
  selected_reco_data.n = 0;
  discarded_reco_data.n = 0;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    auto x_shift = reco_data.x[iPnt] - ring_data.x0[0];
    auto y_shift = reco_data.y[iPnt] - ring_data.y0[0];
    auto current_radius = sqrt(x_shift * x_shift + y_shift * y_shift);
    if (fabs(reco_data.t[0]) > 20)
    {
      discarded_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
      discarded_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
      discarded_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
      discarded_reco_data.n++;
      continue;
    }
    if (fabs(current_radius - ring_data.r0[0]) > radius_cut)
    {
      discarded_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
      discarded_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
      discarded_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
      discarded_reco_data.n++;
      continue;
    }
    selected_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
    selected_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
    selected_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
    selected_reco_data.n++;
  }
  return return_selected ? selected_reco_data : discarded_reco_data;
}

TGraph *plot(recodata reco_data)
{
  auto result = new TGraph();
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    result->SetPoint(iPnt, reco_data.x[iPnt], reco_data.y[iPnt]);
  return result;
}

TEllipse *plot(ringdata ring_data)
{
  auto result = new TEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0]);
  result->SetFillStyle(0);
  result->SetLineColor(kRed);
  result->DrawEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0], 0, 0, 360, 0, "same");
  return result;
}

TEllipse *plot(ringdata ring_data, float tolerance_cut)
{
  auto result = new TEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0]);
  result->SetFillStyle(0);
  result->SetLineColor(kRed);
  result->SetLineWidth(5);
  result->DrawEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0], 0, 0, 360, 0, "same");
  result->SetLineColor(kRed);
  result->SetLineStyle(kDashed);
  result->SetLineWidth(5);
  result->DrawEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0] + tolerance_cut, 0, 0, 360, 0, "same");
  result->DrawEllipse(ring_data.x0[0], ring_data.y0[0], ring_data.r0[0] - tolerance_cut, 0, 0, 360, 0, "same");
  return result;
}

TEllipse *plot(circle_fit_results ring_fit)
{
  auto result = new TEllipse(ring_fit[0][0], ring_fit[1][0], ring_fit[2][0]);
  result->SetFillStyle(0);
  result->SetLineColor(kRed);
  result->SetLineStyle(kDashed);
  result->DrawEllipse(ring_fit[0][0], ring_fit[1][0], ring_fit[2][0], 0, 0, 360, 0, "same");
  return result;
}

TCanvas *plot_standard_canvas(std::string nametitle = "std_canvas")
{
  TCanvas *standard_canvas = new TCanvas(nametitle.c_str(), nametitle.c_str(), 1000, 1000);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->DrawFrame(-100., -100., 100., 100., ";x (mm); y (mm)");
  return standard_canvas;
};

template <bool fix_XY = true>
circle_fit_results
fit_circle(TGraph *gTarget, std::array<float, 3> initial_values)
{
  circle_fit_results result;

  //  Chi2 minimisation for points in a circle
  auto chi2_function = [&](const double *parameters)
  {
    float chi2 = 0;
    for (int iPnt = 0; iPnt < gTarget->GetN(); iPnt++)
    {
      double delta_x = gTarget->GetPointX(iPnt) - parameters[0];
      double delta_y = gTarget->GetPointY(iPnt) - parameters[1];
      double delta_r = parameters[2] - std::sqrt(delta_x * delta_x + delta_y * delta_y);
      chi2 += delta_r * delta_r;
    }
    return chi2;
  };

  // wrap chi2 function in a function object for the fit
  ROOT::Math::Functor fit_function(chi2_function, 3);
  ROOT::Fit::Fitter fitter;

  //  Set initial values and variables names
  double internal_initial_values[3] = {initial_values[0], initial_values[1], initial_values[2]};
  fitter.SetFCN(fit_function, internal_initial_values);
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("R");
  fitter.Config().ParSettings(2).SetLowerLimit(0);
  if (fix_XY)
  {
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
  }

  //  Fitting
  if (!fitter.FitFCN())
  {
    // Error("fit_circle", "Fit failed");
    //  return {{{-2., 0.}, {-2., 0.}, {-2., 0.}}};
  }
  const ROOT::Fit::FitResult &fit_result = fitter.Result();

  auto iTer = -1;
  for (auto current_parameter : fit_result.Parameters())
  {
    iTer++;
    result[iTer][0] = current_parameter;
    result[iTer][1] = fit_result.Errors()[iTer];
  }

  //  Calculate chi2
  double *test = new double[3];
  test[0] = result[0][0];
  test[1] = result[1][0];
  test[2] = result[2][0];
  auto myChi2 = chi2_function(test);

  return result;
}

template <bool fix_XY = true>
circle_fit_results
fit_circles(std::vector<TGraphErrors *> gTargets, std::array<float, 2> center_guess, std::vector<float> radiuses_guesses)
{
  circle_fit_results result;

  //  Chi2 minimisation for points in a circle
  auto chi2_function = [&](const double *parameters)
  {
    float chi2 = 0;
    auto iGraph = -1;
    for (auto current_graph : gTargets)
    {
      iGraph++;
      for (int iPnt = 0; iPnt < current_graph->GetN(); iPnt++)
      {
        double delta_x = current_graph->GetPointX(iPnt) - parameters[0];
        double delta_y = current_graph->GetPointY(iPnt) - parameters[1];
        double delta_r = parameters[2 + iGraph] - std::sqrt(delta_x * delta_x + delta_y * delta_y);
        chi2 += delta_r * delta_r;
      }
    }
    return chi2;
  };

  // wrap chi2 function in a function object for the fit
  ROOT::Math::Functor fit_function(chi2_function, 2 + gTargets.size());
  ROOT::Fit::Fitter fitter;

  //  Set initial values and variables names
  double *internal_initial_values = new double[2 + gTargets.size()];
  internal_initial_values[0] = center_guess[0];
  internal_initial_values[1] = center_guess[1];
  auto iTer = 1;
  for (auto current_radius : radiuses_guesses)
  {
    iTer++;
    internal_initial_values[iTer] = current_radius;
  }
  fitter.SetFCN(fit_function, internal_initial_values);
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");

  iTer = -1;
  for (auto current_radius : radiuses_guesses)
  {
    iTer++;
    if (iTer >= gTargets.size())
      break;
    fitter.Config().ParSettings(2 + iTer).SetName(Form("R_%i", iTer));
    fitter.Config().ParSettings(2 + iTer).SetLowerLimit(0);
  }
  if (fix_XY)
  {
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
  }

  //  Fitting
  if (!fitter.FitFCN())
  {
    // Error("fit_circle", "Fit failed");
    //  return {{{-2., 0.}, {-2., 0.}, {-2., 0.}}};
  }
  const ROOT::Fit::FitResult &fit_result = fitter.Result();

  iTer = -1;
  for (auto current_parameter : fit_result.Parameters())
  {
    iTer++;
    result[iTer][0] = current_parameter;
    result[iTer][1] = fit_result.Errors()[iTer];
  }
  return result;
}

template <bool fix_XY = true>
circle_fit_results
fit_circle_2D(TH2F *hTarget, std::array<float, 3> initial_values)
{

  // TODO

  /*
  circle_fit_results result;

  //  Chi2 minimisation for points in a circle
  auto chi2_function = [&](const double *parameters)
  {
    float chi2 = 0;
    for (int iPnt = 0; iPnt < gTarget->GetN(); iPnt++)
    {
      double delta_x = gTarget->GetPointX(iPnt) - parameters[0];
      double delta_y = gTarget->GetPointY(iPnt) - parameters[1];
      double delta_r = parameters[2] - std::sqrt(delta_x * delta_x + delta_y * delta_y);
      chi2 += delta_r * delta_r;
    }
    return chi2;
  };

  // wrap chi2 function in a function object for the fit
  ROOT::Math::Functor fit_function(chi2_function, 3);
  ROOT::Fit::Fitter fitter;

  //  Set initial values and variables names
  double internal_initial_values[3] = {initial_values[0], initial_values[1], initial_values[2]};
  fitter.SetFCN(fit_function, internal_initial_values);
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("R");
  fitter.Config().ParSettings(2).SetLowerLimit(0);
  if (fix_XY)
  {
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
  }

  //  Fitting
  if (!fitter.FitFCN())
  {
    // Error("fit_circle", "Fit failed");
    //  return {{{-2., 0.}, {-2., 0.}, {-2., 0.}}};
  }
  const ROOT::Fit::FitResult &fit_result = fitter.Result();

  auto iTer = -1;
  for (auto current_parameter : fit_result.Parameters())
  {
    iTer++;
    result[iTer][0] = current_parameter;
    result[iTer][1] = fit_result.Errors()[iTer];
  }

  //  Calculate chi2
  double *test = new double[3];
  test[0] = result[0][0];
  test[1] = result[1][0];
  test[2] = result[2][0];
  auto myChi2 = chi2_function(test);

  return result;
  */
}

circle_fit_results get_best_circle(recodata reco_data, ringdata ring_data, float radius_cut)
{
  std::array<float, 3> initial_values = {ring_data.x0[0], ring_data.y0[0], ring_data.r0[0]};
  auto target_circle = plot(select_points(reco_data, ring_data, radius_cut));
  if (target_circle->GetN() <= 3)
    return {{{-1., 0.}, {-1., 0.}, {-1., 0.}}};
  return fit_circle(target_circle, initial_values);
}

circle_fit_results get_best_circle(recodata reco_data, std::array<float, 3> ring_data, float radius_cut)
{
  auto target_circle = plot(select_points(reco_data, ring_data, radius_cut));
  if (target_circle->GetN() <= 3)
    return {{{-1., 0.}, {-1., 0.}, {-1., 0.}}};
  return fit_circle(target_circle, ring_data);
}

TCanvas *plot(recodata reco_data, ringdata ring_data, float selection_cut = 1.)
{
  auto cDrawResult = plot_standard_canvas();
  cDrawResult->SetName("cDrawResult");

  //  Full graph
  auto full_graph = plot(reco_data);
  full_graph->SetMarkerStyle(5);
  full_graph->SetMarkerColor(kBlack);
  full_graph->Draw("SAME PE");

  //  Selected graph
  auto select_graph = plot(select_points(reco_data, ring_data, selection_cut));
  select_graph->SetMarkerStyle(24);
  select_graph->Draw("SAME PE");

  //  Reference from hough
  auto ring_graph = plot(ring_data);

  //  Reference from fit
  auto ring_fit_rslt = plot(get_best_circle(reco_data, ring_data, selection_cut));

  return cDrawResult;
}

std::vector<float> get_radius_delta(circle_fit_results ring_fit_rslt, recodata reco_data)
{
  std::vector<float> result;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    auto current_dx = reco_data.x[iPnt] - ring_fit_rslt[0][0];
    auto current_dy = reco_data.y[iPnt] - ring_fit_rslt[1][0];
    result.push_back(sqrt(current_dx * current_dx + current_dy * current_dy) - ring_fit_rslt[2][0]);
  }
  return result;
}

bool good_identified_ring(ringdata ring_data, float offset_x = 2.7, float offset_y = 1.5, float sigma_x = 2.1, float sigma_y = 1.4)
{
  auto normalised_x = (ring_data.x0[0] - offset_x) / sigma_x;
  auto normalised_y = (ring_data.y0[0] - offset_y) / sigma_y;
  if (hypot(normalised_x, normalised_y) > 3.)
    return false;
  return true;
}
//  *** *** Put it in a higher level of utility *** ***
//  Functions
TF1 *fPoissonian = new TF1("fPoissonian", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 10);
//  Make array
float *make_array(std::vector<float> vTarget)
{
  float *result = new float[vTarget.size()];
  auto iTer = -1;
  for (auto element : vTarget)
  {
    iTer++;
    vTarget[iTer] = element;
  }
  return result;
}

void plot_circle(std::array<float, 3> parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1)
{
  auto result = new TEllipse(parameters[0], parameters[1], parameters[2]);
  result->SetFillStyle(0);
  result->SetLineColor(line_color);
  result->SetLineStyle(line_style);
  result->SetLineWidth(line_width);
  result->DrawEllipse(parameters[0], parameters[1], parameters[2], 0, 0, 360, 0, "same");
}

void plot_box(std::array<float, 4> parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1)
{
  auto result = new TBox(parameters[0], parameters[1], parameters[2], parameters[3]);
  result->SetFillStyle(0);
  result->SetLineColor(line_color);
  result->SetLineStyle(line_style);
  result->SetLineWidth(line_width);
  result->DrawBox(parameters[0], parameters[1], parameters[2], parameters[3]);
}

void plot_circle(std::array<std::array<float, 2>, 3> parameters, std::array<std::array<int, 3>, 3> plot_options)
{
  std::array<float, 3> reduced_parameters = {parameters[0][0], parameters[1][0], parameters[2][0]};
  plot_circle(reduced_parameters, plot_options[0][0], plot_options[0][1], plot_options[0][2]);
  reduced_parameters = {parameters[0][0], parameters[1][0], parameters[2][0] + parameters[2][1]};
  plot_circle(reduced_parameters, plot_options[1][0], plot_options[1][1], plot_options[1][2]);
  reduced_parameters = {parameters[0][0], parameters[1][0], parameters[2][0] - parameters[2][1]};
  plot_circle(reduced_parameters, plot_options[1][0], plot_options[1][1], plot_options[1][2]);
  std::array<float, 4> box_parameters = {parameters[0][0] - parameters[0][1], parameters[1][0] - parameters[1][1], parameters[0][0] + parameters[0][1], parameters[1][0] + parameters[1][1]};
  plot_box(box_parameters, plot_options[2][0], plot_options[2][1], plot_options[2][2]);
}

template <bool grid_x = true, bool grid_y = true>
TCanvas *standard_canvas(std::string name, std::string title, float nXpixels, float nYpixels)
{
  TCanvas *cResult = new TCanvas(name.c_str(), title.c_str(), nXpixels, nYpixels);
  gStyle->SetOptStat(0);
  gPad->SetGridx(grid_x);
  gPad->SetGridy(grid_y);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  return cResult;
}

template <bool grid_x = true, bool grid_y = true>
TCanvas *standard_canvas_2D(std::string name, std::string title, float nXpixels = 1145, float nYpixels = 1000)
{
  TCanvas *cResult = new TCanvas(name.c_str(), title.c_str(), nXpixels, nYpixels);
  gStyle->SetOptStat(0);
  gPad->SetGridx(grid_x);
  gPad->SetGridy(grid_y);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  return cResult;
}

void plot_prediction_for_run(std::string run_tag, float x_coordinate = 0., float y_coordinate = 0.)
{
  if (!database::is_run_recorded(run_tag))
  {
    cout << "[ERROR] Prediction plot for run " << run_tag << " is undefined, run info is not stored!" << endl;
    return;
  }
  auto beam_polarity = database::get_beam_polarity(run_tag).second;
  if (beam_polarity == 0)
  {
    cout << "[WARNING] Polarity for run " << run_tag << " is undefined, will suppose it is negative beam" << endl;
    beam_polarity = -1;
  }
  bool has_protons = true; //(beam_polarity > 0);
  bool has_aerogel = database::get_has_aerogel(run_tag).second;
  bool has_gas = true; // database::get_has_gas(run_tag).second;
  auto beam_energy = database::get_beam_energy(run_tag).second;

  //  Prediction for Aerogel rings
  if (has_aerogel)
  {
    auto n_aerogel = database::get_n_aerogel(run_tag).second;
    auto arm_length_aerogel = 373 + database::get_aerogel_mirror_z(run_tag).second - 75;

    //  Prediction radius
    auto predicted_radius_aerogel_pions = get_radius(particle_mass["pion+-"], beam_energy, n_aerogel, arm_length_aerogel);
    auto predicted_radius_aerogel_kaons = get_radius(particle_mass["kaon+-"], beam_energy, n_aerogel, arm_length_aerogel);
    auto predicted_radius_aerogel_protons = get_radius(particle_mass["proton"], beam_energy, n_aerogel, arm_length_aerogel);

    // Plot
    plot_circle({x_coordinate, y_coordinate, predicted_radius_aerogel_pions}, kMagenta - 3, 10, 2);
    plot_circle({x_coordinate, y_coordinate, predicted_radius_aerogel_kaons}, kMagenta - 3, 9, 2);
    if (has_protons)
      plot_circle({x_coordinate, y_coordinate, predicted_radius_aerogel_protons}, kMagenta - 3, 5, 2);
  }

  //  Prediction for Gas rings
  if (has_gas)
  {
    //  Prediction for Aerogel rings
    auto n_gas = database::get_n_gas(run_tag).second;
    auto arm_length_gas = 1272 + database::get_aerogel_mirror_z(run_tag).second;

    //  Prediction radius
    auto predicted_radius_gas_pions = get_radius(particle_mass["pion+-"], beam_energy, n_gas, arm_length_gas);
    auto predicted_radius_gas_kaons = get_radius(particle_mass["kaon+-"], beam_energy, n_gas, arm_length_gas);
    auto predicted_radius_gas_protons = get_radius(particle_mass["proton"], beam_energy, n_gas, arm_length_gas);

    // Plot
    plot_circle({x_coordinate, y_coordinate, predicted_radius_gas_pions}, kOrange + 3, 10, 2);
    plot_circle({x_coordinate, y_coordinate, predicted_radius_gas_kaons}, kOrange + 3, 9, 2);
    if (has_protons)
      plot_circle({x_coordinate, y_coordinate, predicted_radius_gas_protons}, kOrange + 3, 5, 2);
  }
}

std::array<std::array<float, 2>, 4> get_prediction(TF1 *fFull)
{
  std::array<std::array<float, 2>, 4> result;

  //  X_{0} prediction
  result[0][0] = -fFull->GetParameter(1) * TMath::Sin(fFull->GetParameter(2));
  result[0][1] = TMath::Sqrt(TMath::Sq(fFull->GetParError(1) * TMath::Sin(fFull->GetParameter(2))) + TMath::Sq(fFull->GetParError(2) * fFull->GetParameter(1) * TMath::Cos(fFull->GetParameter(2))));

  //  Y_{0} prediction
  result[1][0] = fFull->GetParameter(1) * TMath::Cos(fFull->GetParameter(2));
  result[1][1] = TMath::Sqrt(TMath::Sq(fFull->GetParError(1) * TMath::Cos(fFull->GetParameter(2))) + TMath::Sq(fFull->GetParError(2) * fFull->GetParameter(1) * TMath::Sin(fFull->GetParameter(2))));

  //  Radius prediction
  result[2][0] = fFull->GetParameter(0);
  result[2][1] = fFull->GetParError(0);

  //  SPR prediction
  result[3][0] = fFull->GetParameter(4);
  result[3][1] = fFull->GetParError(4);

  return result;
}

std::array<std::array<float, 2>, 4> get_run_coordinates(TH2F *hEfficiency, TH2F *hPersistance2D, std::string plot_save_path = "")
{
  //  Create a local copy
  auto hTarget = (TH2F *)(hEfficiency->Clone("hTarget"));

  //  Plot aid
  TCanvas *cCheckResults = new TCanvas("cCheckResults", "cCheckResults", 1500, 500);
  cCheckResults->Divide(3, 1);
  gStyle->SetOptStat(0);
  TLatex *lLatex = new TLatex();
  TCanvas *dump = new TCanvas("dump", "dump", 500, 500);

  //  Find guess for average radius and spread
  //  --- Target projection
  auto hTarget_Yprojection = hTarget->ProjectionY("hTarget_Yprojection", -1, -1);
  //  --- Functions
  TF1 *full_first_guesstimate_function = new TF1("full_first_guesstimate_function", "[3]+[4]*x+[5]*x*x+[0]*exp(-0.5*((x-[1])/[2])**2)", 0, 1000);
  full_first_guesstimate_function->SetNpx(10000);
  full_first_guesstimate_function->SetLineColor(kBlue);
  full_first_guesstimate_function->SetLineStyle(kSolid);
  TF1 *gaus_first_guesstimate_function = new TF1("gaus_first_guesstimate_function", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0, 1000);
  gaus_first_guesstimate_function->SetNpx(10000);
  gaus_first_guesstimate_function->SetLineColor(kRed);
  gaus_first_guesstimate_function->SetLineStyle(kSolid);
  TF1 *pol2_first_guesstimate_function = new TF1("pol2_first_guesstimate_function", "[0]+[1]*x+[2]*x*x", 0, 1000);
  pol2_first_guesstimate_function->SetNpx(10000);
  pol2_first_guesstimate_function->SetLineColor(kBlue);
  pol2_first_guesstimate_function->SetLineStyle(kDashed);

  full_first_guesstimate_function->SetParameters(hTarget_Yprojection->GetEntries() * 0.1, hTarget_Yprojection->GetMean(), hTarget_Yprojection->GetRMS(), 0.1, 0.1, 0.1);
  hTarget_Yprojection->Fit(full_first_guesstimate_function);
  gaus_first_guesstimate_function->SetParameters(full_first_guesstimate_function->GetParameter(0), full_first_guesstimate_function->GetParameter(1), full_first_guesstimate_function->GetParameter(2));
  pol2_first_guesstimate_function->SetParameters(full_first_guesstimate_function->GetParameter(3), full_first_guesstimate_function->GetParameter(4), full_first_guesstimate_function->GetParameter(5));

  auto normalisation_guesstimate = full_first_guesstimate_function->GetParameter(0);
  auto main_radius_guesstimate = full_first_guesstimate_function->GetParameter(1);
  auto main_radius_spread_guesstimate = full_first_guesstimate_function->GetParameter(2);
  //  --- Canvas
  TCanvas *check_guesstimate = standard_canvas("check_guesstimate", "", 1000, 1000);
  gPad->SetLogy();
  hTarget_Yprojection->Draw("");
  gaus_first_guesstimate_function->DrawCopy("SAME");
  pol2_first_guesstimate_function->DrawCopy("SAME");

  //  --- full function
  TF2 *full_fit = new TF2("full_fit", "[5]+[3]*exp(-0.5*((y-([0]+[1]*TMath::Sin(x-[2])))/[4])**2)", -1000, 1000, -1000, 1000);
  full_fit->SetNpx(1.e7);
  full_fit->SetParLimits(0, 0., 100.);
  full_fit->SetParLimits(1, 0., 100.);
  full_fit->SetParLimits(2, -TMath::Pi(), TMath::Pi());
  full_fit->SetParLimits(3, 0., 100000.);
  full_fit->SetParLimits(4, 0., 10.);
  full_fit->SetParLimits(5, 0., 1000.);
  full_fit->SetLineColor(kRed);
  full_fit->SetLineWidth(2);

  // Assign guesstimates to full fit
  cCheckResults->cd(3);
  full_fit->SetParameter(0, main_radius_guesstimate);
  full_fit->SetParameter(1, main_radius_spread_guesstimate);
  full_fit->SetParameter(2, 0.);
  full_fit->FixParameter(3, normalisation_guesstimate);
  full_fit->FixParameter(4, main_radius_spread_guesstimate);
  full_fit->FixParameter(5, 0.);
  hTarget->Fit(full_fit);
  full_fit->ReleaseParameter(3);
  full_fit->ReleaseParameter(4);
  full_fit->ReleaseParameter(5);
  full_fit->SetParLimits(3, 0., 100000.);
  full_fit->SetParLimits(4, 0., 10.);
  full_fit->SetParLimits(5, 0., 1000.);
  hTarget->Fit(full_fit);

  auto final_prediction = get_prediction(full_fit);
  TCanvas *cDrawResult = standard_canvas_2D("cDrawResult", "cDrawResult");
  hPersistance2D->Draw("COLZ");
  std::array<std::array<int, 3>, 3> plot_options = {{{kBlack, kSolid, 3}, {kRed, kDashed, 2}, {kRed, kSolid, 2}}};
  std::array<std::array<float, 2>, 3> plot_coordinates = {{{final_prediction[0][0], final_prediction[0][1]}, {final_prediction[1][0], final_prediction[1][1]}, {final_prediction[2][0], final_prediction[3][0]}}};
  plot_circle(plot_coordinates, plot_options);

  lLatex->SetTextSize(0.03);
  lLatex->DrawLatexNDC(0.02, 0.965, Form("X_{0} : %.3f#pm%.3f mm", final_prediction[0][0], final_prediction[0][1]));
  lLatex->DrawLatexNDC(0.27, 0.965, Form("Y_{0} : %.3f#pm%.3f mm", final_prediction[1][0], final_prediction[1][1]));
  lLatex->DrawLatexNDC(0.52, 0.965, Form("R : %.3f#pm%.3f mm", final_prediction[2][0], final_prediction[2][1]));
  lLatex->DrawLatexNDC(0.77, 0.965, Form("#sigma_{R} : %.3f#pm%.3f mm", full_fit->GetParameter(4), full_fit->GetParError(4)));

  if (!plot_save_path.empty())
    cDrawResult->SaveAs(plot_save_path.c_str());

  return final_prediction;

  /*
    //   Define functions
    //  --- full function
    //  --- sin function
    TF1 *Sin = new TF1("Sin", "[0]+[1]*TMath::Sin(x-[2])", -1000, 1000);
    Sin->SetNpx(1.e7);
    Sin->SetParLimits(0, 0., 1000.);
    Sin->SetParLimits(1, 0., 1000.);
    Sin->SetParLimits(2, -TMath::Pi(), TMath::Pi());
    Sin->SetLineColor(kRed);
    Sin->SetLineWidth(2);
    //  --- gauss function
    TF1 *Gauss = new TF1("Gauss", "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", -1000, 1000);
    Gauss->SetParLimits(0, 0., 100000.);
    Gauss->SetParLimits(1, 0., 100.);
    Gauss->SetParLimits(2, 0., 10.);
    Gauss->SetParLimits(3, 0., 10000.);
    Gauss->SetNpx(1.e7);
    Gauss->SetLineColor(kRed);
    Gauss->SetLineWidth(2);

    //  Fit slices for gaus guesstimates
    auto select_slice_xBin = 0;
    auto select_slice_entries = 0;
    TGraphErrors *g_gaus_average = new TGraphErrors();
    TGraphErrors *g_gaus_stdv = new TGraphErrors();
    //  --- Fit all slices to rejet bkg
    for (auto xBin = 1; xBin <= hTarget->GetNbinsX(); xBin++)
    {
      auto fit_slice = hTarget->ProjectionY(Form("current_slice_%i", xBin), xBin, xBin);
      //  --- Skip empty slices
      if (fit_slice->GetEntries() < 100)
        continue;
      auto bin_center = hTarget->GetXaxis()->GetBinCenter(xBin);
      //  --- Prepare the fit
      auto gaus_mean = fit_slice->GetMean();
      auto gaus_stdv = fit_slice->GetRMS() / (TMath::Sqrt(12));
      auto gaus_peak = fit_slice->GetBinContent(fit_slice->GetXaxis()->FindBin(gaus_mean));
      auto _bkg_mean = fit_slice->GetBinContent(fit_slice->GetXaxis()->FindBin(gaus_mean - 5 * gaus_stdv));
      //  --- Skip
      auto current_point = g_gaus_average->GetN();
      if ((xBin > 1) && (fabs(g_gaus_average->GetPointY(current_point - 1) - gaus_mean) > 3 * g_gaus_stdv->GetPointY(current_point - 1)))
        continue;
      Gauss->SetParameters(gaus_peak, gaus_mean, gaus_stdv, _bkg_mean);
      dump->cd();
      fit_slice->Fit(Gauss, "IMRESQ");
      g_gaus_average->SetPoint(current_point, bin_center, Gauss->GetParameter(1));
      g_gaus_average->SetPointError(current_point, 0., Gauss->GetParError(1));
      g_gaus_stdv->SetPoint(current_point, bin_center, Gauss->GetParameter(2));
      g_gaus_stdv->SetPointError(current_point, 0., Gauss->GetParError(2));
      break;
    }

    // Assign guesstimates to full fit
    cCheckResults->cd(3);
    full_fit->SetParameter(0, Gauss->GetParameter(1) 70.);
    full_fit->SetParameter(1, g_gaus_average->GetRMS(2));
    full_fit->SetParameter(2, 0.);
    full_fit->FixParameter(3, Gauss->GetParameter(0));
    full_fit->FixParameter(4, Gauss->GetParameter(2));
    full_fit->FixParameter(5, Gauss->GetParameter(3));
    hTarget->Fit(full_fit);
    full_fit->ReleaseParameter(3);
    full_fit->ReleaseParameter(4);
    full_fit->ReleaseParameter(5);
    full_fit->SetParLimits(3, 0., 100000.);
    full_fit->SetParLimits(4, 0., 10.);
    full_fit->SetParLimits(5, 0., 1000.);
    hTarget->Fit(full_fit);

    //  Check the sin fit
    cCheckResults->cd(2);
    g_gaus_average->Draw();
    Sin->SetParameters(full_fit->GetParameter(0), full_fit->GetParameter(1), full_fit->GetParameter(2));
    Sin->SetLineColor(kRed);
    Sin->SetLineStyle(kDashed);
    Sin->SetLineWidth(2);
    Sin->DrawCopy("SAMEL");

    auto final_prediction = get_prediction(full_fit);
    TCanvas *cDrawResult = standard_canvas_2D("cDrawResult", "cDrawResult");
    hPersistance2D->Draw("COLZ");
    std::array<std::array<int, 3>, 3> plot_options = {{{kBlack, kSolid, 5}, {kRed, kDashed, 3}, {kRed, kSolid, 3}}};
    std::array<std::array<float, 2>, 3> plot_coordinates = {{{final_prediction[0][0], final_prediction[0][1]}, {final_prediction[1][0], final_prediction[1][1]}, {final_prediction[2][0], final_prediction[2][1]}}};
    plot_circle(plot_coordinates, plot_options);

    lLatex->SetTextSize(0.03);
    lLatex->DrawLatexNDC(0.02, 0.965, Form("X_{0} : %.3f#pm%.3f mm", final_prediction[0][0], final_prediction[0][1]));
    lLatex->DrawLatexNDC(0.27, 0.965, Form("Y_{0} : %.3f#pm%.3f mm", final_prediction[1][0], final_prediction[1][1]));
    lLatex->DrawLatexNDC(0.52, 0.965, Form("R : %.3f#pm%.3f mm", final_prediction[2][0], final_prediction[2][1]));
    lLatex->DrawLatexNDC(0.77, 0.965, Form("#sigma_{R} : %.3f#pm%.3f mm", full_fit->GetParameter(4), full_fit->GetParError(4)));

    if (!plot_save_path.empty())
      cDrawResult->SaveAs(plot_save_path.c_str());

    return final_prediction;
    */
}

//  Coordinates utility
inline std::array<float, 2> polar_to_cartesian(std::array<float, 2> target, std::array<float, 2> center_shift = {0., 0.})
{
  float r_new_coordinate = target[0] - center_shift[0];
  float phi_new_coordinate = target[1] - center_shift[1];
  float x_coordinate = target[0] * TMath::Cos(target[1]);
  float y_coordinate = target[0] * TMath::Sin(target[1]);
  return {x_coordinate, y_coordinate};
}
inline std::array<float, 2> cartesian_to_polar(std::array<float, 2> target, std::array<float, 2> center_shift = {0., 0.})
{
  float x_new_coordinate = target[0] - center_shift[0];
  float y_new_coordinate = target[1] - center_shift[1];
  float r_coordinate = TMath::Sqrt(x_new_coordinate * x_new_coordinate + y_new_coordinate * y_new_coordinate);
  float phi_coordinate = TMath::ATan2(y_new_coordinate, x_new_coordinate);
  return {r_coordinate, phi_coordinate};
}
bool is_within_ring(std::array<float, 2> cartesian_coordinates, std::array<float, 3> ring_parameters, float ring_radius_sigma = 0.)
{
  //  Improvements:
  //  Make X and Y have different possible values
  //  Check values received for the interval
  bool is_within_ring = false;
  auto target_radius = cartesian_to_polar(cartesian_coordinates, {ring_parameters[0], ring_parameters[1]})[0];
  if ((target_radius >= ring_parameters[2] - ring_radius_sigma) && (target_radius <= ring_parameters[2] + ring_radius_sigma))
    is_within_ring = true;
  return is_within_ring;
}
bool is_within_SiPM(std::array<float, 2> cartesian_coordinates, std::array<float, 2> SiPM_center)
{
  bool is_within_SiPM = false;
  auto x_distance = fabs(cartesian_coordinates[0] - SiPM_center[0]);
  auto y_distance = fabs(cartesian_coordinates[1] - SiPM_center[1]);
  if ((x_distance) <= kSiPM_x_dimension && fabs(y_distance) <= kSiPM_y_dimension)
    is_within_SiPM = true;
  return is_within_SiPM;
}

//  Fill functions
template <bool rnd_smearing = true>
void fill_persistance(TH2F *hTarget, recodata reco_data)
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    hTarget->Fill(rnd_smearing ? gRandom->Uniform(reco_data.x[iPnt] - 1.5, reco_data.x[iPnt] + 1.5) : reco_data.x[iPnt], rnd_smearing ? gRandom->Uniform(reco_data.y[iPnt] - 1.5, reco_data.y[iPnt] + 1.5) : reco_data.y[iPnt]);
}
void fill_time(TH1F *hTarget, recodata reco_data)
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    hTarget->Fill(reco_data.t[iPnt]);
}
void fill_persistance(TH3F *hTarget, recodata reco_data)
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    hTarget->Fill(gRandom->Uniform(reco_data.x[iPnt] - 1.5, reco_data.x[iPnt] + 1.5), gRandom->Uniform(reco_data.y[iPnt] - 1.5, reco_data.y[iPnt] + 1.5), reco_data.t[iPnt]);
}
void fill_polar_efficiency(recodata reco_data, TH2F *hTarget, std::array<float, 2> center_shift = {0., 0.})
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    float new_x = gRandom->Uniform(reco_data.x[iPnt] - 1.5, reco_data.x[iPnt] + 1.5);
    float new_y = gRandom->Uniform(reco_data.y[iPnt] - 1.5, reco_data.y[iPnt] + 1.5);
    auto new_coordinates = cartesian_to_polar({new_x, new_y}, center_shift);
    hTarget->Fill(new_coordinates[1], new_coordinates[0]);
  }
}

float delta_phi(float phi_1, float phi_2, float min_phi = -TMath::Pi(), float max_phi = TMath::Pi())
{
  auto delta_phi = phi_2 - phi_1;
  if (delta_phi < min_phi)
    delta_phi += 2 * TMath::Pi();
  if (delta_phi > max_phi)
    delta_phi -= 2 * TMath::Pi();
  return delta_phi;
}

void fill_correlations(recodata reco_data, TH2F *hTarget)
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    for (int jPnt = 0; jPnt < reco_data.n; jPnt++)
    {
      if (iPnt == jPnt)
        continue;
      auto current_iPhi = TMath::ATan2(reco_data.y[iPnt], reco_data.x[iPnt]);
      auto current_jPhi = TMath::ATan2(reco_data.y[jPnt], reco_data.x[jPnt]);
      hTarget->Fill(reco_data.n, delta_phi(current_iPhi, current_jPhi));
    }
  }
}

void fill_with_ring_coverage(TH2F *target, std::array<float, 3> ring_parameters, float ring_radius_sigma = 0.)
{
  for (auto xBin = 0; xBin < target->GetNbinsX(); xBin++)
  {
    float current_bin_center_x = (float)(target->GetXaxis()->GetBinCenter(xBin));
    for (auto yBin = 0; yBin < target->GetNbinsY(); yBin++)
    {
      float current_bin_center_y = float(target->GetYaxis()->GetBinCenter(yBin));
      if (is_within_ring({current_bin_center_x, current_bin_center_y}, ring_parameters, ring_radius_sigma))
      {
        target->SetBinContent(xBin, yBin, 1);
      }
      else
      {
        target->SetBinContent(xBin, yBin, 0);
      }
    }
  }
}

template <bool allow_overlap = false>
void fill_with_SiPM_coverage(TH2F *target, std::array<float, 2> SiPM_center)
{
  float x_center_bin = target->GetXaxis()->FindBin(SiPM_center[0]);
  float y_center_bin = target->GetYaxis()->FindBin(SiPM_center[1]);
  float x_center_bin_center = target->GetXaxis()->GetBinCenter(x_center_bin);
  float y_center_bin_center = target->GetYaxis()->GetBinCenter(y_center_bin);
  float current_center_bin_content = target->GetBinContent(target->FindBin(x_center_bin, y_center_bin));

  //  Did we set this SiPM already?
  //  * No overlap is permitted
  if ((current_center_bin_content > 0) && !allow_overlap)
    return;

  //  Start with know center bin
  for (auto xBin = 0; xBin < target->GetNbinsX(); xBin++)
  {
    //  Check at least one fill has been done
    bool at_least_one_fill = false;

    float current_bin_center_x_high = (float)(target->GetXaxis()->GetBinCenter(x_center_bin + xBin));
    float current_bin_center_x_low_ = (float)(target->GetXaxis()->GetBinCenter(x_center_bin - xBin));
    for (auto yBin = 0; yBin < target->GetNbinsY(); yBin++)
    {
      float current_bin_center_y_high = (float)(target->GetYaxis()->GetBinCenter(y_center_bin + yBin));
      float current_bin_center_y_low_ = (float)(target->GetYaxis()->GetBinCenter(y_center_bin - yBin));

      //  Set bin contents
      if (is_within_SiPM({current_bin_center_x_high, current_bin_center_y_high}, SiPM_center))
      {
        target->SetBinContent(x_center_bin + xBin, y_center_bin + yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_high, current_bin_center_y_low_}, SiPM_center))
      {
        target->SetBinContent(x_center_bin + xBin, y_center_bin - yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_low_, current_bin_center_y_high}, SiPM_center))
      {
        target->SetBinContent(x_center_bin - xBin, y_center_bin + yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_low_, current_bin_center_y_low_}, SiPM_center))
      {
        target->SetBinContent(x_center_bin - xBin, y_center_bin - yBin, 1);
        at_least_one_fill = true;
      }

      //  return if we stopped filling
      if (!at_least_one_fill)
        break;
    }
    if (!is_within_SiPM({current_bin_center_x_high, y_center_bin_center}, SiPM_center) && !is_within_SiPM({current_bin_center_x_low_, y_center_bin_center}, SiPM_center))
      return;
  }
}

float calculate_efficiency(TH2F *sipm_coverage_map, std::array<float, 3> ring_parameters, float ring_radius_sigma, TH2F *&ring_coverage_map_SiPM_filtered, TH2F *&ring_coverage_map)
{
  ring_coverage_map = (TH2F *)(sipm_coverage_map->Clone("ring_coverage_map"));
  ring_coverage_map_SiPM_filtered = (TH2F *)(sipm_coverage_map->Clone("ring_coverage_map_SiPM_filtered"));
  ring_coverage_map->Clear();
  fill_with_ring_coverage(ring_coverage_map, ring_parameters, ring_radius_sigma);
  ring_coverage_map_SiPM_filtered->Multiply(ring_coverage_map, sipm_coverage_map);
  return ring_coverage_map_SiPM_filtered->Integral(-1., -1., -1., -1.) / ring_coverage_map->Integral(-1., -1., -1., -1.);
}

template <bool return_selected = true>
recodata select_points(recodata reco_data, std::array<float, 2> square_center, std::array<float, 2> xy_square_cut)
{
  recodata selected_reco_data;
  recodata discarded_reco_data;
  selected_reco_data.n = 0;
  discarded_reco_data.n = 0;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    auto x_shift = reco_data.x[iPnt] - square_center[0];
    auto y_shift = reco_data.y[iPnt] - square_center[1];
    if ((fabs(x_shift) > xy_square_cut[0]) || (fabs(y_shift) > xy_square_cut[1]))
    {
      discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
      discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
      discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
      discarded_reco_data.n++;
      continue;
    }
    selected_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
    selected_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
    selected_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
    selected_reco_data.n++;
  }
  return return_selected ? selected_reco_data : discarded_reco_data;
}

template <int coordinate = 3,
          bool return_selected = true>
recodata select_points(recodata reco_data, std::array<float, 2> center_and_cut)
{
  recodata selected_reco_data;
  recodata discarded_reco_data;
  selected_reco_data.n = 0;
  discarded_reco_data.n = 0;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    switch (coordinate)
    {
    case 1:
      if (fabs(reco_data.x[iPnt] - center_and_cut[0]) < center_and_cut[1])
      {
        discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
        discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
        discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
        discarded_reco_data.n++;
        continue;
      }
      break;
    case 2:
      if (fabs(reco_data.y[iPnt] - center_and_cut[0]) < center_and_cut[1])
      {
        discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
        discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
        discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
        discarded_reco_data.n++;
        continue;
      }
      break;
    case 3:
      if (fabs(reco_data.t[iPnt] - center_and_cut[0]) < center_and_cut[1])
      {
        discarded_reco_data.x[discarded_reco_data.n] = reco_data.x[iPnt];
        discarded_reco_data.y[discarded_reco_data.n] = reco_data.y[iPnt];
        discarded_reco_data.t[discarded_reco_data.n] = reco_data.t[iPnt];
        discarded_reco_data.n++;
        continue;
      }
      break;
    default:
      break;
    }
    selected_reco_data.x[selected_reco_data.n] = reco_data.x[iPnt];
    selected_reco_data.y[selected_reco_data.n] = reco_data.y[iPnt];
    selected_reco_data.t[selected_reco_data.n] = reco_data.t[iPnt];
    selected_reco_data.n++;
  }
  return return_selected ? selected_reco_data : discarded_reco_data;
}

bool fill_coordinate_x(TH1F *target_histo, recodata reco_data, float min_val = std::numeric_limits<double>::min(), float max_val = std::numeric_limits<double>::max())
{
  bool filled = false;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    if ((reco_data.y[iPnt] >= min_val) && (reco_data.y[iPnt] <= max_val))
    {
      target_histo->Fill(reco_data.x[iPnt]);
      filled = true;
    }
  }
  return filled;
}

bool fill_coordinate_y(TH1F *target_histo, recodata reco_data, float min_val = std::numeric_limits<double>::min(), float max_val = std::numeric_limits<double>::max())
{
  bool filled = false;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    if ((reco_data.x[iPnt] >= min_val) && (reco_data.x[iPnt] <= max_val))
    {
      target_histo->Fill(reco_data.y[iPnt]);
      filled = true;
    }
  }
  return filled;
}

bool fill_coordinate_R(TH1F *target_histo, recodata reco_data, std::array<float, 2> center_shift = {0., 0.}, float radius_offset = 0.)
{
  bool filled = false;
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
  {
    auto polar_coordinates = cartesian_to_polar({reco_data.x[iPnt], reco_data.y[iPnt]}, center_shift);
    target_histo->Fill(polar_coordinates[0] - radius_offset);
    filled = true;
  }
  return filled;
}

std::vector<std::tuple<float, float, int>> 
find_peaks(TH1F *original_histo, int bin_span = 10, float cutoff_threshold = 0.33)
{
  //  Result
  std::vector<std::tuple<float, float, int>> result;

  //  Clone target to manipulate
  auto target_histo = (TH1F *)(original_histo->Clone("tmp_get_filtered_center"));
  auto target_maximum = target_histo->GetMaximum();

  // Loop over bins
  for (auto iBin = 1 + bin_span; iBin <= target_histo->GetNbinsX() - bin_span; iBin++)
  {
    auto neighbourhood_maximum = 0;
    auto current_bin = target_histo->GetBinContent(iBin);
    for (auto iSpan = 1; iSpan <= bin_span; iSpan++)
    {
      auto current_bin_left = target_histo->GetBinContent(iBin - iSpan);
      auto current_bin_right = target_histo->GetBinContent(iBin + iSpan);
      if ((neighbourhood_maximum < current_bin_left) || (neighbourhood_maximum < current_bin_right))
        neighbourhood_maximum = max(current_bin_left, current_bin_right);
      if (neighbourhood_maximum > current_bin)
        break;
    }
    if ((neighbourhood_maximum < current_bin) && (current_bin > cutoff_threshold * target_maximum))
      result.push_back({target_histo->GetBinCenter(iBin), target_histo->GetBinContent(iBin), iBin});
  }

  delete target_histo;
  return result;
}

std::vector<std::vector<std::array<float, 2>>>
find_rings(TH1F *original_histo_X, TH1F *original_histo_Y, float x_center = -1.5, float y_center = -1.5)
{
  //  Clone target to manipulate
  auto target_histo_X = (TH1F *)(original_histo_X->Clone("tmp_find_rings_X"));
  auto target_histo_Y = (TH1F *)(original_histo_Y->Clone("tmp_find_rings_Y"));

  //  Find peaks
  auto current_peaks_X = find_peaks(target_histo_X);
  auto current_peaks_Y = find_peaks(target_histo_Y);

  //  ==  Assumption is made that the beam is within the smaller circle
  auto n_cirles = 0;
  auto current_circle = 0;
  std::vector<std::vector<std::array<float, 2>>> circles;

  //  First loop on x_findings
  for (auto current_peak : current_peaks_X)
  {
    auto current_X = get<0>(current_peak);
    if (current_X < 0)
    {
      n_cirles++;
      circles.push_back({{current_X, y_center}});
      current_circle++;
    }
    else
    {
      current_circle--;
      if (current_circle < 0)
      {
        current_circle = n_cirles;
        n_cirles++;
        circles.insert(circles.begin(), {{current_X, y_center}});
        continue;
      }
      circles[current_circle].push_back({current_X, y_center});
    }
  }

  //  Second loop on y_findings
  auto n_left_rings = 0;
  for (auto current_peak : current_peaks_Y)
  {
    auto current_Y = get<0>(current_peak);
    if (current_Y < 0)
      n_left_rings++;
  }
  if (n_left_rings == n_cirles)
  {
    current_circle = 0;
    for (auto current_peak : current_peaks_Y)
    {
      auto current_Y = get<0>(current_peak);
      if (current_Y < 0)
      {
        circles[current_circle].push_back({x_center, current_Y});
        current_circle++;
      }
      else
      {
        current_circle--;
        if (current_circle < 0)
        {
          current_circle = n_cirles;
          n_cirles++;
          circles.insert(circles.begin(), {{x_center, current_Y}});
          continue;
        }
        circles[current_circle].push_back({x_center, current_Y});
      }
    }
  }

  std::reverse(circles.begin(), circles.end());
  return circles;
}

std::vector<std::vector<std::array<float, 2>>>
merge_circles(std::vector<std::vector<std::vector<std::array<float, 2>>>> target_circles)
{
  std::vector<std::vector<std::array<float, 2>>> final_circles;
  for (auto current_circle_array : target_circles)
  {
    auto icircle = -1;
    for (auto current_circle : current_circle_array)
    {
      icircle++;
      if (final_circles.size() < icircle + 1)
        final_circles.push_back({});
      for (auto current_point : current_circle)
      {
        final_circles[icircle].push_back(current_point);
      }
    }
  }
  return final_circles;
}

template <bool simultaneous_fit = true>
std::vector<circle_fit_results> fit_multiple_rings(std::vector<std::array<TH1F *, 2>> original_histo_s, std::vector<std::array<float, 2>> centers)
{
  std::vector<std::vector<std::vector<std::array<float, 2>>>> final_circles_array;
  auto iTer = -1;
  for (auto current_xy_pair : original_histo_s)
  {
    iTer++;
    final_circles_array.push_back(find_rings(current_xy_pair[0], current_xy_pair[1], centers[iTer][0], centers[iTer][1]));
  }

  auto final_circles = merge_circles(final_circles_array);

  std::vector<circle_fit_results> final_rings;
  std::vector<TGraphErrors *> gfinal_rings;
  auto current_n_circle = -1;
  for (auto current_circle : final_circles)
  {
    current_n_circle++;
    gfinal_rings.push_back(new TGraphErrors());
    for (auto current_point : current_circle)
    {
      gfinal_rings[current_n_circle]->SetPoint(gfinal_rings[current_n_circle]->GetN(), current_point[0], current_point[1]);
    }
    auto circle_fit = fit_circle<false>(gfinal_rings[current_n_circle], {0, 0, 20});
    final_rings.push_back(circle_fit);
  }

  if (simultaneous_fit)
  {
    auto fit_results = fit_circles<false>(gfinal_rings, {0, 0}, {100, 100, 100, 100, 100});
    auto iRing = -1;
    for (auto &current_ring : final_rings)
    {
      iRing++;
      current_ring[0][0] = fit_results[0][0];
      current_ring[1][0] = fit_results[1][0];
      current_ring[0][1] = fit_results[0][1];
      current_ring[1][1] = fit_results[1][1];
      current_ring[2][0] = fit_results[2 + iRing][0];
      current_ring[2][1] = fit_results[2 + iRing][1];
    }
  }

  return final_rings;
}

TCanvas *plot_check_coordinates(TH2F *hPersistance2D, std::vector<circle_fit_results> found_rings)
{
  TCanvas *cDrawResult = standard_canvas_2D("cDrawResult", "cDrawResult");
  TLatex *lLatex = new TLatex();
  hPersistance2D->Draw("COLZ");
  auto iRing = -1;
  for (auto current_ring : found_rings)
  {
    iRing++;
    std::array<std::array<int, 3>, 3> plot_options = {{{kBlack, kSolid, 3}, {kRed, kDashed, 2}, {kRed, kSolid, 2}}};
    std::array<std::array<float, 2>, 3> plot_coordinates = {{{current_ring[0][0], current_ring[0][1]}, {current_ring[1][0], current_ring[1][1]}, {current_ring[2][0], current_ring[2][1]}}};
    plot_circle(plot_coordinates, plot_options);

    lLatex->SetTextSize(0.03);
    lLatex->DrawLatexNDC(0.02, 0.090 - iRing * 0.032, Form("X_{0} : %.2f#pm%.2f mm", current_ring[0][0], current_ring[0][1]));
    lLatex->DrawLatexNDC(0.25, 0.090 - iRing * 0.032, Form("Y_{0} : %.2f#pm%.2f mm", current_ring[1][0], current_ring[1][1]));
    lLatex->DrawLatexNDC(0.48, 0.090 - iRing * 0.032, Form("R_%i : %.2f#pm%.2f mm", iRing, current_ring[2][0], current_ring[2][1]));
  }
  return cDrawResult;
}