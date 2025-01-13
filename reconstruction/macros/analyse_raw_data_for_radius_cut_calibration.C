#pragma once

#include "../lib/Utility.h"

TF1 *fit_function = new TF1("fit_function", "TMath::Min(TMath::Max(0.,[0]*(x-[1])+[2]),TMath::Max(0.,-[0]*(x-[1])+[2]))+gaus(3)");

void AnalyseResolution(std::string current_TH3_target, std::map<std::string, TH3F *> &_TH3F, std::map<std::string, TGraphErrors *> &_TGraphErrors, std::map<std::string, TH1D *> &_TH1D, TF1 *target_fit, float min, float max)
{
  //  Loop over the cut slices
  auto x_bins = _TH3F[current_TH3_target]->GetNbinsX();
  auto y_bins = _TH3F[current_TH3_target]->GetNbinsY();
  _TGraphErrors[current_TH3_target + "_Full_from_fit"] = new TGraphErrors();
  for (auto xBin = 1; xBin <= x_bins; xBin++)
  {
    _TH3F[current_TH3_target]->GetXaxis()->SetRange(xBin, xBin);
    auto xBin_value = _TH3F[current_TH3_target]->GetXaxis()->GetBinCenter(xBin);
    auto current_graph = current_TH3_target + std::string(Form("_graph_x%i_Rcut_%f", xBin, xBin_value));
    _TGraphErrors[current_graph] = new TGraphErrors();
    for (auto yBin = 1; yBin <= y_bins; yBin++)
    {
      _TH3F[current_TH3_target]->GetYaxis()->SetRange(yBin, yBin);
      auto yBin_value = _TH3F[current_TH3_target]->GetYaxis()->GetBinCenter(yBin);
      auto current_target_name = current_TH3_target + std::string(Form("current_slice_x%i_Rcut_%f_y%i_Nph_%f", xBin, xBin_value, yBin, yBin_value));
      _TH1D[current_target_name] = (TH1D *)((TH1D *)(_TH3F[current_TH3_target]->Project3D("Z"))->Clone(current_target_name.c_str()));
      if (_TH1D[current_target_name]->GetEntries() < 100)
        continue;
      _TH1D[current_target_name]->Rebin(2);
      TF1 *gaus_fnc = new TF1("gaus_fnc", "[0]*exp(-0.5*((x-[1])/[2])**2)");
      gaus_fnc->SetParLimits(0, 0., 100000);
      gaus_fnc->SetParLimits(1, 0., 100000);
      gaus_fnc->SetParLimits(2, 0., 100000);
      gaus_fnc->SetParameter(0, _TH1D[current_target_name]->GetMaximum());
      gaus_fnc->SetParameter(1, _TH1D[current_target_name]->GetMean());
      gaus_fnc->SetParameter(2, _TH1D[current_target_name]->GetRMS());
      _TH1D[current_target_name]->Fit(gaus_fnc, "IMESQ");
      if (target_fit->GetParError(0) / target_fit->GetParameter(0) > 0.1)
        _TH1D[current_target_name]->Fit(gaus_fnc, "IMESQ", "", _TH1D[current_target_name]->GetMean() - _TH1D[current_target_name]->GetRMS(), _TH1D[current_target_name]->GetMean() + _TH1D[current_target_name]->GetRMS());
      auto current_point = _TGraphErrors[current_graph]->GetN();
      _TGraphErrors[current_graph]->SetPoint(current_point, yBin_value, gaus_fnc->GetParameter(2));
      _TGraphErrors[current_graph]->SetPointError(current_point, 0., gaus_fnc->GetParError(2));
    }
    auto current_point = _TGraphErrors[current_TH3_target + "_Full_from_fit"]->GetN();
    _TGraphErrors[current_graph]->Fit(target_fit, "", "", min, max);
    _TGraphErrors[current_TH3_target + "_Full_from_fit"]->SetPoint(current_point, xBin_value, target_fit->GetParameter(0));
    _TGraphErrors[current_TH3_target + "_Full_from_fit"]->SetPointError(current_point, 0., target_fit->GetParError(0));
  }
}

void analyse_raw_data_for_radius_cut_calibration(std::string run_tag = "20231010-084623", std::string output_filename = "")
{
  //  Output
  TObjects output_struct;
  TObjects utility_struct;

  //  Input/Output files
  std::string input_filename = data_dir + "/" + run_tag + "/pre_process_recodata.root";
  if (output_filename.empty())
    output_filename = data_dir + "/" + run_tag + "/raw_prod_analysis.root";

  //  Recover raw data histogram
  TFile *input_file = new TFile(input_filename.c_str());
  utility_struct._TH3F["hDeltaRadius"] = (TH3F *)(input_file->Get("hDeltaRadius"));
  utility_struct._TH3F["hFitRadius"] = (TH3F *)(input_file->Get("hFitRadius"));

  //  Fit function
  TF1 *pol0_fit = new TF1("pol0_fit", "[0]");
  TF1 *sqrt_fit = new TF1("sqrt_fit", "TMath::Sqrt([1]*[1]+[0]*[0]/x)");
  sqrt_fit->SetParLimits(0, 0, 1000);

  //  Loop over the cut slices
  // AnalyseResolution("hDeltaRadius", utility_struct._TH3F, output_struct._TGraphErrors, output_struct._TH1D, pol0_fit, 10., 15.);
  AnalyseResolution("hFitRadius", utility_struct._TH3F, output_struct._TGraphErrors, output_struct._TH1D, sqrt_fit, 3.5, 15.5);

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(output_struct);
  output_file->Close();
}