#include "../lib/Utility.h"
#include "../lib/Database.h"
#include "../lib/mapping.h"

void photon_yield_analysis(std::string run_tag = "20231010-084623", std::string output_filename = "")
{
  //  Input/Output files
  std::string input_filename = data_dir + "/" + run_tag + "/pre_process_recodata.root";
  if (output_filename.empty())
    output_filename = data_dir + "/" + run_tag + "/photon_yield.root";
  TFile *input_file = new TFile(input_filename.c_str());

  //  Output & Utility
  TObjects output_struct;
  TObjects utility_struct;
  output_struct._TGraphErrors["Sig_gaus_component"] = new TGraphErrors();
  output_struct._TGraphErrors["Sig_pois_component"] = new TGraphErrors();
  output_struct._TGraphErrors["Bkg_poi1_component"] = new TGraphErrors();
  output_struct._TGraphErrors["Bkg_poi2_component"] = new TGraphErrors();

  // Functions used
  TF1 *utl_poisson = new TF1("utl_poisson", "[0]*TMath::Poisson(x,[1])+[2]*TMath::Poisson(x,[3])", 0, 100);
  utl_poisson->SetParLimits(0, 0., 100000);
  utl_poisson->SetParLimits(1, 0., 100000);
  utl_poisson->SetParLimits(2, 0., 100000);
  utl_poisson->SetParLimits(3, 0., 100000);
  TF1 *utl_poisson_and_gaus = new TF1("utl_poisson_and_gaus", "[0]*TMath::Poisson(x,[1])+[2]*TMath::Poisson(x,[3])")://[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*TMath::Poisson(x,[4])", 0, 100);
  utl_poisson_and_gaus->SetParLimits(0, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(1, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(2, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(3, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(4, 0., 100000);
  TF1 *utl_gaus = new TF1("utl_gaus", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0, 100);
  utl_gaus->SetParLimits(0, 0., 100000);
  utl_gaus->SetParLimits(1, 0., 100000);
  utl_gaus->SetParLimits(2, 0., 100000);
  utl_gaus->SetLineColor(kBlue);
  utl_gaus->SetLineStyle(kDashed);

  // Recover target signal photons
  output_struct._TH2F["hAveragePhotonEv"] = (TH2F *)(input_file->Get("hAveragePhotonEv"));
  output_struct._TH2F["hAveragePhotonEvDisc"] = (TH2F *)(input_file->Get("hAveragePhotonEvDisc"));

  for (auto iBin = 1; iBin <= output_struct._TH2F["hAveragePhotonEv"]->GetNbinsX(); iBin++)
  {
    auto current_selected_slice = output_struct._TH2F["hAveragePhotonEv"]->ProjectionY("tmp", iBin, iBin);
    auto current_discarded_slice = output_struct._TH2F["hAveragePhotonEvDisc"]->ProjectionY("tmpDisc", iBin, iBin);

    // Prepare Fit
    utl_poisson_and_gaus->SetParameters(current_selected_slice->GetEntries() * 0.5, current_selected_slice->GetMean(), current_selected_slice->GetRMS(), current_selected_slice->GetEntries() * 0.5, current_selected_slice->GetMean());
    utl_poisson->SetParameters(current_discarded_slice->GetEntries(), current_discarded_slice->GetMean(), 1000, 20);
    current_selected_slice->Fit(utl_poisson_and_gaus,"IMRESQ");
    current_discarded_slice->Fit(utl_poisson,"IMRESQ");
    current_selected_slice->Fit(utl_poisson_and_gaus,"IMRESQ");
    current_discarded_slice->Fit(utl_poisson,"IMRESQ");
    current_selected_slice->Fit(utl_poisson_and_gaus,"IMRESQ");
    current_discarded_slice->Fit(utl_poisson,"IMRESQ");
    auto current_point = output_struct._TGraphErrors["Sig_gaus_component"]->GetN();
    auto current_xvalue = output_struct._TH2F["hAveragePhotonEv"]->GetXaxis()->GetBinCenter(iBin);
    output_struct._TGraphErrors["Sig_gaus_component"]->SetPoint(current_point, current_xvalue, utl_poisson_and_gaus->GetParameter(1));
    output_struct._TGraphErrors["Sig_gaus_component"]->SetPointError(current_point, 0, utl_poisson_and_gaus->GetParameter(2));
    output_struct._TGraphErrors["Sig_pois_component"]->SetPoint(current_point, current_xvalue, utl_poisson_and_gaus->GetParameter(4));
    output_struct._TGraphErrors["Sig_pois_component"]->SetPointError(current_point, 0, utl_poisson_and_gaus->GetParError(4));
    output_struct._TGraphErrors["Bkg_poi1_component"]->SetPoint(current_point, current_xvalue, utl_poisson->GetParameter(1));
    output_struct._TGraphErrors["Bkg_poi1_component"]->SetPointError(current_point, 0, utl_poisson->GetParError(1));
    output_struct._TGraphErrors["Bkg_poi2_component"]->SetPoint(current_point, current_xvalue, utl_poisson->GetParameter(3));
    output_struct._TGraphErrors["Bkg_poi2_component"]->SetPointError(current_point, 0, utl_poisson->GetParError(3));
  }

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(output_struct);
  output_file->Close();
}