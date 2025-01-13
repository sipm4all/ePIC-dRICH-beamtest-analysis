#include "../lib/Utility.h"
#include "../lib/Database.h"
#include "../lib/mapping.h"
#include "./produce_raw_data_for_radius_cut_calibration.C"
#include "./analyse_raw_data_for_radius_cut_calibration.C"

std::array<float, 2> get_radius_sigma(std::string run_tag = "20231010-084623", float requested_sigma_cut = 3.)
{
  //  Result
  std::array<float, 2> result = {-1, -1};

  //  Open target file
  std::string input_filename = data_dir + "/" + run_tag + "/raw_prod_analysis.root";
  TFile *target_file = new TFile(input_filename.c_str());

  //  Get the radius_sigma graph
  auto target_graph_err = (TGraphErrors *)(target_file->Get("hFitRadius_Full_from_fit"));

  //  Look for requested sigma
  for (auto iPnt = 0; iPnt < target_graph_err->GetN(); iPnt++)
  {
    float current_x = target_graph_err->GetPointX(iPnt);
    float current_y = target_graph_err->GetPointY(iPnt);
    float current_y_err = target_graph_err->GetErrorY(iPnt);
    if (fabs(current_x - requested_sigma_cut) < 0.1)
    {
      result = {current_y, current_y_err};
      break;
    }
  }

  //  Close target file
  target_file->Close();

  //  Return
  return result;
}

std::array<float, 2> get_from_histogram(std::string run_tag = "20231010-084623", float requested_sigma_cut = 3., std::string target_file_name = "raw_prod", std::string target_histogram_name = "hFitRadiusVal", TF1 *target_function = nullptr)
{
  //  Result
  std::array<float, 2> result = {-1, -1};

  //  Open target file
  std::string input_filename = data_dir + "/" + run_tag + "/" + target_file_name + ".root";
  TFile *target_file = new TFile(input_filename.c_str());

  //  Get the radius_sigma graph
  auto target_histo = (TH2F *)(target_file->Get(target_histogram_name.c_str()));

  //  Look for requested sigma
  auto raget_bin = target_histo->GetXaxis()->FindBin(requested_sigma_cut);
  auto target_profile = target_histo->ProjectionY("tmp", raget_bin, raget_bin);
  result = {static_cast<float>(target_profile->GetMean()), static_cast<float>(target_profile->GetMeanError())};
  if (target_function)
  {
    // TCanvas *c1 = new TCanvas();
    target_function->SetParameter(0, target_profile->GetEntries());
    target_function->SetParameter(1, target_profile->GetMean());
    for (int i = 0; i < 25; i++)
      target_profile->Fit(target_function, "Q");
    // target_profile->DrawCopy();
    // target_function->DrawCopy("SAME");
  }

  //  Close target file
  target_file->Close();

  //  Return
  return result;
}

std::array<float, 2> get_radius(std::string run_tag = "20231010-084623", float requested_sigma_cut = 3., TF1 *target_function = nullptr)
{
  return get_from_histogram(run_tag, requested_sigma_cut, "raw_prod", "hFitRadiusVal", target_function);
}

std::array<float, 2> get_accepted_photons(std::string run_tag = "20231010-084623", float requested_sigma_cut = 3., TF1 *target_function = nullptr)
{
  return get_from_histogram(run_tag, requested_sigma_cut, "raw_prod", "hAveragePhotonEv", target_function);
}

std::array<float, 2> get_rejected_photons(std::string run_tag = "20231010-084623", float requested_sigma_cut = 3., TF1 *target_function = nullptr)
{
  return get_from_histogram(run_tag, requested_sigma_cut, "raw_prod", "hAveragePhotonEvDisc", target_function);
}

void test_macro_for_poster(std::string run_tag = "20231010-084623", std::string output_filename = "")
{
  //  Output
  TObjects output_struct;
  TObjects utility_struct;
  output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_radius_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_radius_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"] = new TGraphErrors();
  output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"] = new TGraphErrors();

  //  Input/Output files
  std::string input_filename = data_dir + "/" + run_tag + "/raw_prod_analysis.root";
  if (output_filename.empty())
    output_filename = data_dir + "/" + run_tag + "/raw_prod_analysis_2.root";

  // produce_raw_data_for_radius_cut_calibration("20231010-225620");

  TF1 *utl_poisson = new TF1("utl_poisson", "[0]*TMath::Poisson(x,[1])");
  utl_poisson->SetParLimits(0, 0., 100000);
  utl_poisson->SetParLimits(1, 0., 100000);

  TF1 *utl_poisson_and_gaus = new TF1("utl_poisson_and_gaus", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*TMath::Poisson(x,[1])");
  utl_poisson_and_gaus->SetParLimits(0, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(1, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(2, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(3, 0., 100000);

  TF1 *utl_gaus = new TF1("utl_gaus", "[0]*exp(-0.5*((x-[1])/[2])**2)");
  utl_poisson_and_gaus->SetParLimits(0, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(1, 0., 100000);
  utl_poisson_and_gaus->SetParLimits(2, 0., 100000);

  auto iPnt = -1;
  for (auto current_run : database::run_lists["mirror_scan_bkw"])
  {
    iPnt++;

    //  Re-run preliminary analysis
    gROOT->SetBatch();
    // produce_raw_data_for_radius_cut_calibration(current_run);
    // analyse_raw_data_for_radius_cut_calibration(current_run);
    gROOT->SetBatch(false);

    //  Data analysis information
    //  --- SPR Radius
    auto SPR_radius = get_radius_sigma(current_run);
    auto sigma_radius_res_val = SPR_radius[0];
    auto sigma_radius_res_err = SPR_radius[1];
    //  --- Radius & Radius sigma
    utl_gaus->SetParameters(1., 1., 1.);
    get_radius(current_run, 3, utl_gaus);
    auto radius_val = utl_gaus->GetParameter(1);
    auto radius_err = utl_gaus->GetParError(1);
    auto sigma_radius_val = utl_gaus->GetParameter(2);
    auto sigma_radius_err = utl_gaus->GetParError(2);
    //  --- Accepted photons
    utl_poisson_and_gaus->SetParameters(1., 1., 1., 1.);
    get_accepted_photons(current_run, 3, utl_poisson_and_gaus);
    auto acc_ph_val = utl_poisson_and_gaus->GetParameter(1);
    auto acc_ph_err = utl_poisson_and_gaus->GetParError(1);
    //  --- Rejected photons
    utl_poisson->SetParameters(1., 1.);
    get_rejected_photons(current_run, 3, utl_poisson);
    auto rej_ph_val = utl_poisson->GetParameter(1);
    auto rej_ph_err = utl_poisson->GetParError(1);
    //  --- Derived quantities
    auto sigma_rel_radius_res_val = 100 * sigma_radius_res_val / radius_val;
    auto sigma_rel_radius_res_err = sigma_rel_radius_res_val * TMath::Sqrt((sigma_radius_res_err / sigma_radius_res_val) * (sigma_radius_res_err / sigma_radius_res_val) + (radius_err / radius_val) * (radius_err / radius_val));
    auto sigma_rel_radius_val = 100 * sigma_radius_val / radius_val;
    auto sigma_rel_radius_err = TMath::Sqrt((sigma_radius_err / sigma_radius_val) * (sigma_radius_err / sigma_radius_val) + (radius_err / radius_val) * (radius_err / radius_val));

    //  Data analysis information
    auto aerogel_mirror_z_position = database::get_aerogel_mirror_z(current_run).second;

    //  Make graphs
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, sigma_radius_res_val);
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetPointError(iPnt, 0, sigma_radius_res_err);
    output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, radius_val);
    output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->SetPointError(iPnt, 0, radius_err);
    output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, sigma_rel_radius_res_val);
    output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->SetPointError(iPnt, 0, sigma_rel_radius_res_err);
    output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, acc_ph_val);
    output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->SetPointError(iPnt, 0, acc_ph_err);
    output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, rej_ph_val);
    output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->SetPointError(iPnt, 0, rej_ph_err);
    output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, acc_ph_val / (acc_ph_val + rej_ph_val));
    output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, rej_ph_val / (acc_ph_val + rej_ph_val));
  }

  iPnt = -1;
  for (auto current_run : database::run_lists["mirror_scan_fwd"])
  {
    iPnt++;
    //  Re-run preliminary analysis
    gROOT->SetBatch();
    //   produce_raw_data_for_radius_cut_calibration(current_run);
    // analyse_raw_data_for_radius_cut_calibration(current_run);
    gROOT->SetBatch(false);

    //  Data analysis information
    //  --- SPR Radius
    auto SPR_radius = get_radius_sigma(current_run);
    auto sigma_radius_res_val = SPR_radius[0];
    auto sigma_radius_res_err = SPR_radius[1];
    //  --- Radius & Radius sigma
    utl_gaus->SetParameters(1., 1., 1.);
    get_radius(current_run, 3, utl_gaus);
    auto radius_val = utl_gaus->GetParameter(1);
    auto radius_err = utl_gaus->GetParError(1);
    auto sigma_radius_val = utl_gaus->GetParameter(2);
    auto sigma_radius_err = utl_gaus->GetParError(2);
    //  --- Accepted photons
    utl_poisson_and_gaus->SetParameters(1., 1., 1., 1.);
    get_accepted_photons(current_run, 3, utl_poisson_and_gaus);
    auto acc_ph_val = utl_poisson_and_gaus->GetParameter(1);
    auto acc_ph_err = utl_poisson_and_gaus->GetParError(1);
    //  --- Rejected photons
    utl_poisson->SetParameters(1., 1.);
    get_rejected_photons(current_run, 3, utl_poisson);
    auto rej_ph_val = utl_poisson->GetParameter(1);
    auto rej_ph_err = utl_poisson->GetParError(1);
    //  --- Derived quantities
    auto sigma_rel_radius_res_val = 100 * sigma_radius_res_val / radius_val;
    auto sigma_rel_radius_res_err = sigma_rel_radius_res_val * TMath::Sqrt((sigma_radius_res_err / sigma_radius_res_val) * (sigma_radius_res_err / sigma_radius_res_val) + (radius_err / radius_val) * (radius_err / radius_val));
    auto sigma_rel_radius_val = 100 * sigma_radius_val / radius_val;
    auto sigma_rel_radius_err = sigma_rel_radius_val * TMath::Sqrt((sigma_radius_err / sigma_radius_val) * (sigma_radius_err / sigma_radius_val) + (radius_err / radius_val) * (radius_err / radius_val));

    //  Data analysis information
    auto aerogel_mirror_z_position = database::get_aerogel_mirror_z(current_run).second;

    //  Make graphs
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, sigma_radius_res_val);
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetPointError(iPnt, 0, sigma_radius_res_err);
    output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, radius_val);
    output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->SetPointError(iPnt, 0, radius_err);
    output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, sigma_rel_radius_res_val);
    output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->SetPointError(iPnt, 0, sigma_rel_radius_res_err);
    output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, acc_ph_val);
    output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->SetPointError(iPnt, 0, acc_ph_err);
    output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, rej_ph_val);
    output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->SetPointError(iPnt, 0, rej_ph_err);
    output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, acc_ph_val / (acc_ph_val + rej_ph_val));
    output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"]->SetPoint(iPnt, aerogel_mirror_z_position, rej_ph_val / (acc_ph_val + rej_ph_val));
  }

  TH1F *frame = new TH1F("frame", "", 1, 0, 145);
  frame->SetBinContent(1, -999);
  gStyle->SetOptStat(0);

  TCanvas *radius_vs_mirror_position = plot_standard_canvas("radius_vs_mirror_position");
  TLegend *leg_radius_vs_mirror_position = new TLegend(0.15, 0.15, 0.50, 0.50);
  leg_radius_vs_mirror_position->SetFillStyle(0);
  leg_radius_vs_mirror_position->SetLineColorAlpha(0., 0.);
  leg_radius_vs_mirror_position->AddEntry(output_struct._TGraphErrors["bkw_radius_vs_mirror_position"], "1^{st} scan", "PE");
  leg_radius_vs_mirror_position->AddEntry(output_struct._TGraphErrors["fwd_radius_vs_mirror_position"], "2^{nd} scan", "PE");
  frame->GetYaxis()->SetRangeUser(68, 74);
  frame->SetTitle(";Aerogel mirror position (mm);Average radius measured (mm)");
  frame->DrawCopy();
  output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->SetMarkerSize(3);
  output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->SetMarkerColor(kRed + 1);
  output_struct._TGraphErrors["bkw_radius_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->SetMarkerSize(3);
  output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->SetMarkerColor(kAzure - 1);
  output_struct._TGraphErrors["fwd_radius_vs_mirror_position"]->Draw("SAME PE");
  leg_radius_vs_mirror_position->Draw("SAME");
  radius_vs_mirror_position->SaveAs("radius_vs_mirror_position.pdf");

  TCanvas *photon_vs_mirror_position = plot_standard_canvas("photon_vs_mirror_position");
  TLegend *leg_photon_vs_mirror_position = new TLegend(0.15, 0.15, 0.9, 0.35);
  leg_photon_vs_mirror_position->SetFillStyle(0);
  leg_photon_vs_mirror_position->SetNColumns(2);
  leg_photon_vs_mirror_position->SetLineColorAlpha(0., 0.);
  leg_photon_vs_mirror_position->AddEntry(output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"], "1^{st} scan accepted photons", "PE");
  leg_photon_vs_mirror_position->AddEntry(output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"], "1^{st} scan rejected photons", "PE");
  leg_photon_vs_mirror_position->AddEntry(output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"], "2^{nd} scan accepted photons", "PE");
  leg_photon_vs_mirror_position->AddEntry(output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"], "2^{nd} scan rejected photons", "PE");
  frame->GetYaxis()->SetRangeUser(1.5, 9.);
  frame->SetTitle(";Aerogel mirror position (mm);Average number of photons considered");
  frame->DrawCopy();
  output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->SetMarkerStyle(20);
  output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->SetMarkerColor(kGreen + 1);
  output_struct._TGraphErrors["bkw_acc_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->SetMarkerStyle(20);
  output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->SetMarkerColor(kRed + 1);
  output_struct._TGraphErrors["bkw_rej_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->SetMarkerColor(kGreen + 3);
  output_struct._TGraphErrors["fwd_acc_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->SetMarkerColor(kRed + 3);
  output_struct._TGraphErrors["fwd_rej_photons_vs_mirror_position"]->Draw("SAME PE");
  leg_photon_vs_mirror_position->Draw("SAME");
  photon_vs_mirror_position->SaveAs("photon_vs_mirror_position.pdf");

  TCanvas *photon_fraction_vs_mirror_position = plot_standard_canvas("photon_fraction_vs_mirror_position");
  frame->GetYaxis()->SetRangeUser(0.1, 0.75);
  frame->SetTitle(";Aerogel mirror position (mm);Fraction of all photons (\%)");
  frame->DrawCopy();
  output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"]->SetMarkerStyle(20);
  output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"]->SetMarkerColor(kGreen + 1);
  output_struct._TGraphErrors["bkw_rel_acc_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"]->SetMarkerStyle(20);
  output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"]->SetMarkerColor(kRed + 1);
  output_struct._TGraphErrors["bkw_rel_rej_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"]->SetMarkerColor(kGreen + 3);
  output_struct._TGraphErrors["fwd_rel_acc_photons_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"]->SetMarkerColor(kRed + 3);
  output_struct._TGraphErrors["fwd_rel_rej_photons_vs_mirror_position"]->Draw("SAME PE");
  leg_photon_vs_mirror_position->Draw("SAME");
  photon_fraction_vs_mirror_position->SaveAs("photon_fraction_vs_mirror_position.pdf");

  TCanvas *radius_sigma_vs_mirror_position = plot_standard_canvas("radius_sigma_vs_mirror_position");
  TLegend *leg_radius_sigma_vs_mirror_position = new TLegend(0.15, 0.15, 0.90, 0.35);
  leg_radius_sigma_vs_mirror_position->SetFillStyle(0);
  leg_radius_sigma_vs_mirror_position->SetNColumns(2);
  leg_radius_sigma_vs_mirror_position->SetLineColorAlpha(0., 0.);
  leg_radius_sigma_vs_mirror_position->AddEntry(output_struct._TGraphErrors["bkw_radius_vs_mirror_position"], "1^{st} scan", "PE");
  leg_radius_sigma_vs_mirror_position->AddEntry(output_struct._TGraphErrors["fwd_radius_vs_mirror_position"], "2^{nd} scan", "PE");
  frame->GetYaxis()->SetRangeUser(1.0, 2.5);
  frame->SetTitle(";Aerogel mirror position (mm);Single-photon radius resolution (mm)");
  frame->DrawCopy();
  output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerColor(kRed + 1);
  output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerColor(kAzure - 1);
  output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->Draw("SAME PE");
  leg_radius_sigma_vs_mirror_position->Draw("SAME");
  radius_sigma_vs_mirror_position->SaveAs("radius_sigma_vs_mirror_position.pdf");

  TCanvas *radius_rel_sigma_vs_mirror_position = plot_standard_canvas("radius_rel_sigma_vs_mirror_position");
  frame->GetYaxis()->SetRangeUser(1.5, 3.5);
  frame->SetTitle(";Aerogel mirror position (mm);Single-photon radius resolution (\%)");
  frame->DrawCopy();
  output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->SetMarkerColor(kRed + 1);
  output_struct._TGraphErrors["bkw_radius_rel_res_vs_mirror_position"]->Draw("SAME PE");
  output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->SetMarkerStyle(21);
  output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->SetMarkerSize(2);
  output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->SetMarkerColor(kAzure - 1);
  output_struct._TGraphErrors["fwd_radius_rel_res_vs_mirror_position"]->Draw("SAME PE");
  leg_radius_sigma_vs_mirror_position->Draw("SAME");
  radius_rel_sigma_vs_mirror_position->SaveAs("radius_rel_sigma_vs_mirror_position.pdf");

  /*

      TCanvas *radius_sigma_vs_mirror_position = plot_standard_canvas("radius_sigma_vs_mirror_position");
    frame->GetYaxis()->SetRangeUser(1.5, 2.5);
    frame->SetTitle(";Aerogel mirror position (mm);Average radius sigma (mm)");
    frame->DrawCopy();
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerStyle(20);
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerSize(2);
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->SetMarkerColor(kBlack);
    output_struct._TGraphErrors["bkw_radius_res_vs_mirror_position"]->Draw("SAME PE");
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerStyle(21);
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerSize(2);
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->SetMarkerColor(kBlack);
    output_struct._TGraphErrors["fwd_radius_res_vs_mirror_position"]->Draw("SAME PE");

    */

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(output_struct);
  output_file->Close();
}