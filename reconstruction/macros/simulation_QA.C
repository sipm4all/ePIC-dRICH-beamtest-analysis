#include "./produce_raw_data_for_radius_cut_calibration.C"

void recover_all_histograms(std::map<std::string, TObjects> &target, std::string run_tag, std::vector<std::string> secondary_check)
{
  // build list of target input files
  std::vector<std::string> all_filenames;
  all_filenames.push_back(run_tag);
  all_filenames.push_back("MC" + run_tag);
  for (auto current_secondary : secondary_check)
  {
    all_filenames.push_back("MC" + run_tag + "-" + current_secondary);
  }

  // Loop over files
  TFile *input_file;
  for (auto current_run_tag : all_filenames)
  {
    std::string input_filename = data_dir + "/" + current_run_tag + "/raw_prod.root";
    if (!current_run_tag.rfind("MC"))
      input_filename = sim_dir + "/" + current_run_tag + "/raw_prod.root";
    input_file = new TFile(input_filename.c_str());
    if (input_file->IsZombie())
      produce_raw_data_for_radius_cut_calibration(current_run_tag);
    load_production_step(target[current_run_tag], input_file, current_run_tag);
  }
}

void recover_all_histograms(std::map<std::string, TObjects> &target, std::vector<std::string> list_of_run_tags)
{
  // Loop over files
  TFile *input_file;
  for (auto current_run_tag : list_of_run_tags)
  {
    std::string input_filename = data_dir + "/" + current_run_tag + "/raw_prod.root";
    if (!current_run_tag.rfind("MC"))
      input_filename = sim_dir + "/" + current_run_tag + "/raw_prod.root";
    input_file = new TFile(input_filename.c_str());
    if (input_file->IsZombie())
      produce_raw_data_for_radius_cut_calibration(current_run_tag);
    load_production_step(target[current_run_tag], input_file, current_run_tag);
  }
}

// TODO:
std::vector<int> kColorMap = {kRed, kBlue + 1, kGreen + 1, kYellow + 1};

TCanvas *
compare_coordinate(std::map<std::string, TObjects> inputContainer, std::string coordinate_tag, float cut_value, std::array<float, 4> draw_limits)
{
  TCanvas *cCompareCoordinates = new TCanvas("cCompareCoordinates", "cCompareCoordinates", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  TLegend *lCompareCoordinates = new TLegend(0.17, 0.93, 0.57, 0.73);
  lCompareCoordinates->SetFillColorAlpha(0., 0.);
  lCompareCoordinates->SetLineColorAlpha(0., 0.);
  TLatex *ltxCompareCoordinates = new TLatex();
  auto iTer = 0;
  for (auto [current_run_tag, current_data_struct] : inputContainer)
  {
    iTer++;
    auto target_bin = current_data_struct._TH2F[coordinate_tag]->GetXaxis()->FindBin(cut_value);
    auto hCurrentCut = current_data_struct._TH2F[coordinate_tag]->ProjectionY(Form("tmp_%s_%s", coordinate_tag.c_str(), current_run_tag.c_str()), target_bin, target_bin);
    hCurrentCut->GetXaxis()->SetRangeUser(draw_limits[0], draw_limits[1]);
    hCurrentCut->GetYaxis()->SetRangeUser(draw_limits[2], draw_limits[3]);
    hCurrentCut->SetLineColor(kColorMap[iTer - 1]);
    hCurrentCut->Draw("SAME");
    ltxCompareCoordinates->SetTextColor(kColorMap[iTer - 1]);
    ltxCompareCoordinates->DrawLatexNDC(0.18, 0.70 - iTer * 0.05, Form("#mu: %.2f #sigma: %.2f", hCurrentCut->GetMean(), hCurrentCut->GetRMS()));
    lCompareCoordinates->AddEntry(hCurrentCut, current_run_tag.c_str());
  }
  lCompareCoordinates->Draw("SAME");
  return cCompareCoordinates;
}

TGraphErrors *
compare_gamma_quantity(std::map<std::string, TObjects> inputContainer, std::string quantity_tag, float cut_value, std::array<float, 4> draw_limits)
{
  TCanvas *cCompareCoordinates = new TCanvas("cCompareCoordinates", "cCompareCoordinates", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  TLegend *lCompareCoordinates = new TLegend(0.17, 0.93, 0.57, 0.73);
  lCompareCoordinates->SetFillColorAlpha(0., 0.);
  lCompareCoordinates->SetLineColorAlpha(0., 0.);
  TLatex *ltxCompareCoordinates = new TLatex();
  TGraphErrors *gTest = new TGraphErrors();
  auto iTer = 0;
  for (auto [current_run_tag, current_data_struct] : inputContainer)
  {
    TCanvas *tmp = new TCanvas();
    iTer++;
    auto target_bin = current_data_struct._TH3F[quantity_tag]->GetXaxis()->FindBin(cut_value);
    current_data_struct._TH3F[quantity_tag]->GetXaxis()->SetRange(target_bin, target_bin);
    auto hCurrentCut2D = (TH2F *)(current_data_struct._TH3F[quantity_tag]->Project3D("zy"));
    auto hCurrentCut = hCurrentCut2D->ProjectionX(Form("tmp_%s_%s", quantity_tag.c_str(), current_run_tag.c_str()), -1, -1, "");
    for (auto iBin = 0; iBin < hCurrentCut2D->GetNbinsX(); iBin++)
    {
      auto current_content = 0.;
      auto current_error = 0.;
      auto hCurrentCutTmp = hCurrentCut2D->ProjectionY(Form("tmp_%s_%s_%i", quantity_tag.c_str(), current_run_tag.c_str(), iBin), iBin, iBin, "");
      if (hCurrentCutTmp->Integral() > 100)
      {
        hCurrentCutTmp->Fit("gaus");
        auto gaus_fit = hCurrentCutTmp->GetFunction("gaus");
        current_content = gaus_fit->GetParameter(2);
        current_error = gaus_fit->GetParError(2);
      }
      hCurrentCut->SetBinContent(iBin, current_content);
      hCurrentCut->SetBinError(iBin, current_error);
    }
    hCurrentCut->SetTitle("");
    hCurrentCut->GetXaxis()->SetRangeUser(draw_limits[0], draw_limits[1]);
    hCurrentCut->GetYaxis()->SetRangeUser(draw_limits[2], draw_limits[3]);
    hCurrentCut->SetLineColor(kColorMap[iTer - 1]);

    // Fit
    TF1 *fit = new TF1("fit", "[0]/TMath::Sqrt(x)");
    hCurrentCut->Fit(fit);
    delete tmp;
    cCompareCoordinates->cd();
    hCurrentCut->Draw("SAME");
    ltxCompareCoordinates->SetTextColor(kColorMap[iTer - 1]);
    ltxCompareCoordinates->DrawLatexNDC(0.57, 0.93 - iTer * 0.05, Form("[0]: %.3f #pm %.3f", fit->GetParameter(0), fit->GetParError(0)));

    auto npoint = gTest->GetN();
    gTest->SetPoint(npoint, std::stoi(current_run_tag.substr(current_run_tag.find("Nbkg") + 5)), fit->GetParameter(0));
    gTest->SetPointError(npoint, 0, fit->GetParError(0));
    lCompareCoordinates->AddEntry(hCurrentCut, current_run_tag.c_str());
  }
  lCompareCoordinates->Draw("SAME");
  TCanvas *all = new TCanvas("all", "all", 1000, 1000);
  gTest->Draw("ALPE");
  return gTest;
}

TCanvas *compare_bkg_evaluation(std::map<std::string, TObjects> inputContainer, float cut_value, std::array<float, 4> draw_limits, TObjects &outputContainer)
{
  std::string quantity_tag = "hAveragePhotonEvDisc";

  TCanvas *cCompare_bkg_evaluation = new TCanvas("cCompare_bkg_evaluation", "cCompare_bkg_evaluation", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->DrawFrame(draw_limits[0], draw_limits[2], draw_limits[1], draw_limits[3], ";Average bkg hits;Background level");

  TLegend *lCompare_bkg_evaluation = new TLegend(0.17, 0.93, 0.57, 0.73);
  lCompare_bkg_evaluation->SetFillColorAlpha(0., 0.);
  lCompare_bkg_evaluation->SetLineColorAlpha(0., 0.);

  TLatex *ltxCompare_bkg_evaluation = new TLatex();

  outputContainer._TGraphErrors["gCompare_bkg_evaluation"] = new TGraphErrors();
  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetMarkerStyle(20);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetMarkerColor(kRed);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"] = new TGraphErrors();
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetMarkerStyle(20);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetMarkerColor(kBlue);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetLineColor(kBlue);
  auto iTer = 0;
  for (auto [current_run_tag, current_data_struct] : inputContainer)
  {
    TCanvas *tmp_dump = new TCanvas();

    iTer++;
    auto target_bin = current_data_struct._TH2F[quantity_tag]->GetXaxis()->FindBin(cut_value);
    auto hCurrentCut = (TH1D *)(current_data_struct._TH2F[quantity_tag]->ProjectionY(Form("tmp_%s_%s_%i", quantity_tag.c_str(), current_run_tag.c_str(), target_bin), target_bin, target_bin));

    fPoissonian->SetParameters(1., 1., 1.);
    fPoissonian->FixParameter(2, 1.);
    hCurrentCut->Fit(fPoissonian, "");
    hCurrentCut->SetTitle("");
    hCurrentCut->SetLineColor(kColorMap[iTer - 1]);

    delete tmp_dump;

    if (!current_run_tag.rfind("MC"))
    {
      auto current_point = outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetN();
      outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetPoint(current_point, std::stoi(current_run_tag.substr(current_run_tag.find("Nbkg") + 5)), fPoissonian->GetParameter(1));
      outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetPointError(current_point, 0, fPoissonian->GetParError(1));
    }
    else
    {
      auto current_point = outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->GetN();
      outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetPoint(current_point, 1, fPoissonian->GetParameter(1));
      outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetPointError(current_point, 100, fPoissonian->GetParError(1));
    }
  }

  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->Draw("SAME PE");
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->Draw("SAME PE");

  lCompare_bkg_evaluation->AddEntry(outputContainer._TGraphErrors["gCompare_bkg_evaluation"], "Sim");
  lCompare_bkg_evaluation->AddEntry(outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"], "Data");
  lCompare_bkg_evaluation->Draw("SAME");

  return cCompare_bkg_evaluation;
}

int compare_bkg_evaluation(std::map<std::string, TObjects> inputContainer, std::array<float, 5> cut_value, std::array<float, 4> draw_limits, TObjects &outputContainer, std::string check_result_image = "")
{
  std::string quantity_tag = "hPersistance";

  TCanvas *cCompare_bkg_evaluation = new TCanvas("cCompare_bkg_evaluation", "cCompare_bkg_evaluation", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->DrawFrame(draw_limits[0], draw_limits[2], draw_limits[1], draw_limits[3], ";Average bkg hits;Background level");

  TLegend *lCompare_bkg_evaluation = new TLegend(0.17, 0.93, 0.57, 0.73);
  lCompare_bkg_evaluation->SetFillColorAlpha(0., 0.);
  lCompare_bkg_evaluation->SetLineColorAlpha(0., 0.);

  TLatex *ltxCompare_bkg_evaluation = new TLatex();

  outputContainer._TGraphErrors["gCompare_bkg_evaluation"] = new TGraphErrors();
  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetMarkerStyle(20);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetMarkerColor(kRed);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"] = new TGraphErrors();
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetMarkerStyle(20);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetMarkerColor(kBlue);
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetLineColor(kBlue);
  auto iTer = 0;
  for (auto [current_run_tag, current_data_struct] : inputContainer)
  {
    TCanvas *tmp_dump = new TCanvas();

    iTer++;
    current_data_struct._TH2F[quantity_tag]->GetXaxis()->SetRangeUser(cut_value[1], cut_value[2]);
    current_data_struct._TH2F[quantity_tag]->GetYaxis()->SetRangeUser(cut_value[3], cut_value[4]);
    auto xLow_bin = current_data_struct._TH2F[quantity_tag]->GetXaxis()->FindBin(cut_value[1]);
    auto xHig_bin = current_data_struct._TH2F[quantity_tag]->GetXaxis()->FindBin(cut_value[2]);
    auto yLow_bin = current_data_struct._TH2F[quantity_tag]->GetYaxis()->FindBin(cut_value[3]);
    auto yHig_bin = current_data_struct._TH2F[quantity_tag]->GetYaxis()->FindBin(cut_value[4]);
    auto current_histo_name = Form("hPersistance_%s_%s", quantity_tag.c_str(), current_run_tag.c_str());
    outputContainer._TH2F[current_histo_name] = (TH2F *)(current_data_struct._TH2F[quantity_tag]->Clone(Form("hPersistance_%s_%s", quantity_tag.c_str(), current_run_tag.c_str())));
    outputContainer._TH1D[current_histo_name] = (current_data_struct._TH2F[quantity_tag]->ProjectionX(Form("hPersistance_%s_%s_2", quantity_tag.c_str(), current_run_tag.c_str()), yLow_bin, yHig_bin));

    std::string inFile_recodata = data_dir + "/" + current_run_tag + "/recodata.root";
    if (!current_run_tag.rfind("MC"))
      inFile_recodata = sim_dir + "/" + current_run_tag + "/recodata.root";
    recodata reco_data;
    auto reco_tree = load_data(inFile_recodata, reco_data);
    auto current_entries = reco_tree->GetEntries();
    outputContainer._TH2F[current_histo_name]->Scale(1. / current_entries);
    outputContainer._TH1D[current_histo_name]->Scale(1. / current_entries);
    auto current_integral_error = 0.;
    auto current_integral = outputContainer._TH2F[current_histo_name]->IntegralAndError(xLow_bin, xHig_bin, yLow_bin, yHig_bin, current_integral_error);

    TF1 *fPol0 = new TF1("fPol0", "[0]");
    outputContainer._TH1D[current_histo_name]->Fit(fPol0, "QI");
    current_integral = fPol0->GetParameter(0);
    current_integral_error = fPol0->GetParError(0);

    delete tmp_dump;

    if (!current_run_tag.rfind("MC"))
    {
      auto current_point = outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetN();
      outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetPoint(current_point, std::stoi(current_run_tag.substr(current_run_tag.find("Nbkg") + 5)), current_integral);
      outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->SetPointError(current_point, 0, current_integral_error);
    }
    else
    {
      auto current_point = outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->GetN();
      outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetPoint(current_point, 1, current_integral);
      outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->SetPointError(current_point, 100, current_integral_error);
    }
  }

  outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->Draw("SAME PE");
  outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->Draw("SAME PE");

  lCompare_bkg_evaluation->AddEntry(outputContainer._TGraphErrors["gCompare_bkg_evaluation"], "Sim");
  lCompare_bkg_evaluation->AddEntry(outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"], "Data");
  lCompare_bkg_evaluation->Draw("SAME");

  if (!check_result_image.empty())
    cCompare_bkg_evaluation->SaveAs(check_result_image.c_str());

  auto result = -1;
  auto target = outputContainer._TGraphErrors["gCompare_bkg_evaluation_data"]->GetPointY(0);
  for (auto iPnt = 0; outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetN(); iPnt++)
  {
    if (outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetPointY(iPnt) < target)
      continue;
    result = outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetPointX(iPnt);
    if (fabs(outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetPointY(iPnt) - target) > fabs(outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetPointY(iPnt - 1) - target))
      result = outputContainer._TGraphErrors["gCompare_bkg_evaluation"]->GetPointX(iPnt - 1);
    break;
  }
  return result;
}

void simulation_QA_(std::string run_tag = "20231010-084623", std::vector<std::string> secondary_check = {"all"}, std::string output_filename = "out.root")
{
  //  Output
  std::map<std::string, TObjects> inputContainer;
  TObjects outputContainter;

  //  Recover raw data histogram
  recover_all_histograms(inputContainer, run_tag, secondary_check);

  //  Compare coordinates
  compare_coordinate(inputContainer, "hAverageX0", 5., {-10, 10, 0., 1400});
  compare_coordinate(inputContainer, "hAverageY0", 5., {-10, 10, 0., 1400});
  compare_coordinate(inputContainer, "hAverageR", 5., {60, 80, 0., 2000});
  compare_gamma_quantity(inputContainer, "hFitRadius", 5., {4, 30, 0, 3.5});

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(outputContainter);
  output_file->Close();
}

TCanvas *plot_squares(TObjects inputContainer, std::vector<std::array<float, 4>> perimetral_values, std::array<float, 4> draw_limits)
{
  TCanvas *persistance_map_with_squares = new TCanvas("persistance_map_with_squares", "persistance_map_with_squares", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->DrawFrame(draw_limits[0], draw_limits[2], draw_limits[1], draw_limits[3], ";X (mm); Y (mm)");
  inputContainer._TH2F["hPersistance"]->Draw("COLZ");
  /*
    for (auto current_values : perimetral_values)
    {
      auto x_low = current_values[0];
      auto x_hig = current_values[1];
      auto y_low = current_values[2];
      auto y_hig = current_values[3];
      TBox *current_box = new TBox();
    }
  */
  return persistance_map_with_squares;
}

std::array<std::array<float, 2>, 3> get_prediction(TF1 *fSin)
{
  std::array<std::array<float, 2>, 3> result;

  //  X_{0} prediction
  result[0][0] = -fSin->GetParameter(1) * TMath::Sin(fSin->GetParameter(2));
  result[0][1] = TMath::Sqrt(TMath::Sq(fSin->GetParError(1) * TMath::Sin(fSin->GetParameter(2))) + TMath::Sq(fSin->GetParError(2) * fSin->GetParameter(1) * TMath::Cos(fSin->GetParameter(2))));

  //  Y_{0} prediction
  result[1][0] = fSin->GetParameter(1) * TMath::Cos(fSin->GetParameter(2));
  result[1][1] = TMath::Sqrt(TMath::Sq(fSin->GetParError(1) * TMath::Cos(fSin->GetParameter(2))) + TMath::Sq(fSin->GetParError(2) * fSin->GetParameter(1) * TMath::Sin(fSin->GetParameter(2))));

  //  Radius prediction
  result[2][0] = fSin->GetParameter(0);
  result[2][1] = fSin->GetParError(0);

  cout << "Prediction > X0: " << result[0][0] << "+-" << result[0][1] << " Y0: " << result[1][0] << "+-" << result[1][1] << " R: " << result[2][0] << "+-" << result[2][1] << endl;
  return result;
}

float calculate_efficiency()
{
  return 1.;
}

std::array<std::array<float, 2>, 3> analyse_efficiency(TObjects hTarget_input)
{
  //  Create a local copy
  auto hTarget = (TH2F *)(hTarget_input._TH2F["hEfficiencyMap"]->Clone("hTarget"));

  //
  TCanvas *cCheckResults = new TCanvas("cCheckResults", "cCheckResults", 1000, 500);
  cCheckResults->Divide(2, 1);
  gStyle->SetOptStat(0);
  TCanvas *tmp = new TCanvas("tmp", "tmp", 1000, 500);

  //   Define functions
  TF2 *full_fit = new TF2("full_fit", "[5]+[3]*exp(-0.5*((y-([0]+[1]*TMath::Sin(x-[2])))/[4])**2)", -1000, 1000, -1000, 1000);
  full_fit->SetNpx(1.e7);
  full_fit->SetParLimits(0, 0., 1000.);
  full_fit->SetParLimits(1, 0., 1000.);
  full_fit->SetParLimits(2, -TMath::Pi(), TMath::Pi());
  full_fit->SetParLimits(3, 0., 100000.);
  full_fit->SetParLimits(4, 0., 1000.);
  full_fit->SetParLimits(5, 0., 100000.);
  full_fit->SetLineColor(kRed);
  full_fit->SetLineWidth(2);
  TF1 *Sin = new TF1("Sin", "[0]+[1]*TMath::Sin(x-[2])", -1000, 1000);
  Sin->SetNpx(1.e7);
  Sin->SetParLimits(0, 0., 1000.);
  Sin->SetParLimits(1, 0., 1000.);
  Sin->SetParLimits(2, -TMath::Pi(), TMath::Pi());
  Sin->SetLineColor(kRed);
  Sin->SetLineWidth(2);
  TF1 *Gauss = new TF1("Gauss", "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", -1000, 1000);
  Gauss->SetParLimits(0, 0., 100000.);
  Gauss->SetParLimits(1, 0., 1000.);
  Gauss->SetParLimits(3, 0., 1000.);
  Gauss->SetParLimits(4, 0., 100000.);
  Gauss->SetNpx(1.e7);
  Gauss->SetLineColor(kRed);
  Gauss->SetLineWidth(2);

  auto continue_loop = true;
  while (continue_loop)
  {
    continue_loop = false;

    //  Fit Profile (only sin dependence)
    cCheckResults->cd(1);
    auto efficiency_profile = hTarget->ProfileX("efficiency_profile", -1, -1);
    Sin->SetParameters(efficiency_profile->GetMean(2), 1., 1.);
    efficiency_profile->Fit(Sin);

    //  Fit slice for gaus guesstimate
    auto select_slice_xBin = 0;
    auto select_slice_entries = 0;
    for (auto xBin = 1; xBin <= hTarget->GetNbinsX(); xBin++)
    {
      auto current_slice = hTarget->ProjectionY(Form("current_slice_%i", xBin), xBin, xBin);
      if (current_slice->GetEntries() > select_slice_entries)
      {
        select_slice_xBin = xBin;
        select_slice_entries = current_slice->GetEntries();
      }
    }
    auto fit_slice = hTarget->ProjectionY(Form("current_slice_%i", select_slice_xBin), select_slice_xBin, select_slice_xBin);
    Gauss->SetParameters(select_slice_entries, fit_slice->GetMean(), fit_slice->GetRMS(), select_slice_entries * 0.1);
    tmp->cd();
    fit_slice->Fit(Gauss);

    // Assign guesstimates to full fit
    cCheckResults->cd(2);
    full_fit->FixParameter(0, Sin->GetParameter(0));
    full_fit->FixParameter(1, Sin->GetParameter(1));
    full_fit->FixParameter(2, Sin->GetParameter(2));
    full_fit->SetParameter(3, Gauss->GetParameter(0));
    full_fit->SetParameter(4, Gauss->GetParameter(2));
    full_fit->SetParameter(5, Gauss->GetParameter(3));
    hTarget->Fit(full_fit);
    full_fit->ReleaseParameter(0);
    full_fit->ReleaseParameter(1);
    full_fit->ReleaseParameter(2);
    hTarget->Fit(full_fit);

    //  Check the sin fit
    cCheckResults->cd(1);
    Sin->SetParameters(full_fit->GetParameter(0), full_fit->GetParameter(1), full_fit->GetParameter(2));
    Sin->SetLineColor(kBlue);
    Sin->SetLineStyle(kDashed);
    Sin->SetLineWidth(2);
    Sin->DrawCopy("SAMEL");

    for (auto xBin = 1; xBin <= hTarget->GetNbinsX(); xBin++)
    {
      auto current_x = hTarget->GetXaxis()->GetBinCenter(xBin);
      auto current_y = Sin->Eval(current_x);
      auto current_z = full_fit->Eval(current_x, current_y);
      auto current_slice = hTarget->ProjectionY(Form("current_slice_%i", xBin), xBin, xBin);
      auto current_peak = current_slice->GetBinContent(current_slice->GetMaximumBin());
      if (current_peak < current_z * 0.1)
      {
        if (current_peak > 0)
          continue_loop = true;
        for (auto yBin = 1; yBin <= hTarget->GetNbinsY(); yBin++)
          hTarget->SetBinContent(xBin, yBin, 0.);
      }
    }
    if (!continue_loop)
    {
      //  Check the sin fit
      cCheckResults->cd(1);
      Sin->SetParameters(full_fit->GetParameter(0), full_fit->GetParameter(1), full_fit->GetParameter(2));
      Sin->SetLineColor(kRed);
      Sin->SetLineStyle(kDashed);
      Sin->SetLineWidth(2);
      Sin->DrawCopy("SAMEL");
    }
  }

  auto final_prediction = get_prediction(full_fit);

  TCanvas *cDrawResult = standard_canvas_2D("cDrawResult", "cDrawResult");
  hTarget_input._TH2F["hPersistance2D"]->Draw("COLZ");
  std::array<std::array<int, 3>, 3> plot_options = {{{kBlack, kSolid, 5}, {kRed, kDashed, 3}, {kRed, kSolid, 3}}};
  plot_circle(final_prediction, plot_options);

  TLatex *lLatex = new TLatex();
  lLatex->SetTextSize(0.03);
  lLatex->DrawLatexNDC(0.05, 0.965, Form("X_{0} : %.3f#pm%.3f mm", final_prediction[0][0], final_prediction[0][1]));
  lLatex->DrawLatexNDC(0.30, 0.965, Form("Y_{0} : %.3f#pm%.3f mm", final_prediction[1][0], final_prediction[1][1]));
  lLatex->DrawLatexNDC(0.55, 0.965, Form("R : %.3f#pm%.3f mm", final_prediction[2][0], final_prediction[2][1]));
  lLatex->DrawLatexNDC(0.80, 0.965, Form("#sigma_{R} : %.3f#pm%.3f mm", full_fit->GetParameter(4), full_fit->GetParError(4)));
  return final_prediction;
}

void simulation_QA(std::vector<std::string> list_of_run_tags = {}, std::string output_filename = "out.root")
{
  /*
  cout << "ff : " << get_radius(particle_mass["proton"], 10, 1.02, 353) << endl;
  cout << "ff : " << get_radius(particle_mass["kaon+-"], 10, 1.02, 353) << endl;
  cout << "ff : " << get_radius(particle_mass["pion+-"], 10, 1.02, 353) << endl;
  cout << "ff : " << get_radius(particle_mass["proton"], 10, 1.0008, 353) << endl;
  cout << "ff : " << get_radius(particle_mass["kaon+-"], 10, 1.0008, 353) << endl;
  cout << "ff : " << get_radius(particle_mass["pion+-"], 10, 1.0008, 353) << endl;

  cout << "Expected radius : " << get_radius(particle_mass["pion+-"], 10, 1.02, 353) << endl;
  cout << "Expected ctheta : " << get_ctheta(get_beta(particle_mass["pion+-"], 10), 1.02) << " Calculated ctheta: " << measure_ctheta(get_radius(particle_mass["pion+-"], 10, 1.02, 353), 353) << endl;
  cout << "Expected beta : " << get_beta(particle_mass["pion+-"], 10) << " Calculated beta: " << measure_beta(1.02, get_radius(particle_mass["pion+-"], 10, 1.02, 353), 353) << endl;
  cout << "Expected mass : " << particle_mass["pion+-"] << " Calculated mass: " << get_mass_hypothesis(10., 1.02, get_radius(particle_mass["pion+-"], 10, 1.02, 353), 353) << endl;

  cout << "temperature: " << database::get_temperature<true>("20231010-084623").second << "C " << endl;
  cout << "temperature: " << database::get_temperature<false>("20231010-084623").second << "K " << endl;
  cout << "Energy: " << database::get_beam_energy("20231010-084623").second << endl;
  cout << "Polarity: " << database::get_beam_polarity("20231010-084623").second << endl;
  cout << "Vbias: " << database::get_Vbias("20231010-084623").second << endl;
  cout << "Gas mirror: " << database::get_gas_mirror_z("20231010-084623").second << endl;
  cout << "Aerogel mirror: " << database::get_gas_aerogel_z("20231010-084623").second << endl;
  cout << "Has gas: " << database::get_has_gas("20231010-084623").second << endl;
  cout << "n_gas: " << database::get_n_gas("20231010-084623").second << endl;
  cout << "n_aerogel: " << database::get_n_aerogel("20231010-084623").second << endl;
*/

  cout << "TEST" << endl;
  list_of_run_tags.push_back("20231010-084623");
  // list_of_run_tags.push_back("MC.real.Nsig=20.Nbkg=15");
  for (auto Nbkg = 0; Nbkg <= 20; Nbkg++)
  {
    // list_of_run_tags.push_back("MC.full.Nsig=20.Nbkg=" + std::to_string(Nbkg));
    // list_of_run_tags.push_back("MC.real.Nsig=20.Nbkg=" + std::to_string(Nbkg));
  }
  cout << "TEST" << endl;

  //  Output
  std::map<std::string, TObjects> inputContainer;
  TObjects outputContainter;
  cout << "TEST" << endl;

  //  Recover raw data histogram
  recover_all_histograms(inputContainer, list_of_run_tags);
  cout << "TEST" << endl;

  //  Compare coordinates
  // cout << compare_bkg_evaluation(inputContainer, {0, -25, 25, -55, -35}, {0, 20, 0, 0.035}, outputContainter, "_SOUTH.pdf") << endl;
  // cout << compare_bkg_evaluation(inputContainer, {0, -25, 25, 35, 55}, {0, 20, 0, 0.035}, outputContainter, "_NORTH.pdf") << endl;
  // cout << compare_bkg_evaluation(inputContainer, {0, -60, -20, -25, 25}, {0, 20, 0, 0.095}, outputContainter, "_EAST.pdf") << endl;
  // cout << compare_bkg_evaluation(inputContainer, {0, 20, 60, -25, 25}, {0, 20, 0, 0.095}, outputContainter, "_WEST.pdf") << endl;
  // plot_squares(inputContainer["20231010-084623"], {{-25, 25, -55, -35}, {-25, 25, 35, 55}, {-60, -20, -25, 25}, {20, 60, -25, 25}}, {-100, 100, -100, 100});
  analyse_efficiency(inputContainer["20231010-084623"]);
  // analyse_efficiency(inputContainer["MC.real.Nsig=20.Nbkg=15"]);

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(outputContainter);
  output_file->Close();
}

void ____()
{
  std::string trg_dir = "/Users/nrubini/Downloads/xnicola_2/";
  for (auto Nbkg = 0; Nbkg <= 20; Nbkg++)
  {
    std::string test_string = "full.Nsig=20.Nbkg=" + std::to_string(Nbkg);
    gROOT->ProcessLine(Form(".! mkdir -p %s/MC.%s", sim_dir.c_str(), test_string.c_str()));
    gROOT->ProcessLine(Form(".! cp %s/hough.%s.root %s/MC.%s/hough.root ", trg_dir.c_str(), test_string.c_str(), sim_dir.c_str(), test_string.c_str()));
    gROOT->ProcessLine(Form(".! cp %s/fastmc.%s.root %s/MC.%s/recodata.root ", trg_dir.c_str(), test_string.c_str(), sim_dir.c_str(), test_string.c_str()));
  }
  for (auto Nbkg = 0; Nbkg <= 20; Nbkg++)
  {
    std::string test_string = "real.Nsig=20.Nbkg=" + std::to_string(Nbkg);
    gROOT->ProcessLine(Form(".! mkdir -p %s/MC.%s", sim_dir.c_str(), test_string.c_str()));
    gROOT->ProcessLine(Form(".! cp %s/hough.%s.root %s/MC.%s/hough.root ", trg_dir.c_str(), test_string.c_str(), sim_dir.c_str(), test_string.c_str()));
    gROOT->ProcessLine(Form(".! cp %s/fastmc.%s.root %s/MC.%s/recodata.root ", trg_dir.c_str(), test_string.c_str(), sim_dir.c_str(), test_string.c_str()));
  }
}