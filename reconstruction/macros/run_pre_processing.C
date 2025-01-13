#pragma once
#include "../lib/Utility.h"
#include "../lib/mapping.h"

template <bool fix_center = true>
void run_pre_processing(std::string run_tag = "20231010-084623", std::string output_filename = "")
{
  //  Output & Utility
  TObjects output_struct;
  TObjects utility_struct;

  //  Input/Output files
  bool running_on_MC = !run_tag.rfind("MC");
  std::string inFile_recodata = data_dir + "/" + run_tag + "/recodata.root";
  //  ->  Legacy hough circle finder
  // std::string inFile_ringdata = data_dir + "/" + run_tag + "/hough.root";
  if (running_on_MC)
  {
    inFile_recodata = sim_dir + "/" + run_tag + "/recodata.root";
    //  ->  Legacy hough circle finder
    // inFile_ringdata = sim_dir + "/" + run_tag + "/hough.root";
    if (output_filename.empty())
      output_filename = sim_dir + "/" + run_tag + "/pre_process_recodata.root";
  }
  if (output_filename.empty())
    output_filename = data_dir + "/" + run_tag + "/pre_process_recodata.root";

  //  Link TTree to local data instance
  recodata reco_data;
  check_recodata(run_tag);
  auto reco_tree = load_data(inFile_recodata, reco_data);

  //  Requested cuts on the radius
  auto ntests = 24;
  auto start = 1.5;
  auto stop = 7.5;
  auto bin_size = (1. / ntests) * (stop - start);

  //  ==  Output objects  ==
  //  ==  ==  Ring(s) discovery & QA
  //  ==  ==  TODO: Generalise w/ a 2D histogram to make a free choice of slices to use
  output_struct._TH1F["hXcoordinateFinder_0"] = new TH1F("hXcoordinateFinder_0", ";X (mm)", 99, -99, 99);
  output_struct._TH1F["hYcoordinateFinder_0"] = new TH1F("hYcoordinateFinder_0", ";Y (mm)", 99, -99, 99);
  output_struct._TH1F["hXcoordinateFinder_1"] = new TH1F("hXcoordinateFinder_1", ";X (mm)", 99, -99, 99);
  output_struct._TH1F["hYcoordinateFinder_1"] = new TH1F("hYcoordinateFinder_1", ";Y (mm)", 99, -99, 99);
  output_struct._TH1F["hXcoordinateFinder_2"] = new TH1F("hXcoordinateFinder_2", ";X (mm)", 99, -99, 99);
  output_struct._TH1F["hYcoordinateFinder_2"] = new TH1F("hYcoordinateFinder_2", ";Y (mm)", 99, -99, 99);
  output_struct._TH1F["hRadius_check_diffcenters"] = new TH1F("hRadius_check_diffcenters", ";R (mm)", 99, 0, 99);
  output_struct._TH1F["hRadius_check_samecenters"] = new TH1F("hRadius_check_samecenters", ";R (mm)", 99, 0, 99);
  output_struct._TH1F["hFirstGuessCoordinates"]; // Defined below
  //  ==  == Persistance maps and acceptance maps
  output_struct._TH2F["hPersistance2D"] = new TH2F("hPersistance2D", ";X (mm);Y (mm)", 396, -99, 99, 396, -99, 99);
  output_struct._TH2F["hMap_selected_rings"] = new TH2F("hMap_selected_rings", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_available_SiPM"] = new TH2F("hMap_available_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_fullsetup_SiPM"] = new TH2F("hMap_fullsetup_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_available_SiPM_acceptance"] = new TH2F("hMap_available_SiPM_acceptance", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_fullsetup_SiPM_acceptance"] = new TH2F("hMap_fullsetup_SiPM_acceptance", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);

  //  --- Radius analysis
  output_struct._TH3F["hFitRadius"] = new TH3F("hFitRadius", ";Radius cut (n#sigma);n_{#gamma};R from fit (mm)", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 30, 3.5, 33.5, 300, 40, 100);
  output_struct._TH3F["hDeltaRadius"] = new TH3F("hDeltaRadius", ";Radius cut (n#sigma);n_{#gamma};#Delta R (mm)", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 30, 3.5, 33.5, 500, -20, 20);
  //  --- QA
  //  --- --- Signal selection
  output_struct._TH2F["hAveragePhotonEv"] = new TH2F("hAveragePhotonEv", ";Radius cut (n#sigma);Selected photons", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 40, 0, 40);
  output_struct._TH2F["hAveragePhotonEvDisc"] = new TH2F("hAveragePhotonEvDisc", ";Radius cut (n#sigma);Discarded photons", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 40, 0, 40);
  //  --- --- Coordinates
  output_struct._TH2F["hFitXcoordinate"] = new TH2F("hFitXcoordinate", ";Radius cut (n#sigma);X_{0} (mm)", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 1000, -10, 10);
  output_struct._TH2F["hFitYcoordinate"] = new TH2F("hFitYcoordinate", ";Radius cut (n#sigma);Y_{0} (mm)", ntests + 1, start - 0.5 * bin_size, stop + 0.5 * bin_size, 1000, -10, 10);
  //  --- --- Correlations
  output_struct._TH2F["hCorrelation_on_sel"] = new TH2F("hCorrelation_on_sel", ";Selected photons;#varphi (rad)", 36, -0.5, 35.5, 72, -TMath::Pi(), TMath::Pi());
  output_struct._TH2F["hCorrelation_on_sel_test"] = new TH2F("hCorrelation_on_sel_test", ";Selected photons;#varphi (rad)", 4, -0.5, 35.5, 72, -TMath::Pi(), TMath::Pi());
  //  --- --- Sort me
  output_struct._TH1F["hTimingFull"] = new TH1F("hTimingFull", ";X (mm);Y (mm)", 99, -99, 99);
  output_struct._TH1F["hTimingLow_"] = new TH1F("hTimingLow_", ";X (mm);Y (mm)", 99, -99, 99);
  output_struct._TH1F["hTimingHigh"] = new TH1F("hTimingHigh", ";X (mm);Y (mm)", 99, -99, 99);
  output_struct._TH2F["hPolarEfficiencyMap"] = new TH2F("hEfficiencyMap", ";#phi (rad); R (mm)", 360, -TMath::Pi(), TMath::Pi(), 150, 0, 150);
  output_struct._TH2F["hPersistance2D_Nphoton_nocut"] = new TH2F("hPersistance2D_Nphoton_nocut", ";X (mm);Y (mm)", 99, -99, 99, 99, -99, 99);
  output_struct._TH2F["hPersistance2D_Nphoton_lowcut"] = new TH2F("hPersistance2D_Nphoton_lowcut", ";X (mm);Y (mm)", 99, -99, 99, 99, -99, 99);
  output_struct._TH2F["hPersistance2D_Nphoton_higcut"] = new TH2F("hPersistance2D_Nphoton_higcut", ";X (mm);Y (mm)", 99, -99, 99, 99, -99, 99);
  output_struct._TH3F["hPersistance3D"] = new TH3F("hPersistance3D", ";X (mm);Y (mm); t (ns)", 396, -99, 99, 396, -99, 99, 160, -80, 80);

  std::vector<std::array<float, 2>> list_of_available_SiPMs;

  //  First loop on events
  cout << "[INFO] Start of preliminary loop for X_{0}, Y_{0} and R_{0}" << endl;
  for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
  {
    //  Recover recodata entry form tree
    reco_tree->GetEntry(iEv);

    if (iEv % 1000 == 0)
      cout << "[INFO] event: " << iEv << endl;

    //  Fill the X-Y finder of central point
    //  TODO: Generalise with string in input fill_coordinate_x
    fill_coordinate_x(output_struct._TH1F["hXcoordinateFinder_0"], reco_data, -3, 0);
    fill_coordinate_y(output_struct._TH1F["hYcoordinateFinder_0"], reco_data, -3, 0);
    fill_coordinate_x(output_struct._TH1F["hXcoordinateFinder_1"], reco_data, -10, -7);
    fill_coordinate_y(output_struct._TH1F["hYcoordinateFinder_1"], reco_data, -10, -7);
    fill_coordinate_x(output_struct._TH1F["hXcoordinateFinder_2"], reco_data, 7, 10);
    fill_coordinate_y(output_struct._TH1F["hYcoordinateFinder_2"], reco_data, 7, 10);

    //  Persistance plot
    fill_persistance(output_struct._TH2F["hPersistance2D"], reco_data);
    fill_persistance(output_struct._TH3F["hPersistance3D"], reco_data);

    //  Efficiency plot
    fill_polar_efficiency(reco_data, output_struct._TH2F["hPolarEfficiencyMap"]);

    //  Loop on hits
    for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    {
      if (!count(list_of_available_SiPMs.begin(), list_of_available_SiPMs.end(), std::array<float, 2>({reco_data.x[iPnt], reco_data.y[iPnt]})))
        //  Add sensor sensitive area
        list_of_available_SiPMs.push_back({reco_data.x[iPnt], reco_data.y[iPnt]});
    }
  }

  for (auto iPDU = 0; iPDU < 8; iPDU++)
  {
    for (auto iCol = 0; iCol < 16; iCol++)
    {
      for (auto iRow = 0; iRow < 16; iRow++)
      {
        fill_with_SiPM_coverage(output_struct._TH2F["hMap_fullsetup_SiPM"], sipm4eic::get_position({iPDU, iCol, iRow}));
      }
    }
  }

  auto found_rings = fit_multiple_rings<false>({{output_struct._TH1F["hXcoordinateFinder_0"], output_struct._TH1F["hYcoordinateFinder_0"]}, {output_struct._TH1F["hXcoordinateFinder_1"], output_struct._TH1F["hYcoordinateFinder_1"]}, {output_struct._TH1F["hXcoordinateFinder_2"], output_struct._TH1F["hYcoordinateFinder_2"]}}, {{-1.5, -1.5}, {-8.5, -8.5}, {8.5, 8.5}});
  auto check_canvas = plot_check_coordinates(output_struct._TH2F["hPersistance2D"], found_rings);
  check_canvas->SaveAs(Form("%s/%s/check_coordinate_finding.pdf", data_dir.c_str(), run_tag.c_str()));

  //  Store initial guesses
  output_struct._TH1F["hFirstGuessCoordinates"] = new TH1F("hRadius_check_samecenters", ";R (mm)", found_rings.size() * 4, 0, found_rings.size() * 4);

  //  Second loop on events
  cout << "[INFO] Start of radius check loop for X_{0} and Y_{0}" << endl;
  auto iRing = -1;
  for (auto current_ring : found_rings)
  {
    iRing++;
    for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
    {
      //  Recover recodata entry form tree
      reco_tree->GetEntry(iEv);

      if (iEv % 1000 == 0)
        cout << "[INFO] event: " << iEv << endl;
      std::array<float, 3> ring_coordinates = {current_ring[0][0], current_ring[1][0], current_ring[2][0]};
      auto current_ring_region = select_points(reco_data, ring_coordinates, 10.);
      fill_coordinate_R(output_struct._TH1F["hRadius_check_samecenters"], current_ring_region);
      fill_coordinate_R(output_struct._TH1F["hRadius_check_diffcenters"], current_ring_region, {current_ring[0][0], current_ring[1][0]});
    }

    TF1 *Gauss = new TF1("Gauss", "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", -1000, 1000);
    Gauss->SetParLimits(0, 0., 100000.);
    Gauss->SetParLimits(1, 0., 100.);
    Gauss->SetParLimits(2, 0., 10.);
    Gauss->SetParLimits(3, 0., 10000.);
    Gauss->SetParameter(1, current_ring[2][0]);
    Gauss->SetParameter(2, 1.);
    output_struct._TH1F["hRadius_check_diffcenters"]->Fit(Gauss, "IMRES", "", current_ring[2][0] - 10, current_ring[2][0] + 10);

    output_struct._TH1F["hFirstGuessCoordinates"]->SetBinContent(iRing * 4 + 1, current_ring[0][0]);
    output_struct._TH1F["hFirstGuessCoordinates"]->SetBinContent(iRing * 4 + 2, current_ring[1][0]);
    output_struct._TH1F["hFirstGuessCoordinates"]->SetBinContent(iRing * 4 + 3, current_ring[2][0]);
    output_struct._TH1F["hFirstGuessCoordinates"]->SetBinContent(iRing * 4 + 4, Gauss->GetParameter(2));
  }

  //  Final loop on events
  cout << "[INFO] Start of final loop for X_{0} and Y_{0}" << endl;
  for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
  {
    //  Recover entry form tree
    reco_tree->GetEntry(iEv);

    if (iEv % 1000 == 0)
      cout << "[INFO] event: " << iEv << endl;

    //  Loop on requested cuts
    auto iRing = -1;
    for (auto current_ring : found_rings)
    {
      iRing++;
      for (auto iCut = 0; iCut <= ntests; iCut++)
      {
        //  Calculate cut and fit
        std::array<float, 3> ring_data = {current_ring[0][0], current_ring[1][0], current_ring[2][0]};
        auto current_cut = start + iCut * bin_size;
        auto circle_best_fit = get_best_circle(reco_data, ring_data, current_cut);

        //  Exclude errors
        if (circle_best_fit[2][0] < 0)
          continue;

        //  Calculate info
        auto radius_delta = get_radius_delta(circle_best_fit, reco_data);
        auto selected_points = select_points(reco_data, ring_data, current_cut);
        auto discarded_points = select_points<false>(reco_data, ring_data, current_cut);
        auto n_selected_photons = selected_points.n;
        auto n_discarded_photons = discarded_points.n;

        //  Fill analysis histograms
        for (auto current_delta : radius_delta)
          output_struct._TH3F["hDeltaRadius"]->Fill(current_cut, n_selected_photons, current_delta);
        output_struct._TH3F["hFitRadius"]->Fill(current_cut, n_selected_photons, circle_best_fit[2][0]);
        output_struct._TH2F["hAveragePhotonEv"]->Fill(current_cut, n_selected_photons);
        output_struct._TH2F["hAveragePhotonEvDisc"]->Fill(current_cut, n_discarded_photons);
      }
    }
  }

  // Make acceptance
  // calculate_efficiency(output_struct._TH2F["hMap_available_SiPM"], ring_data, 3 * run_coordinates[3][0], output_struct._TH2F["hMap_available_SiPM_acceptance"], output_struct._TH2F["hMap_selected_rings"]);
  // calculate_efficiency(output_struct._TH2F["hMap_fullsetup_SiPM"], ring_data, 3 * run_coordinates[3][0], output_struct._TH2F["hMap_fullsetup_SiPM_acceptance"], output_struct._TH2F["hMap_selected_rings"]);

  output_struct._TH1F["hCorrelation_on_sel_test_1"] = (TH1F *)(output_struct._TH2F["hCorrelation_on_sel_test"]->ProjectionY("tmp1", 1, 1));
  output_struct._TH1F["hCorrelation_on_sel_test_1"]->Scale(1. / (output_struct._TH1F["hCorrelation_on_sel_test_1"]->GetEntries()));
  output_struct._TH1F["hCorrelation_on_sel_test_2"] = (TH1F *)(output_struct._TH2F["hCorrelation_on_sel_test"]->ProjectionY("tmp2", 2, 2));
  output_struct._TH1F["hCorrelation_on_sel_test_2"]->Scale(1. / (output_struct._TH1F["hCorrelation_on_sel_test_2"]->GetEntries()));
  output_struct._TH1F["hCorrelation_on_sel_test_3"] = (TH1F *)(output_struct._TH2F["hCorrelation_on_sel_test"]->ProjectionY("tmp3", 3, 3));
  output_struct._TH1F["hCorrelation_on_sel_test_3"]->Scale(1. / (output_struct._TH1F["hCorrelation_on_sel_test_3"]->GetEntries()));
  output_struct._TH1F["hCorrelation_on_sel_test_4"] = (TH1F *)(output_struct._TH2F["hCorrelation_on_sel_test"]->ProjectionY("tmp4", 4, 4));
  output_struct._TH1F["hCorrelation_on_sel_test_4"]->Scale(1. / (output_struct._TH1F["hCorrelation_on_sel_test_4"]->GetEntries()));

  //  Save to file
  TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
  write_all(output_struct);
  outfile->Close();

  return;
}