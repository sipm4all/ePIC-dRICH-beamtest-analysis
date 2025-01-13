#include "../lib/Utility.h"
#include "../lib/Database.h"
#include "../lib/mapping.h"


void test_macro()
{
  //  Output & Utility
  TObjects output_struct;
  TObjects utility_struct;
  output_struct._TH2F["hMap_selected_rings"] = new TH2F("hMap_selected_rings", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_available_SiPM"] = new TH2F("hMap_available_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_fullsetup_SiPM"] = new TH2F("hMap_fullsetup_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_available_SiPM_acceptance"] = new TH2F("hMap_available_SiPM_acceptance", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
  output_struct._TH2F["hMap_fullsetup_SiPM_acceptance"] = new TH2F("hMap_fullsetup_SiPM_acceptance", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);

  //  Load recodata
  std::string run_tag = "20231010-084623";
  std::string inFile_recodata = data_dir + "/" + run_tag + "/recodata.root";
  recodata reco_data;
  auto reco_tree = load_data(inFile_recodata, reco_data);

  //   Loop on events
  for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
  {
    //  Recover recodata entry form tree
    reco_tree->GetEntry(iEv);

    if (iEv % 200 == 0)
      cout << "iEv: " << iEv << endl;

    if (iEv == 200)
      break;

    //  Loop on hits
    for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    {
      //  Add sensor sensitive area
      fill_with_SiPM_coverage(output_struct._TH2F["map_of_SiPM_coverage_data"], {reco_data.x[iPnt], reco_data.y[iPnt]});
    }
  }

  for (auto iPDU = 0; iPDU < 8; iPDU++)
  {
    for (auto iCol = 0; iCol < 16; iCol++)
    {
      for (auto iRow = 0; iRow < 16; iRow++)
      {
        fill_with_SiPM_coverage(output_struct._TH2F["map_of_SiPM_coverage"], sipm4eic::get_position({iPDU, iCol, iRow}));
      }
    }
  }

  cout << calculate_efficiency(output_struct._TH2F["map_of_SiPM_coverage"], {2.86829, 1.18674, 73.0474}, 5 * 1.8854, output_struct._TH2F["map_of_acceptance"], output_struct._TH2F["map_of_ring_coverage"]) << endl;
  cout << calculate_efficiency(output_struct._TH2F["map_of_SiPM_coverage_data"], {2.86829, 1.18674, 73.0474}, 5 * 1.8854, output_struct._TH2F["map_of_acceptance_data"], output_struct._TH2F["map_of_ring_coverage"]) << endl;

  TCanvas *c1 = new TCanvas("", "", 1200, 1200);
  gStyle->SetOptStat(0);
  c1->Divide(2, 2);
  c1->cd(2);
  output_struct._TH2F["map_of_acceptance"]->Draw("COLZ");
  c1->cd(1);
  output_struct._TH2F["map_of_SiPM_coverage"]->Draw("COLZ");
  c1->cd(3);
  output_struct._TH2F["map_of_SiPM_coverage_data"]->Draw("COLZ");
  c1->cd(4);
  output_struct._TH2F["map_of_acceptance_data"]->Draw("COLZ");

}

/*
void test_macro(std::string run_list = "mirror_scan_bkw", std::string output_filename = "")
{
  //  Output
  TObjects output_struct;
  TObjects utility_struct;

  //  Input/Output files
  std::string input_filename_template = data_dir + "/%s/raw_prod_analysis.root";
  if (output_filename.empty())
    output_filename = data_dir + "/analysis_" + run_list + ".root";

  //  Set Batch
  gROOT->SetBatch(kTRUE);

  //  Run pre-process on all runs
  for (auto current_run_tag : database::run_lists[run_list])
  {
    //continue;
    produce_raw_data_for_radius_cut_calibration(current_run_tag);
    analyse_raw_data_for_radius_cut_calibration(current_run_tag);
  }

  //  Set Batch
  gROOT->SetBatch(kFALSE);

  for (auto current_run_tag : database::run_lists[run_list])
  {
    TFile *current_input_file = new TFile(Form(input_filename_template.c_str(), current_run_tag.c_str()));
    output_struct._TGraphErrors[Form("mirror_at_%.0fmm_%s", database::get_aerogel_mirror_z(current_run_tag).second,current_run_tag.c_str())] = (TGraphErrors *)(current_input_file->Get("hFitRadius_Full_from_fit"));
    current_input_file->Close();
  }

  //  Parabola fit
  TF1 *parabola_fit = new TF1("parabola_fit", "[0]*x*(x-[1])+[2]");

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  write_all(output_struct);
  output_file->Close();
}
*/