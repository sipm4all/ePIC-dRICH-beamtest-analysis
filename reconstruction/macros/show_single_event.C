#include "../lib/Utility.h"
#include "../lib/Database.h"
#include "../lib/Histogram_handling.h"
#include "../macros/recowriter.C"
//#include "../macros/produce_raw_data_for_radius_cut_calibration.C"

void show_single_event(std::string run_tag = "20231011-015609", int i_event = 28519)
{
  //  Output
  TObjects output_struct;
  TObjects utility_struct;

  //  Recover recodata
  std::string inFile_lightdata = data_dir + "/" + run_tag + "/lightdata.root";
  std::string inFile_recodata = data_dir + "/" + run_tag + "/recodata.root";
  // recowriter(inFile_lightdata, inFile_recodata);
  recodata reco_data;
  auto reco_tree = load_data(inFile_recodata, reco_data);

  reco_tree->GetEntry(i_event);
  // produce_raw_data_for_radius_cut_calibration(run_tag);

  //  Plot Graph
  auto current_graph = plot(reco_data);
  current_graph->SetMarkerSize(2.5);
  current_graph->SetMarkerColor(kBlack);
  current_graph->SetMarkerStyle(20);
  auto current_graph_2 = plot(reco_data);
  current_graph_2->SetMarkerSize(2.5);
  current_graph_2->SetMarkerColor(kWhite);
  current_graph_2->SetMarkerStyle(89);

  //  Plot Graph with selected points
  std::array<float, 3> run_coordinates = {4.432, -0.323, 71.544};
  auto current_graph_cut = plot(select_points(reco_data, run_coordinates, 6 * 2.359));
  current_graph_cut->SetMarkerSize(2.5);
  current_graph_cut->SetMarkerColor(2);
  current_graph_cut->SetMarkerStyle(89);

  //  Plot the fit
  auto ring_fit_rslt = plot(get_best_circle(reco_data, run_coordinates, 6 * 2.359));
  ring_fit_rslt->SetLineWidth(4);

  // Plot map
  TCanvas *c1 = standard_canvas_2D("", "");
  get_persistance_plot(run_tag)->Draw("COLZ");
  current_graph->Draw("SAMEPE");
  current_graph_2->Draw("SAMEPE");
  current_graph_cut->Draw("SAMEPE");
  ring_fit_rslt->Draw("SAME");
  c1->SaveAs(Form("%s_%i.pdf", run_tag.c_str(), i_event));
}