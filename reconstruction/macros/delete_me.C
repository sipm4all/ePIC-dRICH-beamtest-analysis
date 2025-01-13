#pragma once
#include "../lib/Utility.h"
#include "../lib/Histogram_handling.h"
#include "../macros/produce_raw_data_for_radius_cut_calibration.C"
#include "../macros/analyse_raw_data_for_radius_cut_calibration.C"

void delete_me(std::string run_tag = "20231011-015609", std::string output_filename = "")
{
  //  Output & Utility
  TObjects output_struct;
  TObjects utility_struct;

  //  Input/Output files
  bool running_on_MC = !run_tag.rfind("MC");
  std::string inFile_recodata = data_dir + "/" + run_tag + "/recodata.root";
  std::string inFile_ringdata = data_dir + "/" + run_tag + "/hough.root";
  if (running_on_MC)
  {
    inFile_recodata = sim_dir + "/" + run_tag + "/recodata.root";
    inFile_ringdata = sim_dir + "/" + run_tag + "/hough.root";
  }

  //  Link TTree to local data instance
  recodata reco_data;
  auto reco_tree = load_data(inFile_recodata, reco_data);

  //  produce_raw_data_for_radius_cut_calibration(run_tag);

  //  Create distribution of time
  output_struct._TH1F["delta_time"] = new TH1F("delta_time", "; #Delta t_{hit-trg} (ns);Entries", 51, -80, 80);

  //  First loop on events
  for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
  {
    //  Recover recodata entry form tree
    reco_tree->GetEntry(iEv);

    for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    {
      //  Time distribution
      output_struct._TH1F["delta_time"]->Fill(reco_data.t[iPnt]);
    }
  }
  analyse_raw_data_for_radius_cut_calibration(run_tag);

  /*
    TCanvas *c1 = plot_standard_canvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    output_struct._TH1F["delta_time"]->Draw();
    output_struct._TH1F["delta_time"]->GetXaxis()->SetTitleSize(0.05);
    output_struct._TH1F["delta_time"]->GetYaxis()->SetTitleSize(0.05);
    c1->SaveAs("c1.pdf");
    delete c1;

    auto plot_radius = get_radius_vs_nsigma_plot(run_tag);
    c1 = standard_canvas_2D<false, false>("dd", "");
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    plot_radius->Draw("COLZ");
    plot_radius->GetYaxis()->SetRangeUser(55, 80);
    plot_radius->GetXaxis()->SetTitleSize(0.05);
    plot_radius->GetYaxis()->SetTitleSize(0.05);
    c1->SaveAs("c2.pdf");
    delete c1;

    auto slice = plot_radius->ProjectionY("", 55, 55);
    c1 = standard_canvas("", "", 1000, 1000);
    gPad->SetLogy();
    slice->Draw("");
    slice->GetXaxis()->SetTitleSize(0.05);
    slice->GetYaxis()->SetTitleSize(0.05);
    c1->SaveAs("c3.pdf");
    delete c1;

    auto plot_radius_2 = get_radius_vs_nsigma_vs_photons_plot(run_tag);
    c1 = standard_canvas_2D<false, false>("dd", "");
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    plot_radius_2->Draw("COLZ");
    plot_radius_2->GetYaxis()->SetRangeUser(55, 80);
    plot_radius_2->GetXaxis()->SetTitleSize(0.05);
    plot_radius_2->GetYaxis()->SetTitleSize(0.05);
    c1->SaveAs("c2.pdf");
    delete c1;

    plot_radius_2->GetXaxis()->SetRange(26, 26);
    plot_radius_2->GetZaxis()->SetRangeUser(65, 75);
    auto slice_2 = (TH2F *)(plot_radius_2->Project3D("ZY"));
    slice_2->SetTitle("");
    c1 = standard_canvas("", "", 1000, 1000);
    gPad->SetLogz();
    slice_2->Draw("COLZ");
    slice_2->GetXaxis()->SetTitleSize(0.05);
    slice_2->GetYaxis()->SetTitleSize(0.05);
    c1->SaveAs("c3.pdf");
    delete c1;

    c1 = standard_canvas("", "", 1000, 1000);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    // gPad->SetLogz();
    auto slice_44 = slice_2->ProjectionY("", 3, 3);
    slice_44->Draw("");
    slice_44->GetYaxis()->SetRangeUser(0, 450);
    slice_44->GetXaxis()->SetTitleSize(0.05);
    slice_44->GetYaxis()->SetTitleSize(0.05);
    slice_44->Fit("gaus");
    c1->SaveAs("c4.pdf");
    gPad->Modified(); gPad->Update();
    delete c1;

    c1 = standard_canvas("", "", 1000, 1000);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    // gPad->SetLogz();
    slice_44 = slice_2->ProjectionY("", 4, 4);
    slice_44->Draw("");
    slice_44->GetYaxis()->SetRangeUser(0, 450);
    slice_44->GetXaxis()->SetTitleSize(0.05);
    slice_44->GetYaxis()->SetTitleSize(0.05);
    slice_44->Fit("gaus");
    c1->SaveAs("c5.pdf");
    gPad->Modified(); gPad->Update();
    delete c1;

    c1 = standard_canvas("", "", 1000, 1000);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    // gPad->SetLogz();
    slice_44 = slice_2->ProjectionY("", 5, 5);
    slice_44->Draw("");
    slice_44->GetYaxis()->SetRangeUser(0, 450);
    slice_44->GetXaxis()->SetTitleSize(0.05);
    slice_44->GetYaxis()->SetTitleSize(0.05);
    slice_44->Fit("gaus");
    gPad->Modified(); gPad->Update();
    c1->SaveAs("c6.pdf");
    delete c1;
    */

  auto SPR = get_Rres_vs_photons_plot(run_tag);
  auto c1 = standard_canvas_2D<false, false>("dd", "");
  TLatex * lLatex = new TLatex();
  gStyle->SetOptStat(0);
  gPad->SetLogz();
  SPR->Draw("ALPE");
  SPR->SetMarkerStyle(20);
  SPR->SetMarkerSize(2);
  SPR->GetYaxis()->SetRangeUser(0.4, 1.2);
  SPR->GetXaxis()->SetTitleSize(0.05);
  SPR->GetYaxis()->SetTitleSize(0.05);
  SPR->GetXaxis()->SetTitle("Number of photons (N_{#gamma})");
  SPR->GetYaxis()->SetTitle("Radius resolution (mm)");
  lLatex->DrawLatexNDC(0.45,0.9,Form("#sigma_{SPR} : %.3f #pm %.3f",2.0263967,0.0090082669));
  c1->SaveAs("c2.pdf");
  delete c1;
}