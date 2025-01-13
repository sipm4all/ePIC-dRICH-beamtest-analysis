//  --- --- ---
//  Histogram handling utility
//  author: nicola.rubini@bo.infn.it
//  --- --- ---
#pragma once
//  Include
#include "Database.h"
#include "Utility.h"
#include "../macros/produce_raw_data_for_radius_cut_calibration.C"
//  --- --- ---
//
std::map<std::string, std::string> histogram_lookup_table_position =
    {{"hPersistance2D", "raw_prod"},
     {"hFitRadiusVal", "raw_prod"},
     {"hFitRadius", "raw_prod"},
     {"hFitRadius_graph_x27_Rcut_3.000000", "raw_prod_analysis"}};

void run_analysis_step(std::string run_tag, std::string target_file)
{
  if (target_file == "raw_prod")
    produce_raw_data_for_radius_cut_calibration(run_tag);
}

template <typename plot_type>
plot_type *get_plot(std::string run_tag, std::string target_plot)
{
  // TODO: Return error if not found histo position in database
  std::string input_filename = data_dir + "/" + run_tag + "/" + histogram_lookup_table_position[target_plot] + ".root";
  TFile *input_file = new TFile(input_filename.c_str());
  if (input_file->IsZombie())
  {
    run_analysis_step(run_tag, histogram_lookup_table_position[target_plot]);
    return get_plot<plot_type>(run_tag, target_plot);
  }
  auto target_from_file = (plot_type *)(input_file->Get(target_plot.c_str()));
  // if (thx)
  // target_from_file->SetDirectory(0);
  input_file->Close();
  return target_from_file;
}

TH2F *get_persistance_plot(std::string run_tag)
{
  return get_plot<TH2F>(run_tag, "hPersistance2D");
}

TH2F *get_radius_vs_nsigma_plot(std::string run_tag)
{
  return get_plot<TH2F>(run_tag, "hFitRadiusVal");
}

TH3F *get_radius_vs_nsigma_vs_photons_plot(std::string run_tag)
{
  return get_plot<TH3F>(run_tag, "hFitRadius");
}

TGraphErrors *get_Rres_vs_photons_plot(std::string run_tag)
{
  return get_plot<TGraphErrors>(run_tag, "hFitRadius_graph_x27_Rcut_3.000000");
}