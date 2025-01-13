#include "../lib/Utility.h"
#include "../lib/Database.h"
#include "../macros/run_pre_processing.C"
#include "../macros/analyse_raw_data_for_radius_cut_calibration.C"

std::vector<std::string> fields = {"run", "quantity"};

std::vector<std::array<std::string, 2>> read_file(std::string input_file)
{
  // Output
  std::vector<std::array<std::string, 2>> result;

  //  Start reading the file
  std::ifstream data_stream(input_file);
  std::string current_line;
  auto iLine = -1;
  while (std::getline(data_stream, current_line))
  {
    //  Skip comment characters
    if (current_line[0] == '#' || current_line[0] == ' ')
      continue;

    // Increment line counter
    iLine++;
    result.push_back({});

    //  Read database
    std::stringstream string_in_stream(current_line);
    std::string current_data;
    auto iField = -1;
    for (auto current_field : fields)
    {
      string_in_stream >> current_data;
      iField++;
      result[iLine][iField] = current_data;
    }
  }
  return result;
}

std::array<float, 2> get_radius_sigma(std::string run_tag, float requested_sigma_cut = 3.)
{
  //  Result
  std::array<float, 2> result = {-1, -1};

  //  Open target file
  std::string input_filename = data_dir + "/" + run_tag + "/raw_prod_analysis.root";
  TFile *target_file = new TFile(input_filename.c_str());
  if (target_file->IsZombie() )
  {
    analyse_raw_data_for_radius_cut_calibration(run_tag);
    target_file = new TFile(input_filename.c_str());
  }

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

void run_scan(std::string input_filename = "./lists/mirror_scan_fwd.list", std::string variable_name = "variable (a.u.)")
{
  // retrieve run list and
  auto run_data = read_file(input_filename);

  // Output graph
  TGraphErrors *result_graph = new TGraphErrors();

  //  Loop over runs requested
  for (auto current_run_info : run_data)
  {
    TFile *current_run_pp_file = new TFile(Form("%s/%s/pre_process_recodata.root", data_dir.c_str(), current_run_info[0].c_str()));
    if (current_run_pp_file->IsZombie() )
    {
      run_pre_processing(current_run_info[0]);
      current_run_pp_file = new TFile(Form("%s/%s/pre_process_recodata.root", data_dir.c_str(), current_run_info[0].c_str()));
    }
    delete current_run_pp_file;
    auto current_resolution = get_radius_sigma(current_run_info[0]);
    auto current_point = result_graph->GetN();
    result_graph->SetPoint(current_point, std::stof(current_run_info[1]), current_resolution[0]);
    result_graph->SetPointError(current_point, 0, current_resolution[1]);
  }

  TCanvas *cDrawResult = new TCanvas();
  result_graph->GetXaxis()->SetTitle(variable_name.c_str());
  result_graph->GetYaxis()->SetTitle("Resolution (mm)");
  result_graph->Draw("ALPE");
}