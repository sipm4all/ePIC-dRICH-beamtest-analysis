#pragma once

//  TODO make setters and getters for the parameter of the clustering algorithm
//  TODO cache value for cluster etc.
//  TODO Export DBscan in a secondary class

using fit_circle_result = std::array<std::array<float, 2>, 3>;

namespace testbeam
{
  std::vector<std::array<std::array<float, 2>, 2>> measure_resolution_in_mult_photons(TH2F *radii_vs_photons, int min_entries);

  //  --- Data structure
  struct hit
  {
    UShort_t n;
    float x[256], y[256], t[256]; // Arrays for x, y, t values
  };

  class data
  {

  public:
    // Variables
    testbeam::hit current_hit;
    TTree *current_tree = nullptr;
    TFile *current_file = nullptr;
    float clustering_point_max_distance;
    float clustering_r_point_max_distance;
    float clustering_cluster_max_distance;
    std::array<float, 2> common_center;
    std::array<float, 2> timing_center_sigma;
    std::vector<float> common_radii;
    std::map<std::array<float, 2>, int> available_SiPMs;

    //  --- Getters ---
    UShort_t get_n() const { return current_hit.n; }
    float get_x(int index) const { return current_hit.x[index]; }
    float get_y(int index) const { return current_hit.y[index]; }
    float get_t(int index) const { return current_hit.t[index]; }
    float get_radius(int index, std::array<float, 2> center = {0., 0.}) const;
    float get_average_radius(std::vector<int> cluster, std::array<float, 2> center = {0., 0.}) const;
    TTree *get_tree() const { return current_tree; }
    TFile *get_file() const { return current_file; }
    float get_clustering_point_max_distance() const { return clustering_point_max_distance; }
    float get_clustering_cluster_max_distance() const { return clustering_cluster_max_distance; }
    std::array<float, 2> get_common_center() const { return common_center; }
    float get_common_center_x() const { return common_center[0]; }
    float get_common_center_y() const { return common_center[1]; }
    std::vector<float> get_common_radii() const { return common_radii; }
    float get_common_radius(int index = 0) const { return common_radii[index]; }
    std::array<float, 2> get_timing_center_sigma() const { return timing_center_sigma; }
    float get_timing_center() const { return timing_center_sigma[0]; }
    float get_timing_sigma() const { return timing_center_sigma[1]; }
    std::map<std::array<float, 2>, int> get_available_SiPMs() const { return available_SiPMs; }

    //  --- Setters ---
    void set_n(UShort_t value) { current_hit.n = value; }
    void set_x(int index, float value) { current_hit.x[index] = value; }
    void set_y(int index, float value) { current_hit.y[index] = value; }
    void set_t(int index, float value) { current_hit.t[index] = value; }
    void set_tree(TTree *tree) { current_tree = tree; }
    void set_file(TFile *file) { current_file = file; }
    void set_clustering_point_max_distance(float value) { clustering_point_max_distance = value; }
    void set_clustering_cluster_max_distance(float value) { clustering_cluster_max_distance = value; }
    void set_common_center(const std::array<float, 2> &center) { common_center = center; }
    void set_common_center_x(float x_value) { common_center[0] = x_value; }
    void set_common_center_y(float y_value) { common_center[1] = y_value; }
    void set_common_radii(std::vector<float> value) { common_radii = value; }
    void set_timing_center_sigma(std::array<float, 2> value) { timing_center_sigma = value; }
    void set_timing_center(float c_value) { timing_center_sigma[0] = c_value; }
    void set_timing_sigma(float s_value) { timing_center_sigma[1] = s_value; }
    void set_available_SiPMs(std::map<std::array<float, 2>, int> value) { available_SiPMs = value; }
    void set_available_SiPM(std::array<float, 2> value) { available_SiPMs[value] = 1; }

    //  --- Selection
    std::array<std::vector<int>, 2> select_points(std::array<float, 4> center, std::array<float, 4> tolerance);
    bool is_dcr(int index) { return fabs(get_t(index) - get_timing_center()) > 4 * get_timing_sigma(); }

    //  --- Fit
    fit_circle_result fit_circle(std::vector<int> selected_points, std::array<float, 3> initial_values = {0., 0., 0.}, bool fix_XY = false);
    fit_circle_result fit_circle(TGraphErrors *target_graph, std::array<float, 3> initial_values = {0., 0., 0.}, bool fix_XY = false);
    std::vector<float> diff_circle(std::vector<int> selected_points, std::array<float, 3> initial_values = {0., 0., 0.}, bool fix_XY = false);

    //  --- Data I/O ---
    std::pair<TTree *, TFile *> load_testbeam_data_tree(const std::string &infilename);

    //  --- DBSCAN
    float get_distance(int i_index, int j_index, std::array<float, 3> distance_weight = {1., 1., 1.});
    float get_r_distance(int i_index, int j_index, std::array<float, 2> center = {0., 0.});
    float get_distance(std::vector<int> first_cluster, std::vector<int> second_cluster) { return fabs(get_average_radius(first_cluster) - get_average_radius(second_cluster)); }
    void dfs_search(int node, std::map<int, std::map<int, float>> &admissible_pairs, std::set<int> &visited, std::vector<int> &cluster);
    std::vector<std::vector<int>> find_clusters(std::map<int, std::map<int, float>> &admissible_pairs);
    std::vector<std::vector<int>> get_clusters(float max_distance = -1);
    std::vector<std::vector<int>> get_r_clusters(float max_distance = -1);
    std::vector<std::vector<int>> merge_clusters(float max_distance = -1);
    std::vector<std::vector<int>> get_merged_clusters();
    std::vector<std::vector<int>> get_ring_seed();

    //  --- Plotting tools
    TCanvas *get_standard_canvas();
    TGraphErrors *data_graph(int marker_style = 20, int marker_color = kBlack);
    TGraphErrors *data_graph(std::vector<int> selected_points, int marker_style = 24, int marker_color = kRed);
    TGraphErrors *data_graph(std::vector<std::vector<int>> selected_points, int marker_style = 24, int marker_color = kRed);
    TGraphErrors *add_data_graph(TGraphErrors *plot_graph, int marker_style = 20, int marker_color = kBlack);
    TGraphErrors *add_data_graph(TGraphErrors *plot_graph, std::vector<int> selected_points, int marker_style = 24, int marker_color = kRed);
    TGraphErrors *add_data_graph(TGraphErrors *plot_graph, std::vector<std::vector<int>> selected_points, int marker_style = 24, int marker_color = kRed);
    TEllipse *plot_circle(fit_circle_result parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1);

    //  --- Utility
    void set_up_global_variables();

    //  --- Constructor ---
    data(const std::string &infilename)
    {
      auto current_tree_and_file = load_testbeam_data_tree(infilename);
      current_tree = std::get<0>(current_tree_and_file);
      current_file = std::get<1>(current_tree_and_file);
      clustering_point_max_distance = 8;
      clustering_cluster_max_distance = 8;
      clustering_r_point_max_distance = 8;
      common_center = {0., 0.};
      set_up_global_variables();
    }
  };

} // namespace testbeam

//  --- testbeam
std::vector<std::array<std::array<float, 2>, 2>> testbeam::measure_resolution_in_mult_photons(TH2F *radii_vs_photons, int min_entries)
{
  std::vector<std::array<std::array<float, 2>, 2>> result;

  //  Radii
  bool found_first_full_slice = false;
  int second_ybin = 1;
  for (auto ybin = 1; ybin <= radii_vs_photons->GetNbinsY(); ybin++)
  {
    //  Get current slice under test
    auto current_slice = radii_vs_photons->ProjectionX(Form("proj_%i", ybin), ybin, ybin);
    if (!found_first_full_slice && current_slice->GetEntries() == 0)
      continue;
    else
      found_first_full_slice = true;

    auto ybin_sec = ybin;
    if (current_slice->GetEntries() < min_entries)
    {
      ybin_sec++;
      for (; ybin_sec <= radii_vs_photons->GetNbinsY(); ybin_sec++)
      {
        auto current_slice = radii_vs_photons->ProjectionX(Form("proj_%i", ybin), ybin, ybin_sec);
        if (current_slice->GetEntries() >= min_entries)
          break;
      }
    }

    //  y-fit
    TF1 *gaus = new TF1("gaus", "gaus");
    current_slice->Fit(gaus);

    // x-determination
    auto y_low_edge = radii_vs_photons->GetYaxis()->GetBinLowEdge(ybin);
    auto y_high_edge = radii_vs_photons->GetYaxis()->GetBinLowEdge(ybin_sec + 1);

    // result
    result.push_back({{{(float)(0.5 * (y_high_edge + y_low_edge)), (float)(0.5 * (y_high_edge - y_low_edge))}, {(float)gaus->GetParameter(2), (float)gaus->GetParError(2)}}});

    ybin = ybin_sec;
  }
  return result;
}

//  --- --- testbeam::data
//  --- --- --- Selection
//  [0] selected points [1] rejected points
std::array<std::vector<int>, 2> testbeam::data::select_points(std::array<float, 4> center, std::array<float, 4> tolerance)
{
  std::array<std::vector<int>, 2> result;
  for (auto ipnt = 0; ipnt < get_n(); ipnt++)
  {
    if (fabs(get_x(ipnt) - center[0]) > tolerance[0])
    {
      result[1].push_back(ipnt);
      continue;
    }
    if (fabs(get_y(ipnt) - center[1]) > tolerance[1])
    {
      result[1].push_back(ipnt);
      continue;
    }
    if (fabs(get_t(ipnt) - center[3]) > tolerance[3])
    {
      result[1].push_back(ipnt);
      continue;
    }
    if (fabs(get_radius(ipnt, {center[0], center[1]}) - center[2]) > tolerance[2])
    {
      result[1].push_back(ipnt);
      continue;
    }
    result[0].push_back(ipnt);
  }
  return result;
}

//  --- --- --- Fit
fit_circle_result testbeam::data::fit_circle(std::vector<int> selected_points, std::array<float, 3> initial_values, bool fix_XY)
{
  auto target_graph = data_graph(selected_points);

  std::array<float, 3> null_initial_values = {0., 0., 0.};
  if (initial_values == null_initial_values)
  {
    initial_values = {common_center[0], common_center[1], get_average_radius(selected_points)};
    fix_XY = true;
  }

  return fit_circle(target_graph, initial_values, fix_XY);
}
fit_circle_result testbeam::data::fit_circle(TGraphErrors *target_graph, std::array<float, 3> initial_values, bool fix_XY)
{
  fit_circle_result result;
  //  Chi2 minimisation for points in a circle
  auto chi2_function = [&](const double *parameters)
  {
    auto chi2 = 0.;
    for (int ipnt = 0; ipnt < target_graph->GetN(); ipnt++)
    {
      auto delta_x = target_graph->GetPointX(ipnt) - parameters[0];
      auto delta_y = target_graph->GetPointY(ipnt) - parameters[1];
      auto delta_r = sqrt(delta_x * delta_x + delta_y * delta_y) - parameters[2];
      chi2 += delta_r * delta_r;
    }
    return chi2;
  };

  // wrap chi2 function in a function object for the fit
  ROOT::Math::Functor fit_function(chi2_function, 3);
  ROOT::Fit::Fitter fitter;

  //  Set initial values and variables names
  double internal_initial_values[3] = {initial_values[0], initial_values[1], initial_values[2]};
  fitter.SetFCN(fit_function, internal_initial_values);
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("R");
  fitter.Config().ParSettings(2).SetLowerLimit(0);
  if (fix_XY)
  {
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
  }

  if (!fitter.FitFCN())
  {
    Error("fit_circle", "Fit failed");
    //  return {{{-2., 0.}, {-2., 0.}, {-2., 0.}}};
  }

  const ROOT::Fit::FitResult &fit_circle_result = fitter.Result();
  auto iTer = -1;
  for (auto current_parameter : fit_circle_result.Parameters())
  {
    iTer++;
    result[iTer][0] = current_parameter;
    result[iTer][1] = fit_circle_result.Errors()[iTer];
  }
  return result;
}
std::vector<float> testbeam::data::diff_circle(std::vector<int> selected_points, std::array<float, 3> initial_values = {0., 0., 0.}, bool fix_XY = false)
{
  std::vector<float> result;
  auto best_fit = fit_circle(selected_points);
  for (auto current_point : selected_points)
  {
    auto point_radius = get_radius(current_point, {best_fit[0][0], best_fit[1][0]});
    result.push_back(point_radius - best_fit[2][0]);
  }
  return result;
}

//  --- --- --- Data I/O
std::pair<TTree *, TFile *> testbeam::data::load_testbeam_data_tree(const std::string &infilename)
{
  auto input_file = TFile::Open(infilename.c_str());
  auto input_tree = (TTree *)input_file->Get("recodata");
  input_tree->SetBranchAddress("n", &this->current_hit.n);
  input_tree->SetBranchAddress("x", &this->current_hit.x);
  input_tree->SetBranchAddress("y", &this->current_hit.y);
  input_tree->SetBranchAddress("t", &this->current_hit.t);
  return {input_tree, input_file};
}

//  --- --- --- DBSCAN
float testbeam::data::get_distance(int i_index, int j_index, std::array<float, 3> distance_weight = {1., 1., 1.})
{
  auto x_distance = (get_x(i_index) - get_x(j_index)) * (distance_weight[0]);
  auto y_distance = (get_y(i_index) - get_y(j_index)) * (distance_weight[1]);
  auto t_distance = (get_t(i_index) - get_t(j_index)) * (distance_weight[2]);
  auto distance = sqrt(x_distance * x_distance + y_distance * y_distance + t_distance * t_distance);
  return distance;
}
float testbeam::data::get_r_distance(int i_index, int j_index, std::array<float, 2> center = {0., 0.})
{
  std::array<float, 2> null_center = {0., 0.};
  if (center == null_center)
    center = common_center;
  auto r_distance = (get_radius(i_index) - get_radius(j_index));
  auto t_distance = (get_t(i_index) - get_t(j_index));
  auto distance = sqrt(r_distance * r_distance + t_distance * t_distance);
  return distance;
}
float testbeam::data::get_radius(int index, std::array<float, 2> center = {0., 0.}) const
{
  std::array<float, 2> null_center = {0., 0.};
  if (center == null_center)
    center = common_center;
  return sqrt((get_x(index) - center[0]) * (get_x(index) - center[0]) + (get_y(index) - center[1]) * (get_y(index) - center[1]));
}
float testbeam::data::get_average_radius(std::vector<int> cluster, std::array<float, 2> center = {0., 0.}) const
{
  std::array<float, 2> null_center = {0., 0.};
  if (center == null_center)
    center = common_center;
  auto result = 0.;
  for (auto ipnt : cluster)
    result += get_radius(ipnt, center);
  return result / (cluster.size());
}
void testbeam::data::dfs_search(int node, std::map<int, std::map<int, float>> &admissible_pairs, std::set<int> &visited, std::vector<int> &cluster)
{
  std::stack<int> stack; // Use a stack for iterative DFS
  stack.push(node);      // Start the DFS from the given node

  while (!stack.empty())
  {
    int current = stack.top(); // Get the current node
    stack.pop();               // Remove it from the stack

    // If the current node hasn't been visited yet
    if (visited.count(current) == 0)
    {
      visited.insert(current);    // Mark the node as visited
      cluster.push_back(current); // Add it to the current cluster

      // Iterate over all neighbors of the current node
      for (const auto &neighbor : admissible_pairs[current])
        if (visited.count(neighbor.first) == 0)
          stack.push(neighbor.first); // Add unvisited neighbors to the stack
    }
  }
}
std::vector<std::vector<int>> testbeam::data::find_clusters(std::map<int, std::map<int, float>> &admissible_pairs)
{
  std::set<int> visited;
  std::vector<std::vector<int>> clusters;
  for (const auto &pair : admissible_pairs)
  {
    int node = pair.first;

    if (visited.count(node) == 0)
    {
      std::vector<int> cluster;
      dfs_search(node, admissible_pairs, visited, cluster);
      clusters.push_back(cluster);
    }
  }
  return clusters;
}
std::vector<std::vector<int>> testbeam::data::get_clusters(float max_distance)
{
  max_distance = max_distance < 0 ? clustering_point_max_distance : max_distance;
  std::map<int, std::map<int, float>> admissible_pairs;
  for (auto ipnt = 0; ipnt < get_n(); ipnt++)
  {
    for (auto jpnt = ipnt + 1; jpnt < get_n(); jpnt++)
    {
      auto current_distance = get_distance(ipnt, jpnt);
      if (current_distance > max_distance)
        continue;
      admissible_pairs[ipnt][jpnt] = current_distance;
    }
  }
  auto result = find_clusters(admissible_pairs);
  std::sort(result.begin(), result.end(), [](const std::vector<int> &a, const std::vector<int> &b)
            { return a.size() > b.size(); });
  return result;
}
std::vector<std::vector<int>> testbeam::data::get_r_clusters(float max_distance)
{
  max_distance = max_distance < 0 ? clustering_r_point_max_distance : max_distance;
  std::map<int, std::map<int, float>> admissible_pairs;
  for (auto ipnt = 0; ipnt < get_n(); ipnt++)
  {
    for (auto jpnt = ipnt + 1; jpnt < get_n(); jpnt++)
    {
      auto current_distance = get_r_distance(ipnt, jpnt);
      if (current_distance > max_distance)
        continue;
      admissible_pairs[ipnt][jpnt] = current_distance;
    }
  }
  auto result = find_clusters(admissible_pairs);
  std::sort(result.begin(), result.end(), [](const std::vector<int> &a, const std::vector<int> &b)
            { return a.size() > b.size(); });
  return result;
}
std::vector<std::vector<int>> testbeam::data::merge_clusters(float max_distance)
{
  max_distance = max_distance < 0 ? clustering_cluster_max_distance : max_distance;
  auto target_clusters = get_clusters();
  std::map<int, std::map<int, float>> admissible_pairs;
  for (auto icluster = 0; icluster < target_clusters.size(); icluster++)
    for (auto jcluster = icluster + 1; jcluster < target_clusters.size(); jcluster++)
    {
      auto current_distance = get_distance(target_clusters[icluster], target_clusters[jcluster]);
      if (current_distance > max_distance)
        continue;
      admissible_pairs[icluster][jcluster] = current_distance;
    }
  return find_clusters(admissible_pairs);
}
std::vector<std::vector<int>> testbeam::data::get_merged_clusters()
{
  std::vector<std::vector<int>> result = {};
  auto target_clusters = get_clusters();
  auto merged_clusters = merge_clusters();
  auto current_cluster_list_iter = -1;
  for (auto current_cluster_list : merged_clusters)
  {
    current_cluster_list_iter++;
    result.push_back({});
    for (auto current_cluster : current_cluster_list)
      for (auto current_point : target_clusters[current_cluster])
        result[current_cluster_list_iter].push_back(current_point);
  }
  return result;
}
std::vector<std::vector<int>> testbeam::data::get_ring_seed()
{
  // std::vector<std::vector<int>> rings;
  // for (auto)
  return {};
}

//  --- --- --- Plotting tools
TCanvas *testbeam::data::get_standard_canvas()
{
  TCanvas *result = new TCanvas("std_canvas", "", 1000, 1000);
  gStyle->SetOptStat(0);
  result->SetMargin(0.13, 0.05, 0.13, 0.05);
  result->DrawFrame(-100., -100., 100., 100., ";x (mm); y (mm)");
  return result;
}
TGraphErrors *testbeam::data::data_graph(int marker_style, int marker_color)
{
  std::vector<int> selected_points;
  for (auto ipnt = 0; ipnt < get_n(); ipnt++)
    selected_points.push_back(ipnt);
  return data_graph(selected_points, marker_style, marker_color);
}
TGraphErrors *testbeam::data::data_graph(std::vector<int> selected_points, int marker_style, int marker_color)
{
  std::vector<std::vector<int>> vector_selected_points = {selected_points};
  return data_graph(vector_selected_points, marker_style, marker_color);
}
TGraphErrors *testbeam::data::data_graph(std::vector<std::vector<int>> selected_points, int marker_style, int marker_color)
{
  TGraphErrors *plot_graph = new TGraphErrors();
  auto i_pnt = -1;
  for (auto cluster_points : selected_points)
    for (auto ipnt : cluster_points)
    {
      i_pnt++;
      plot_graph->SetPoint(i_pnt, get_x(ipnt), get_y(ipnt));
    }
  plot_graph->SetMarkerStyle(marker_style);
  plot_graph->SetMarkerColor(marker_color);
  return plot_graph;
}
TGraphErrors *testbeam::data::add_data_graph(TGraphErrors *plot_graph, int marker_style, int marker_color)
{
  std::vector<int> selected_points;
  for (auto ipnt = 0; ipnt < get_n(); ipnt++)
    selected_points.push_back(ipnt);
  return add_data_graph(plot_graph, selected_points, marker_style, marker_color);
}
TGraphErrors *testbeam::data::add_data_graph(TGraphErrors *plot_graph, std::vector<int> selected_points, int marker_style, int marker_color)
{
  std::vector<std::vector<int>> vector_selected_points = {selected_points};
  return add_data_graph(plot_graph, vector_selected_points, marker_style, marker_color);
}
TGraphErrors *testbeam::data::add_data_graph(TGraphErrors *plot_graph, std::vector<std::vector<int>> selected_points, int marker_style, int marker_color)
{
  for (auto cluster_points : selected_points)
    for (auto ipnt : cluster_points)
    {
      auto i_pnt = plot_graph->GetN();
      plot_graph->SetPoint(i_pnt, get_x(ipnt), get_y(ipnt));
    }
  plot_graph->SetMarkerStyle(marker_style);
  plot_graph->SetMarkerColor(marker_color);
  return plot_graph;
}
TEllipse *testbeam::data::plot_circle(fit_circle_result parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1)
{
  auto result = new TEllipse(parameters[0][0], parameters[1][0], parameters[2][0]);
  result->SetFillStyle(0);
  result->SetLineColor(line_color);
  result->SetLineStyle(line_style);
  result->SetLineWidth(line_width);
  return result;
}

//  --- --- --- Utility
void testbeam::data::set_up_global_variables()
{
  //  iterator for events
  auto iev = 0;

  //  Graph to collect all clusters of points
  TGraphErrors *data_graph_target;
  //  Histogram to set the time cut limits
  TH1F *time_cut = new TH1F("time_cut", ";#Delta_{t_{hit}-t_{ref}}", 200, -100, 100);

  //  First loop to find first good event
  for (; iev < current_tree->GetEntries(); iev++)
  {
    //  Recover event
    current_tree->GetEntry(iev);
    //  Fill time histo & set available SiPM
    for (auto ihit = 0; ihit < get_n(); ihit++)
    {
      time_cut->Fill(get_t(ihit));
      set_available_SiPM({get_x(ihit), get_y(ihit)});
    }
    //  Clustering algorithm
    auto merged_cluster_points = get_merged_clusters();
    //  at least one found cluster
    if (merged_cluster_points.size() != 1)
      continue;
    data_graph_target = data_graph(merged_cluster_points[0]);
    break;
  }
  //  Once the first graph is built, look in all data the good events
  for (; iev < current_tree->GetEntries(); iev++)
  {
    //  Recover event
    current_tree->GetEntry(iev);
    //  Fill time histo
    for (auto ihit = 0; ihit < get_n(); ihit++)
    {
      time_cut->Fill(get_t(ihit));
      set_available_SiPM({get_x(ihit), get_y(ihit)});
    }
    //  Clustering algorithm
    auto merged_cluster_points = get_merged_clusters();
    //  at least one found cluster
    if (merged_cluster_points.size() != 1)
      continue;
    add_data_graph(data_graph_target, merged_cluster_points[0]);
  }

  //  Fit the cumulativeevent clusters
  auto fit_circle_result = fit_circle(data_graph_target);

  //  Fit the time distribution
  time_cut->Fit("gaus", "IMESQ");
  auto gaus_fit = time_cut->GetFunction("gaus");

  //  Set the gloabl variables
  set_common_center({fit_circle_result[0][0], fit_circle_result[1][0]});
  set_common_radii({fit_circle_result[2][0]});
  set_timing_center_sigma({(float)gaus_fit->GetParameter(1), (float)gaus_fit->GetParameter(2)});
}
