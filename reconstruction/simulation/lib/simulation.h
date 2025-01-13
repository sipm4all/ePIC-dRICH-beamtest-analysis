// Includes
#include "../../macros/produce_raw_data_for_radius_cut_calibration.C"

//  Constants
//  --- Physical constants
const float planck_constant = 4.1357 * 1.e-12;   // MeV * ns
const float clight_constant = 2.99792458 * 1.e8; // nm / ns
//  --- --- Conversion constants
//  --- --- --- Energy
const float GeV = 1.e-3;
const float MeV = 1.;
const float keV = 1.e3;
const float eV = 1.e6;
//  --- --- --- Distance
const float m = 1.e-9;
const float mm = 1.e-6;
const float um = 1.e-3;
const float nm = 1.;
//  --- --- --- Time
const float s = 1.e-9;
const float ms = 1.e-6;
const float us = 1.e-3;
const float ns = 1.;
//  --- Simulation constants
//  --- --- Detector
const float DETECTOR_AVG_POSITION_Z = -16.;
const float DETECTOR_TOLERANCE_Z = 1.;
//  --- --- Aerogel
const float AEROGEL_MID_POINT_Z = -100.;
const float AEROGEL_TOLERANCE_Z = 25.;
const float AEROGEL_MIRROR_RADIUS = 700.;
const float AEROGEL_MIRROR_CENTER = -380.;
const float AEROGEL_MIRROR_START_POSITION_Z = 317.;
const float AEROGEL_MIRROR_END_POSITION_Z = 322.;
//  --- --- Gas
// const float GAS_END_POSITION_Z = 61.;
const float GAS_MIRROR_RADIUS = 2420.;
const float GAS_MIRROR_CENTER = -380.;
const float GAS_MIRROR_START_POSITION_Z = 1189.;
const float GAS_MIRROR_END_POSITION_Z = 1191.;
//  --- PDE approximation
std::map<std::string, std::array<float, 3>> PDE_approx_params = {
    {"HPK S13360-3050VS", {229.072, 466.963, 70.2685}}};

struct simdata
{
  //  Event hits
  int nhits;
  //  particle ID
  int id[16384];
  //  particle energy
  double E[16384];
  //  particle production vertex
  double vtx_x[16384], vtx_y[16384], vtx_z[16384];
  //  particle momentum
  double px[16384], py[16384], pz[16384];
  //  particle position
  double hit_x[16384], hit_y[16384], hit_z[16384], hit_t[16384];
};

void read_tree(TTree *input_tree, simdata &read_data_struct)
{
  input_tree->SetBranchAddress("nhits", &read_data_struct.nhits);
  input_tree->SetBranchAddress("trackE", &read_data_struct.E);
  input_tree->SetBranchAddress("avg_x", &read_data_struct.hit_x);
  input_tree->SetBranchAddress("avg_y", &read_data_struct.hit_y);
  input_tree->SetBranchAddress("avg_z", &read_data_struct.hit_z);
  input_tree->SetBranchAddress("avg_t", &read_data_struct.hit_t);
  input_tree->SetBranchAddress("px", &read_data_struct.px);
  input_tree->SetBranchAddress("py", &read_data_struct.py);
  input_tree->SetBranchAddress("pz", &read_data_struct.pz);
  input_tree->SetBranchAddress("vx", &read_data_struct.vtx_x);
  input_tree->SetBranchAddress("vy", &read_data_struct.vtx_y);
  input_tree->SetBranchAddress("vz", &read_data_struct.vtx_z);
  input_tree->SetBranchAddress("id", &read_data_struct.id);
}

std::pair<TTree *, TFile *> read_tree(std::string input_filename, simdata &read_data_struct)
{
  TFile *input_file = new TFile(input_filename.c_str());
  TTree *input_tree = (TTree *)(input_file->Get("htree"));
  read_tree(input_tree, read_data_struct);
  return {input_tree, input_file};
}

bool is_hit_vertex_in_aerogel(simdata data, int iTer)
{
  if (fabs(data.vtx_z[iTer] - AEROGEL_MID_POINT_Z) < AEROGEL_TOLERANCE_Z)
    return true;
  return false;
}

bool is_photon(simdata data, int iTer)
{
  return data.id[iTer] == 0;
}

bool is_hit_on_detection_plane(simdata data, int iTer)
{
  if (fabs(data.hit_z[iTer] - DETECTOR_AVG_POSITION_Z) < DETECTOR_TOLERANCE_Z)
    return true;
  return false;
}

bool is_aerogel_photon_on_detection_plane(simdata data, int iTer)
{
  return ((is_hit_vertex_in_aerogel(data, iTer) && is_photon(data, iTer)) && is_hit_on_detection_plane(data, iTer));
}

bool is_gas_photon_on_detection_plane(simdata data, int iTer)
{
  return ((!is_hit_vertex_in_aerogel(data, iTer) && is_photon(data, iTer)) && is_hit_on_detection_plane(data, iTer));
}

float get_wavelenght(float photon_energy) { return planck_constant * clight_constant / (photon_energy); }

float get_approximate_PDE(std::string target_sensor, float wavelength)
{
  return PDE_approx_params[target_sensor][0]*TMath::Landau(wavelength,PDE_approx_params[target_sensor][1],PDE_approx_params[target_sensor][2])/100.;
}
