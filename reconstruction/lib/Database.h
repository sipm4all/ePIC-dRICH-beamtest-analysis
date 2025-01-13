//  --- --- ---
//  RUN DATABASE UTILITY & CONSTANTS
//  author: nicola.rubini@bo.infn.it
//  --- --- ---

#pragma once

//  Constants
//  --- Conversion constants
const float GeV = 1.;
const float MeV = 1.e3;
const float keV = 1.e6;
const float eV = 1.e9;
//  --- Physical constants
std::map<std::string, float> particle_mass = {
    {"proton", 938.2720882 / MeV},
    {"pion+-", 139.5704 / MeV},
    {"kaon+-", 493.678 / MeV}};
//  --- Sensor
const float kSiPM_x_dimension = 1.5;
const float kSiPM_y_dimension = 1.5;

namespace database
{
  //  run info database
  //  [run_tag][quantity] -> value
  std::map<std::string, std::map<std::string, std::string>> database = {
      {"20231010-084623", {{"temperature", "253.15K"}}},
      {"20231015-150049", {{"polarity", "positive"}}},
      {"20231010-231817", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "135mm"}}},
      {"20231010-235919", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "125mm"}}},
      {"20231011-015609", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "115mm"}}},
      {"20231011-021453", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "105mm"}}},
      {"20231011-023619", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "95mm"}}},
      {"20231011-041853", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "85mm"}}},
      {"20231011-044520", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "75mm"}}},
      {"20231011-051622", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "65mm"}}},
      {"20231011-053156", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "55mm"}}},
      {"20231011-055656", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "45mm"}}},
      {"20231011-061339", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "35mm"}}},
      {"20231011-064602", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "25mm"}}},
      {"20231011-070309", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "15mm"}}},
      {"20231011-071607", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "5mm"}}},
      {"20231010-163636", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "5mm"}}},
      {"20231010-171133", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "15mm"}}},
      {"20231010-173650", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "25mm"}}},
      {"20231010-181846", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "35mm"}}},
      {"20231010-190826", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "45mm"}}},
      {"20231010-193642", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "55mm"}}},
      {"20231010-195910", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "65mm"}}},
      {"20231010-201928", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "75mm"}}},
      {"20231010-204934", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "85mm"}}},
      {"20231010-211228", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "95mm"}}},
      {"20231010-213426", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "105mm"}}},
      {"20231010-221114", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "115mm"}}},
      {"20231010-223412", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "125mm"}}},
      {"20231010-225620", {{"temperature", "253.15K"}, {"V_bias", "52V"}, {"aerogel_mirror", "135mm"}}}};

  //  return standard values if not specified
  std::map<std::string, std::string> standard_values = {
      //  Environment conditions
      {"temperature", "243.15K"},
      //  Beam conditions
      {"polarity", "negative"},
      {"energy", "10GeV"},
      //  Detector conditions
      {"V_bias", "53V"},
      {"aerogel_mirror", "75mm"},
      {"gas_mirror", "0mm"},
      //  Radiator conditions
      {"has_gas", "false"},
      {"has_aerogel", "true"},
      {"n_gas", "1.0008"},
      {"n_aerogel", "1.02"}};

  //  general requests
  //  --- checkers
  inline bool is_run_recorded(std::string run_tag)
  {
    return database.find(run_tag) != database.end();
  }
  inline bool is_info_available(std::string run_tag, std::string quantity)
  {
    if (is_run_recorded(run_tag))
      return database[run_tag].find(quantity) != database[run_tag].end();
    return false;
  }
  //  --- getters
  std::string get_info(std::string run_tag, std::string quantity)
  {
    if (is_run_recorded(run_tag))
    {
      if (is_info_available(run_tag, quantity))
      {
        return database[run_tag][quantity];
      }
      else
      {
        return standard_values[quantity];
      }
    }
    return "Run not found!";
  }
  // --- --- Environment conditions
  template <bool return_celsius = true>
  std::pair<bool, float> get_temperature(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto energy_string = get_info(run_tag, "temperature");
    result.second = return_celsius ? std::stof(energy_string) - 273.15 : std::stof(energy_string);
    return result;
  };
  // --- --- Beam conditions
  std::pair<bool, float> get_beam_energy(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto energy_string = get_info(run_tag, "energy");
    result.second = std::stof(energy_string);
    return result;
  };
  std::pair<bool, int> get_beam_polarity(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto polarity_string = get_info(run_tag, "polarity");
    result.second = 0;
    if ((polarity_string == "positive"))
      result.second = +1;
    if ((polarity_string == "negative"))
      result.second = -1;
    return result;
  };
  //  Detector conditions
  std::pair<bool, float> get_Vbias(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto info_string = get_info(run_tag, "V_bias");
    result.second = std::stof(info_string);
    return result;
  };
  std::pair<bool, float> get_gas_mirror_z(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto info_string = get_info(run_tag, "gas_mirror");
    result.second = std::stof(info_string);
    return result;
  };
  std::pair<bool, float> get_aerogel_mirror_z(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto info_string = get_info(run_tag, "aerogel_mirror");
    result.second = std::stof(info_string);
    return result;
  };
  //  Radiator conditions
  std::pair<bool, bool> get_has_gas(std::string run_tag)
  {
    std::pair<bool, bool> result;
    result.first = is_run_recorded(run_tag);
    auto has_gas_string = get_info(run_tag, "has_gas");
    result.second = !has_gas_string.empty() && (strcasecmp(has_gas_string.c_str(), "true") == 0 || atoi(has_gas_string.c_str()) != 0);
    return result;
  };
  std::pair<bool, bool> get_has_aerogel(std::string run_tag)
  {
    std::pair<bool, bool> result;
    result.first = is_run_recorded(run_tag);
    auto has_gas_string = get_info(run_tag, "has_aerogel");
    result.second = !has_gas_string.empty() && (strcasecmp(has_gas_string.c_str(), "true") == 0 || atoi(has_gas_string.c_str()) != 0);
    return result;
  };
  std::pair<bool, float> get_n_gas(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto energy_string = get_info(run_tag, "n_gas");
    result.second = std::stof(energy_string);
    return result;
  };
  std::pair<bool, float> get_n_aerogel(std::string run_tag)
  {
    std::pair<bool, float> result;
    result.first = is_run_recorded(run_tag);
    auto energy_string = get_info(run_tag, "n_aerogel");
    result.second = std::stof(energy_string);
    return result;
  };
  //  --- setters

  //  Run lists
  std::map<std::string, std::vector<std::string>> run_lists = {
      {"mirror_scan_bkw", {"20231010-231817", "20231010-235919", "20231011-015609", "20231011-021453", "20231011-023619", "20231011-041853", "20231011-044520", "20231011-051622", "20231011-053156", "20231011-055656", "20231011-061339", "20231011-064602", "20231011-070309", "20231011-071607"}},
      {"mirror_scan_fwd", {"20231010-163636", "20231010-171133", "20231010-173650", "20231010-181846", "20231010-190826", "20231010-193642", "20231010-195910", "20231010-201928", "20231010-204934", "20231010-211228", "20231010-213426", "20231010-221114", "20231010-223412", "20231010-225620"}}};

}