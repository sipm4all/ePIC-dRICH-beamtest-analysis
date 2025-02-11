#pragma once

namespace alcor
{
  // --- Constants

  // --- Conversion
  const double coarse_to_s = 3.1250000e-09;
  const double coarse_to_ms = 3.1250000e-06;
  const double coarse_to_us = 3.1250000e-03;
  const double coarse_to_ns = 3.1250000;

  const int rollover_to_coarse = 32768;
  const double rollover_to_s = 0.0001024;
  const double rollover_to_ms = 0.1024;
  const double rollover_to_us = 102.4;
  const double rollover_to_ns = 102400.0;

  // --- Data structure
  struct hit
  {
    int device, fifo, type, counter, column, pixel, tdc, rollover, coarse, fine;
    double time = -1;
  };

  enum hit_type
  {
    alcor_hit = 1,
    trigger_tag = 9,
    start_spill = 7,
    end_spill = 15
  };

  class data
  {
  private:
    // Variables
    alcor::hit current_hit;
    //  [0] IF, OFF, CUT [1] IF, OFF
    char16_t alcor_mode_flag;
    TTree *current_tree = nullptr;
    TFile *current_file = nullptr;

  public:
    // --- Constructor
    data() {}
    data(const std::string &infilename)
    {
      auto current_tree_and_file = load_alcor_data_tree(infilename);
      current_tree = std::get<0>(current_tree_and_file);
      current_file = std::get<1>(current_tree_and_file);
    }

    // --- Getters & Setters

    // --- Getters
    int get_device() const { return current_hit.device; }
    int get_fifo() const { return current_hit.fifo; }
    int get_type() const { return current_hit.type; }
    int get_counter() const { return current_hit.counter; }
    int get_column() const { return current_hit.column; }
    int get_pixel() const { return current_hit.pixel; }
    int get_tdc() const { return current_hit.tdc; }
    int get_rollover() const { return current_hit.rollover; }
    int get_coarse() const { return current_hit.coarse; }
    int get_fine() const { return current_hit.fine; }
    double get_time() const { return current_hit.time; }
    TTree *get_tree() const { return current_tree; }
    TFile *get_file() const { return current_file; }
    char16_t get_alcor_mode_flag() const { return alcor_mode_flag; }

    // --- Derived Getters
    int get_chip() const { return get_fifo() / 4; }
    int get_eo_channel() const { return get_pixel() + 4 * get_column(); }
    int get_calibration_index() const { return get_tdc() + 4 * get_eo_channel() + 128 * get_chip(); }
    int get_global_coarse() const { return get_coarse() + get_rollover() * alcor::rollover_to_coarse; }
    double get_coarse_time_ns() const { return get_global_coarse() * alcor::coarse_to_ns; }

    // --- Setters
    void set_device(int value) { current_hit.device = value; }
    void set_fifo(int value) { current_hit.fifo = value; }
    void set_type(int value) { current_hit.type = value; }
    void set_counter(int value) { current_hit.counter = value; }
    void set_column(int value) { current_hit.column = value; }
    void set_pixel(int value) { current_hit.pixel = value; }
    void set_tdc(int value) { current_hit.tdc = value; }
    void set_rollover(int value) { current_hit.rollover = value; }
    void set_coarse(int value) { current_hit.coarse = value; }
    void set_fine(int value) { current_hit.fine = value; }
    void set_time(double value) { current_hit.time = value; }
    void set_tree(TTree *tree) { current_tree = tree; }
    void set_file(TFile *file) { current_file = file; }
    void get_alcor_mode_flag(char16_t value) { alcor_mode_flag = value; }

    // --- Data flags
    bool is_alcor_hit() const { return current_hit.type == alcor_hit; }
    bool is_trigger_tag() const { return current_hit.type == trigger_tag; }
    bool is_start_spill() const { return current_hit.type == start_spill; }
    bool is_end_spill() const { return current_hit.type == end_spill; }

    // --- Data I/O
    std::pair<TTree *, TFile *> load_alcor_data_tree(const std::string &infilename);

    // --- Operators
    bool operator<(const data &rhs) const { return current_hit.time < rhs.get_time(); }
    bool operator>(const data &rhs) const { return current_hit.time > rhs.get_time(); }
    bool operator<=(const data &rhs) const { return current_hit.time <= rhs.get_time(); }
    bool operator>=(const data &rhs) const { return current_hit.time >= rhs.get_time(); }
  };

  using calibration_unit = std::array<float, 3>; // [0] OFF, [1] IF, [2] CUT

  class fine_data : public data
  {
  private:
    // --- Variables ---
    // Calibration scheme:
    // calibration_scheme[device][cindex] stores:
    // [0] OFF / Intercept
    // [1] IF / Slope
    // [2] CUT / Unused
    const std::vector<std::string> calibration_csv_fields = {"device", "cindex", "OFF", "IF", "CUT"};
    std::map<std::string, std::map<int, calibration_unit>> calibration_database;
    enum ALCOR_v
    {
      tb_2022,
      tb_2023,
      tb_2024
    };

  public:
    // --- Constructors ---
    fine_data() {}

    // --- Getters ---
    float get_OFF(const std::string &device, int cindex) const { return calibration_database.at(device).at(cindex)[0]; }
    float get_IF(const std::string &device, int cindex) const { return calibration_database.at(device).at(cindex)[1]; }
    float get_CUT(const std::string &device, int cindex) const { return calibration_database.at(device).at(cindex)[2]; }

    // --- Setters ---
    void set_OFF(const std::string &device, int cindex, float value) { calibration_database[device][cindex][0] = value; }
    void set_IF(const std::string &device, int cindex, float value) { calibration_database[device][cindex][1] = value; }
    void set_CUT(const std::string &device, int cindex, float value) { calibration_database[device][cindex][2] = value; }

    // --- Calculate phase ---
    float get_phase(const std::string &device, int cindex, float fine) const;
    float get_phase() const { return get_phase(std::to_string(get_device()), get_calibration_index(), get_fine()); }

    // --- I/O calibration scheme from file ---
    void read_calibration(const std::string &filename);
    void write_calibration(const std::string &filename);
  };

  struct vec_hit
  {
    unsigned char device[1024], index[1024], tdc[1024], coarse[1024], fine[1024];
    double time[1024] = {-1};
  };

  class trkdata : public fine_data
  {
  private:
    alcor::vec_hit current_vec_hit;
    TTree *current_tree = nullptr;
    TFile *current_file = nullptr;
    unsigned short current_spill;
    unsigned int current_frame;
    unsigned short current_in_frame_hits;
    unsigned short current_trgmask;
    float current_tref;

  public:
    // --- Constructor
    trkdata(const std::string &infilename)
    {
      auto current_tree_and_file = load_alcor_data_tree(infilename);
      current_tree = std::get<0>(current_tree_and_file);
      current_file = std::get<1>(current_tree_and_file);
    }

    // Getters
    int get_channel(int index) const { return (current_vec_hit.index[index]) % 32; }
    int get_chip(int index) const { return (current_vec_hit.index[index] / (32)) % 6; }
    unsigned short get_spill() const { return current_spill; }
    unsigned int get_frame() const { return current_frame; }
    unsigned short get_in_frame_hits() const { return current_in_frame_hits; }
    unsigned short get_trgmask() const { return current_trgmask; }
    float get_tref() const { return current_tref; }
    unsigned char get_device(int index) const { return current_vec_hit.device[index]; }
    unsigned char get_tdc(int index) const { return current_vec_hit.tdc[index]; }
    unsigned char get_coarse(int index) const { return current_vec_hit.coarse[index]; }
    unsigned char get_fine(int index) const { return current_vec_hit.fine[index]; }
    unsigned char get_index(int index) const { return current_vec_hit.index[index]; }
    double get_time(int index) const { return current_vec_hit.time[index]; }
    alcor::vec_hit &get_current_vec_hit() { return current_vec_hit; }
    TTree *get_tree() const { return current_tree; }
    TFile *get_file() const { return current_file; }

    // Setters
    void set_spill(unsigned short value) { current_spill = value; }
    void set_frame(unsigned int value) { current_frame = value; }
    void set_in_frame_hits(unsigned short value) { current_in_frame_hits = value; }
    void set_trgmask(unsigned short value) { current_trgmask = value; }
    void set_tref(float value) { current_tref = value; }
    void set_device(int index, unsigned char value) { current_vec_hit.device[index] = value; }
    void set_tdc(int index, unsigned char value) { current_vec_hit.tdc[index] = value; }
    void set_coarse(int index, unsigned char value) { current_vec_hit.coarse[index] = value; }
    void set_fine(int index, unsigned char value) { current_vec_hit.fine[index] = value; }
    void set_index(int index, unsigned char value) { current_vec_hit.index[index] = value; }
    void set_time(int index, double value) { current_vec_hit.time[index] = value; }
    void set_current_vec_hit(const alcor::vec_hit &vec) { current_vec_hit = vec; }
    void set_tree(TTree *tree) { current_tree = tree; }
    void set_file(TFile *file) { current_file = file; }

    // --- Data I/O
    std::pair<TTree *, TFile *> load_alcor_data_tree(const std::string &infilename);
  };

} // namespace alcor

std::pair<TTree *, TFile *>
alcor::data::load_alcor_data_tree(const std::string &infilename)
{
  auto input_file = TFile::Open(infilename.c_str());
  auto input_tree = (TTree *)input_file->Get("alcor");
  input_tree->SetBranchAddress("fifo", &this->current_hit.fifo);
  input_tree->SetBranchAddress("type", &this->current_hit.type);
  input_tree->SetBranchAddress("counter", &this->current_hit.counter);
  input_tree->SetBranchAddress("column", &this->current_hit.column);
  input_tree->SetBranchAddress("pixel", &this->current_hit.pixel);
  input_tree->SetBranchAddress("tdc", &this->current_hit.tdc);
  input_tree->SetBranchAddress("rollover", &this->current_hit.rollover);
  input_tree->SetBranchAddress("coarse", &this->current_hit.coarse);
  input_tree->SetBranchAddress("fine", &this->current_hit.fine);
  return {input_tree, input_file};
}

float alcor::fine_data::get_phase(const std::string &device, int cindex, float fine) const
{
  return get_fine() * get_IF(device, cindex) + get_OFF(device, cindex);
}

void alcor::fine_data::read_calibration(const std::string &filename)
{
  std::ifstream data_stream(filename);
  std::string current_line;
  while (std::getline(data_stream, current_line))
  {
    //  Skip comment characters
    if (current_line[0] == '#' || current_line[0] == ' ')
      continue;
    //  Read database
    std::stringstream string_in_stream(current_line);
    std::map<std::string, std::string> data_by_field;
    std::string current_data;
    for (auto current_field : alcor::fine_data::calibration_csv_fields)
    {
      string_in_stream >> current_data;
      data_by_field[current_field] = current_data;
    }
    this->set_OFF(data_by_field["device"], std::stoi(data_by_field["cindex"]), std::stof(data_by_field["OFF"]));
    this->set_IF(data_by_field["device"], std::stoi(data_by_field["cindex"]), std::stof(data_by_field["IF"]));
    this->set_CUT(data_by_field["device"], std::stoi(data_by_field["cindex"]), std::stof(data_by_field["CUT"]));
  }
};

void alcor::fine_data::write_calibration(const std::string &filename)
{
  std::ofstream data_stream(filename);
  data_stream << "#\tCalibration file for Fine Tune\n";
  auto current_time = std::chrono::system_clock::now();
  std::time_t current_time_stamp = std::chrono::system_clock::to_time_t(current_time);
  data_stream << "#\tGenerated on " << std::ctime(&current_time_stamp);
  data_stream << "#\t## Device ## Calibration index ## Intercept ## Slope ##\n";
  for (auto [device, calibration_of_device] : calibration_database)
  {
    for (auto [cindex, current_calibration_unit] : calibration_of_device)
    {
      data_stream << device << "\t" << cindex << "\t" << current_calibration_unit[0] << "\t" << current_calibration_unit[1] << "\n";
    }
  }
};

std::pair<TTree *, TFile *>
alcor::trkdata::load_alcor_data_tree(const std::string &infilename)
{
  auto input_file = TFile::Open(infilename.c_str());
  auto input_tree = (TTree *)input_file->Get("recodata");
  input_tree->SetBranchAddress("spill", &this->current_spill);
  input_tree->SetBranchAddress("frame", &this->current_frame);
  input_tree->SetBranchAddress("tref", &this->current_tref);
  input_tree->SetBranchAddress("trgmask", &this->current_trgmask);
  input_tree->SetBranchAddress("n", &this->current_in_frame_hits);
  input_tree->SetBranchAddress("device", &this->current_vec_hit.device);
  input_tree->SetBranchAddress("index", &this->current_vec_hit.index);
  input_tree->SetBranchAddress("coarse", &this->current_vec_hit.coarse);
  input_tree->SetBranchAddress("fine", &this->current_vec_hit.fine);
  input_tree->SetBranchAddress("tdc", &this->current_vec_hit.tdc);
  return {input_tree, input_file};
}

/*
  //  const int data::eo2do[32] = {22, 20, 18, 16, 24, 26, 28, 30, 25, 27, 29, 31, 23, 21, 19, 17, 9, 11, 13, 15, 7, 5, 3, 1, 6, 4, 2, 0, 8, 10, 12, 14};

  const int data::rollover_to_clock = 32768;
  const double data::coarse_to_ns = 3.125;
  const double data::rollover_to_ns = 102400.;

  double data::fine_min[768] = {0.};
  double data::fine_max[768] = {0.};
  double data::fine_off[768] = {0.};

  void data::link_to_tree(TTree *t)
  {
    if (!t)
      return;
    t->SetBranchAddress("device", &device);
    t->SetBranchAddress("fifo", &fifo);
    t->SetBranchAddress("type", &type);
    t->SetBranchAddress("counter", &counter);
    t->SetBranchAddress("column", &column);
    t->SetBranchAddress("pixel", &pixel);
    t->SetBranchAddress("tdc", &tdc);
    t->SetBranchAddress("rollover", &rollover);
    t->SetBranchAddress("coarse", &coarse);
    t->SetBranchAddress("fine", &fine);
  }

  bool data::load_fine_calibration(std::string filename)
  {
    std::cout << " --- loading fine calibration: " << filename << std::endl;
    auto fin = TFile::Open(filename.c_str());
    auto hFine_min = (TH1 *)fin->Get("hFine_min");
    auto hFine_max = (TH1 *)fin->Get("hFine_max");
    auto hFine_off = (TH1 *)fin->Get("hFine_off");
    int found = 0;
    for (int i = 0; i < 768; ++i)
    {
      if (hFine_min->GetBinError(i + 1) <= 0. ||
          hFine_max->GetBinError(i + 1) <= 0.)
        coninput_treeue;
      fine_min[i] = hFine_min ? hFine_min->GetBinContent(i + 1) : 0.;
      fine_max[i] = hFine_max ? hFine_max->GetBinContent(i + 1) : 0.;
      fine_off[i] = hFine_off ? hFine_off->GetBinContent(i + 1) : 0.;
      found++;
    }
    std::cout << " --- loaded fine calibration: found " << found << " channels " << std::endl;
    return true;
  }

  //
  //  [0] Slope,  [1] Intercept
  typedef std::array<float, 3> calibration_unit;

  //  Class for the fine calibration
  class finecalibration
  {
  public:
    //  --- Variables ---
    //  calibration_scheme[device][cindex]  [0] OFF / Intercept
    //                                      [1] IF  / Slope
    //                                      [2] CUT
    const std::vector<std::string> calibration_csv_fields = {"device", "cindex", "OFF", "IF", "CUT"};
    std::map<std::string, std::map<int, calibration_unit>> calibration_database;

    //  --- Generators ---
    finecalibration() {};
    finecalibration(std::string filename) { read_calibration(filename); };

    //  --- Methods ---
    //  --- Getters
    float get_OFF(std::string device, int cindex) { return calibration_database[device][cindex][0]; }
    float get_IF(std::string device, int cindex) { return calibration_database[device][cindex][1]; }
    float get_CUT(std::string device, int cindex) { return calibration_database[device][cindex][2]; }
    //  --- Setters
    void set_intercept(std::string device, int cindex, float intercept) { calibration_database[device][cindex][3] = intercept; }
    void set_slope(std::string device, int cindex, float slope) { calibration_database[device][cindex][4] = slope; }
    //  --- Calculate phase
    float get_phase(std::string device, int cindex, float fine) { return fine * get_slope(device, cindex) + get_intercept(device, cindex) }
    //  --- I/O calibration scheme from file
    void read_calibration(std::string filename);
    void write_calibration(std::string filename);

  private:
  };
}





*/