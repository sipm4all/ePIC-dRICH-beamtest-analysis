#pragma once

#include "lightdata.h"

namespace sipm4eic
{

  class lightio
  {

  public:
    static const int frame_size = 256;
    static const int max_devices = 256;
    static const int max_frames = 65534;   // maximum number of frames in a spill
    static const int max_triggers = 65534; // maximum number of triggers in a spill
    static const int max_hits = 262144;    // maximum number of hits in a spill
    //
    unsigned char part_n;
    unsigned char part_device[max_devices];
    unsigned int part_mask[max_devices];
    //
    unsigned char dead_n;
    unsigned char dead_device[max_devices];
    unsigned int dead_mask[max_devices];
    //
    unsigned short frame_n;         // number of frames
    unsigned int frame[max_frames]; // frame index
    //
    unsigned int trigger0_size;           // trigger hits in spill
    unsigned char trigger0_n[max_frames]; // trigger hits in frame
    unsigned char trigger0_coarse[max_triggers];
    //
    unsigned int timing_size;            // timing hits in spill
    unsigned short timing_n[max_frames]; // timing hits in frame
    unsigned char timing_device[max_hits];
    unsigned char timing_index[max_hits];
    unsigned char timing_coarse[max_hits];
    unsigned char timing_fine[max_hits];
    unsigned char timing_tdc[max_hits];
    //
    unsigned int cherenkov_size;            // cherenkov hits in spill
    unsigned short cherenkov_n[max_frames]; // cherenkov hits in frame
    unsigned char cherenkov_device[max_hits];
    unsigned char cherenkov_index[max_hits];
    unsigned char cherenkov_coarse[max_hits];
    unsigned char cherenkov_fine[max_hits];
    unsigned char cherenkov_tdc[max_hits];
    //
    std::map<int, std::map<int, int>> active_device_fifos; // Decode map indicating active fifos
    std::map<int, std::map<int, int>> dead_device_fifos;   // Decode map indicating active fifos
    //
    //  --- Constructor
    lightio() = default;
    //
    //  --- Getters
    // --- Getters for Static Constants ---
    static int get_frame_size() { return frame_size; }
    static int get_max_devices() { return max_devices; }
    static int get_max_frames() { return max_frames; }
    static int get_max_triggers() { return max_triggers; }
    static int get_max_hits() { return max_hits; }
    // --- Getters for Part Variables ---
    unsigned char get_part_n() const { return part_n; }
    unsigned char get_part_device(int index) const { return part_device[index]; }
    unsigned int get_part_mask(int index) const { return part_mask[index]; }
    // --- Getters for Dead Variables ---
    unsigned char get_dead_n() const { return dead_n; }
    unsigned char get_dead_device(int index) const { return dead_device[index]; }
    unsigned int get_dead_mask(int index) const { return dead_mask[index]; }
    // --- Getters for Frame Variables ---
    unsigned short get_frame_n() const { return frame_n; }
    unsigned int get_frame(int index) const { return frame[index]; }
    // --- Getters for Trigger0 Variables ---
    unsigned int get_trigger0_size() const { return trigger0_size; }
    unsigned char get_trigger0_n(int index) const { return trigger0_n[index]; }
    unsigned char get_trigger0_coarse(int index) const { return trigger0_coarse[index]; }
    // --- Getters for Timing Variables ---
    unsigned int get_timing_size() const { return timing_size; }
    unsigned short get_timing_n(int index) const { return timing_n[index]; }
    unsigned char get_timing_device(int index) const { return timing_device[index]; }
    unsigned char get_timing_index(int index) const { return timing_index[index]; }
    unsigned char get_timing_coarse(int index) const { return timing_coarse[index]; }
    unsigned char get_timing_fine(int index) const { return timing_fine[index]; }
    unsigned char get_timing_tdc(int index) const { return timing_tdc[index]; }
    // --- Getters for Cherenkov Variables ---
    unsigned int get_cherenkov_size() const { return cherenkov_size; }
    unsigned short get_cherenkov_n(int index) const { return cherenkov_n[index]; }
    unsigned char get_cherenkov_device(int index) const { return cherenkov_device[index]; }
    unsigned char get_cherenkov_index(int index) const { return cherenkov_index[index]; }
    unsigned char get_cherenkov_coarse(int index) const { return cherenkov_coarse[index]; }
    unsigned char get_cherenkov_fine(int index) const { return cherenkov_fine[index]; }
    unsigned char get_cherenkov_tdc(int index) const { return cherenkov_tdc[index]; }
    //  --- Getters for active and dead masks
    std::map<int, std::map<int, int>> get_active_fifos() { return active_device_fifos; };
    std::map<int, int> get_active_fifos(int device) { return active_device_fifos[device]; };
    std::map<int, std::map<int, int>> get_dead_fifos() { return dead_device_fifos; };
    std::map<int, int> get_dead_fifos(int device) { return dead_device_fifos[device]; };
    //
    //  --- Setters
    // --- Setters for Part Variables ---
    void set_part_n(unsigned char value) { part_n = value; }
    void set_part_device(int index, unsigned char value) { part_device[index] = value; }
    void set_part_mask(int index, unsigned int value) { part_mask[index] = value; }
    // --- Setters for Dead Variables ---
    void set_dead_n(unsigned char value) { dead_n = value; }
    void set_dead_device(int index, unsigned char value) { dead_device[index] = value; }
    void set_dead_mask(int index, unsigned int value) { dead_mask[index] = value; }
    // --- Setters for Frame Variables ---
    void set_frame_n(unsigned short value) { frame_n = value; }
    void set_frame(int index, unsigned int value) { frame[index] = value; }
    // --- Setters for Trigger0 Variables ---
    void set_trigger0_size(unsigned int value) { trigger0_size = value; }
    void set_trigger0_n(int index, unsigned char value) { trigger0_n[index] = value; }
    void set_trigger0_coarse(int index, unsigned char value) { trigger0_coarse[index] = value; }
    // --- Setters for Timing Variables ---
    void set_timing_size(unsigned int value) { timing_size = value; }
    void set_timing_n(int index, unsigned short value) { timing_n[index] = value; }
    void set_timing_device(int index, unsigned char value) { timing_device[index] = value; }
    void set_timing_index(int index, unsigned char value) { timing_index[index] = value; }
    void set_timing_coarse(int index, unsigned char value) { timing_coarse[index] = value; }
    void set_timing_fine(int index, unsigned char value) { timing_fine[index] = value; }
    void set_timing_tdc(int index, unsigned char value) { timing_tdc[index] = value; }
    // --- Setters for Cherenkov Variables ---
    void set_cherenkov_size(unsigned int value) { cherenkov_size = value; }
    void set_cherenkov_n(int index, unsigned short value) { cherenkov_n[index] = value; }
    void set_cherenkov_device(int index, unsigned char value) { cherenkov_device[index] = value; }
    void set_cherenkov_index(int index, unsigned char value) { cherenkov_index[index] = value; }
    void set_cherenkov_coarse(int index, unsigned char value) { cherenkov_coarse[index] = value; }
    void set_cherenkov_fine(int index, unsigned char value) { cherenkov_fine[index] = value; }
    void set_cherenkov_tdc(int index, unsigned char value) { cherenkov_tdc[index] = value; }

    void new_spill(unsigned int ispill);
    void new_frame(unsigned int iframe);
    void add_part(unsigned char device, unsigned int mask);
    void add_dead(unsigned char device, unsigned int mask);
    void add_trigger0(unsigned char coarse);
    void add_timing(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc);
    void add_cherenkov(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc);
    void add_frame() { ++frame_n; };
    void fill();
    void write_and_close();
    void write_to_tree(std::string filename, std::string treename = "lightdata");
    void write_to_tree(TTree *t);

    void read_from_tree(std::string filename, std::string treename = "lightdata");
    void read_from_tree(TTree *t);
    bool next_spill();
    bool next_frame();
    void reset() { spill_current = frame_current = 0; };
    void decode_part_device_fifos();
    void decode_dead_device_fifos();
    void decode_device_fifos()
    {
      decode_part_device_fifos();
      decode_dead_device_fifos();
    };

    std::vector<lightdata> &get_trigger0_vector() { return trigger0_vector; };
    std::vector<lightdata> &get_timing_vector() { return timing_vector; };
    std::vector<lightdata> &get_cherenkov_vector() { return cherenkov_vector; };

    std::map<std::array<unsigned char, 2>, std::vector<lightdata>> &get_timing_map() { return timing_map; };
    std::map<std::array<unsigned char, 2>, std::vector<lightdata>> &get_cherenkov_map() { return cherenkov_map; };

    TTree *get_tree() { return tree; };

  private:
    TFile *file = nullptr;
    TTree *tree = nullptr;

    int spill_current = 0;
    int frame_current = 0;

    int trigger0_offset = 0;
    int timing_offset = 0;
    int cherenkov_offset = 0;

    std::vector<lightdata> trigger0_vector;
    std::vector<lightdata> timing_vector;
    std::vector<lightdata> cherenkov_vector;

    std::map<std::array<unsigned char, 2>, std::vector<lightdata>> timing_map;
    std::map<std::array<unsigned char, 2>, std::vector<lightdata>> cherenkov_map;
  };

  void
  lightio::new_spill(unsigned int ispill)
  {
    std::cout << " --- new spill: " << ispill << std::endl;
    part_n = 0;
    dead_n = 0;
    frame_n = 0;
    trigger0_size = 0;
    timing_size = 0;
    cherenkov_size = 0;
  };

  void
  lightio::new_frame(unsigned int iframe)
  {
    frame[frame_n] = iframe;
    trigger0_n[frame_n] = 0;
    timing_n[frame_n] = 0;
    cherenkov_n[frame_n] = 0;
  }

  void
  lightio::add_part(unsigned char device, unsigned int mask)
  {
    part_device[part_n] = device;
    part_mask[part_n] = mask;
    ++part_n;
  }

  void
  lightio::add_dead(unsigned char device, unsigned int mask)
  {
    dead_device[dead_n] = device;
    dead_mask[dead_n] = mask;
    ++dead_n;
  }

  void
  lightio::add_trigger0(unsigned char coarse)
  {
    trigger0_coarse[trigger0_size] = coarse;
    ++trigger0_n[frame_n];
    ++trigger0_size;
  }

  void
  lightio::add_timing(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc)
  {
    timing_device[timing_size] = device;
    timing_index[timing_size] = index;
    timing_coarse[timing_size] = coarse;
    timing_fine[timing_size] = fine;
    timing_tdc[timing_size] = tdc;
    ++timing_n[frame_n];
    ++timing_size;
  }

  void
  lightio::add_cherenkov(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc)
  {
    cherenkov_device[cherenkov_size] = device;
    cherenkov_index[cherenkov_size] = index;
    cherenkov_coarse[cherenkov_size] = coarse;
    cherenkov_fine[cherenkov_size] = fine;
    cherenkov_tdc[cherenkov_size] = tdc;
    ++cherenkov_n[frame_n];
    ++cherenkov_size;
  }

  void
  lightio::fill()
  {
    std::cout << " --- fill tree: trigger0_size = " << trigger0_size << std::endl;
    std::cout << "                  timing_size = " << timing_size << std::endl;
    std::cout << "                     cherenkov_size = " << cherenkov_size << std::endl;
    tree->Fill();
  };

  void
  lightio::write_and_close()
  {
#if 0
  auto n_spills = tree->GetEntries();
  auto n_frames = 0;
  for (int ispill = 0; ispill < n_spills; ++ispill) {
    tree->GetEvent(ispill);
    n_frames += frame_n;
  }
  std::cout << " --- write and close " << std::endl;
  std::cout << " --- collected " << n_spills << " spills " << std::endl;
  std::cout << "               " << n_frames << " frames " << std::endl;
#endif
    file->cd();
    tree->Write();
    file->Close();
  }

  void
  lightio::write_to_tree(std::string filename, std::string treename)
  {
    file = TFile::Open(filename.c_str(), "RECREATE");
    tree = new TTree(treename.c_str(), treename.c_str());
    write_to_tree(tree);
  }

  void
  lightio::write_to_tree(TTree *t)
  {
    t->Branch("part_n", &part_n, "part_n/b");
    t->Branch("part_device", &part_device, "part_device[part_n]/b");
    t->Branch("part_mask", &part_mask, "part_mask[part_n]/i");
    t->Branch("dead_n", &dead_n, "dead_n/b");
    t->Branch("dead_device", &dead_device, "dead_device[dead_n]/b");
    t->Branch("dead_mask", &dead_mask, "dead_mask[dead_n]/i");
    t->Branch("frame_n", &frame_n, "frame_n/s");
    t->Branch("frame", &frame, "frame[frame_n]/i");
    t->Branch("trigger0_size", &trigger0_size, "trigger0_size/i");
    t->Branch("trigger0_n", &trigger0_n, "trigger0_n[frame_n]/b");
    t->Branch("trigger0_coarse", &trigger0_coarse, "trigger0_coarse[trigger0_size]/b");
    t->Branch("timing_size", &timing_size, "timing_size/i");
    t->Branch("timing_n", &timing_n, "timing_n[frame_n]/s");
    t->Branch("timing_device", &timing_device, "timing_device[timing_size]/b");
    t->Branch("timing_index", &timing_index, "timing_index[timing_size]/b");
    t->Branch("timing_coarse", &timing_coarse, "timing_coarse[timing_size]/b");
    t->Branch("timing_fine", &timing_fine, "timing_fine[timing_size]/b");
    t->Branch("timing_tdc", &timing_tdc, "timing_tdc[timing_size]/b");
    t->Branch("cherenkov_size", &cherenkov_size, "cherenkov_size/i");
    t->Branch("cherenkov_n", &cherenkov_n, "cherenkov_n[frame_n]/s");
    t->Branch("cherenkov_device", &cherenkov_device, "cherenkov_device[cherenkov_size]/b");
    t->Branch("cherenkov_index", &cherenkov_index, "cherenkov_index[cherenkov_size]/b");
    t->Branch("cherenkov_coarse", &cherenkov_coarse, "cherenkov_coarse[cherenkov_size]/b");
    t->Branch("cherenkov_fine", &cherenkov_fine, "cherenkov_fine[cherenkov_size]/b");
    t->Branch("cherenkov_tdc", &cherenkov_tdc, "cherenkov_tdc[cherenkov_size]/b");
  }

  void
  lightio::read_from_tree(std::string filename, std::string treename)
  {
    file = TFile::Open(filename.c_str());
    tree = (TTree *)file->Get(treename.c_str());
    read_from_tree(tree);
  }

  void
  lightio::read_from_tree(TTree *t)
  {
    t->SetBranchAddress("part_n", &part_n);
    t->SetBranchAddress("part_device", &part_device);
    t->SetBranchAddress("part_mask", &part_mask);
    t->SetBranchAddress("dead_n", &dead_n);
    t->SetBranchAddress("dead_device", &dead_device);
    t->SetBranchAddress("dead_mask", &dead_mask);
    t->SetBranchAddress("frame_n", &frame_n);
    t->SetBranchAddress("frame", &frame);
    t->SetBranchAddress("trigger0_size", &trigger0_size);
    t->SetBranchAddress("trigger0_n", &trigger0_n);
    t->SetBranchAddress("trigger0_coarse", &trigger0_coarse);
    t->SetBranchAddress("timing_size", &timing_size);
    t->SetBranchAddress("timing_n", &timing_n);
    t->SetBranchAddress("timing_device", &timing_device);
    t->SetBranchAddress("timing_index", &timing_index);
    t->SetBranchAddress("timing_coarse", &timing_coarse);
    t->SetBranchAddress("timing_fine", &timing_fine);
    t->SetBranchAddress("timing_tdc", &timing_tdc);
    t->SetBranchAddress("cherenkov_size", &cherenkov_size);
    t->SetBranchAddress("cherenkov_n", &cherenkov_n);
    t->SetBranchAddress("cherenkov_device", &cherenkov_device);
    t->SetBranchAddress("cherenkov_index", &cherenkov_index);
    t->SetBranchAddress("cherenkov_coarse", &cherenkov_coarse);
    t->SetBranchAddress("cherenkov_fine", &cherenkov_fine);
    t->SetBranchAddress("cherenkov_tdc", &cherenkov_tdc);
  }

  bool
  lightio::next_spill()

  {
    if (spill_current >= tree->GetEntries())
      return false;

    tree->GetEntry(spill_current);
    frame_current = 0;
    trigger0_offset = 0;
    timing_offset = 0;
    cherenkov_offset = 0;

    decode_device_fifos();

    ++spill_current;
    return true;
  }

  bool
  lightio::next_frame()
  {
    if (frame_current >= frame_n)
      return false;

    // fill trigger0 vector
    trigger0_vector.clear();
    for (int i = 0; i < trigger0_n[frame_current]; ++i)
    {
      auto ii = trigger0_offset + i;
      trigger0_vector.push_back(lightdata(0, 0, trigger0_coarse[ii], 0, 0));
    }
    trigger0_offset += trigger0_n[frame_current];

    // fill timing vector and map
    timing_vector.clear();
    timing_map.clear();
    for (int i = 0; i < timing_n[frame_current]; ++i)
    {
      auto ii = timing_offset + i;
      timing_vector.push_back(lightdata(timing_device[ii], timing_index[ii], timing_coarse[ii], timing_fine[ii], timing_tdc[ii]));
      timing_map[{timing_device[ii], timing_index[ii]}].push_back(lightdata(timing_device[ii], timing_index[ii], timing_coarse[ii], timing_fine[ii], timing_tdc[ii]));
    }
    timing_offset += timing_n[frame_current];

    // fill cherenkov vector and map
    cherenkov_vector.clear();
    cherenkov_map.clear();
    for (int i = 0; i < cherenkov_n[frame_current]; ++i)
    {
      auto ii = cherenkov_offset + i;
      cherenkov_vector.push_back(lightdata(cherenkov_device[ii], cherenkov_index[ii], cherenkov_coarse[ii], cherenkov_fine[ii], cherenkov_tdc[ii]));
      cherenkov_map[{cherenkov_device[ii], cherenkov_index[ii]}].push_back(lightdata(cherenkov_device[ii], cherenkov_index[ii], cherenkov_coarse[ii], cherenkov_fine[ii], cherenkov_tdc[ii]));
    }
    cherenkov_offset += cherenkov_n[frame_current];

    ++frame_current;
    return true;
  }

  void
  lightio::decode_part_device_fifos()
  {
    for (auto i_device = 0; i_device < part_n; i_device++)
    {
      auto target_device = (int)part_device[i_device];
      auto device_fifo_mask = (int)part_mask[i_device];
      for (int i_bit = 0; device_fifo_mask > 0; i_bit++)
      {
        if (device_fifo_mask & 1)
          active_device_fifos[target_device][i_bit] = 1;
        else
          active_device_fifos[target_device][i_bit] = 0;
        device_fifo_mask >>= 1;
      }
    }
  }

  void
  lightio::decode_dead_device_fifos()
  {
    for (auto i_device = 0; i_device <= dead_n; i_device++)
    {
      auto target_device = dead_device[i_device];
      auto device_fifo_mask = dead_mask[i_device];
      for (int i_bit = 0; device_fifo_mask > 0; i_bit++)
      {
        if (device_fifo_mask & 1)
          dead_device_fifos[target_device][i_bit] = 1;
        else
          dead_device_fifos[target_device][i_bit] = 0;
        device_fifo_mask >>= 1;
      }
    }
  }

} /** namespace sipm4eic **/
