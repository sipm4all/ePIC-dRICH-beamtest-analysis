#pragma once

#include "lightdata.h"

#define TESTBEAM2023
// #define TESTBEAM2024

namespace sipm4eic
{

#ifdef TESTBEAM2023
  // maps device and chip to pdu and matrix, generated from /etc/drich/drich_readout.conf
  std::map<std::array<int, 2>, std::array<int, 2>> pdu_matrix_map = {
      {{195, 0}, {1, 3}},
      {{195, 1}, {1, 3}},
      {{195, 2}, {1, 4}},
      {{195, 3}, {1, 4}},
      {{193, 0}, {2, 1}},
      {{193, 1}, {2, 1}},
      {{193, 2}, {2, 2}},
      {{193, 3}, {2, 2}},
      {{192, 0}, {2, 3}},
      {{192, 1}, {2, 3}},
      {{192, 2}, {2, 4}},
      {{192, 3}, {2, 4}},
      {{192, 4}, {3, 1}},
      {{192, 5}, {3, 1}},
      {{193, 4}, {3, 2}},
      {{193, 5}, {3, 2}},
      {{194, 4}, {3, 3}},
      {{194, 5}, {3, 3}},
      {{194, 1}, {1, 1}},
      {{194, 2}, {1, 2}},
      {{195, 4}, {3, 4}},
      {{195, 5}, {3, 4}},
      {{196, 0}, {4, 1}},
      {{196, 1}, {4, 1}},
      {{196, 2}, {4, 2}},
      {{196, 3}, {4, 2}},
      {{197, 0}, {4, 3}},
      {{197, 1}, {4, 3}},
      {{197, 2}, {4, 4}},
      {{197, 3}, {4, 4}},
      {{196, 4}, {5, 2}},
      {{196, 5}, {5, 2}},
      {{197, 4}, {6, 2}},
      {{197, 5}, {6, 2}},
      {{198, 0}, {7, 4}},
      {{198, 1}, {7, 4}},
      {{198, 2}, {8, 4}},
      {{198, 3}, {8, 4}}};
#endif
#ifdef TESTBEAM2024
  // automatically generated
  // cat /etc/drich/drich_readout.conf | grep -v "^#" | awk {'print "{ { " substr($4,7,3) " , " $5 " } ,  {" $1 " , " $3 " } } ," '}
  // cat /etc/drich/drich_readout.conf | grep -v "^#" | awk {'print "{ { " substr($4,7,3) " , " $6 " } ,  {" $1 " , " $3 " } } ," '}
  std::map<std::array<int, 2>, std::array<int, 2>> pdu_matrix_map = {
      {{193, 0}, {1, 1}},
      {{193, 2}, {1, 2}},
      {{194, 0}, {1, 3}},
      {{194, 2}, {1, 4}},
      {{195, 0}, {2, 1}},
      {{195, 2}, {2, 2}},
      {{196, 0}, {2, 3}},
      {{196, 2}, {2, 4}},
      {{197, 0}, {3, 1}},
      {{197, 2}, {3, 2}},
      {{198, 0}, {3, 3}},
      {{198, 2}, {3, 4}},
      {{193, 4}, {4, 1}},
      {{194, 4}, {4, 2}},
      {{195, 4}, {4, 3}},
      {{196, 4}, {4, 4}},
      {{197, 4}, {5, 1}},
      {{198, 4}, {5, 2}},
      {{199, 0}, {5, 3}},
      {{199, 2}, {5, 4}},
      {{201, 0}, {6, 1}},
      {{201, 2}, {6, 2}},
      {{202, 0}, {6, 3}},
      {{202, 2}, {6, 4}},
      {{203, 0}, {7, 1}},
      {{203, 2}, {7, 2}},
      {{192, 0}, {7, 3}},
      {{192, 2}, {7, 4}},
      {{201, 4}, {8, 1}},
      {{202, 4}, {8, 2}},
      {{203, 4}, {8, 3}},
      {{192, 4}, {8, 4}},
      {{193, 1}, {1, 1}},
      {{193, 3}, {1, 2}},
      {{194, 1}, {1, 3}},
      {{194, 3}, {1, 4}},
      {{195, 1}, {2, 1}},
      {{195, 3}, {2, 2}},
      {{196, 1}, {2, 3}},
      {{196, 3}, {2, 4}},
      {{197, 1}, {3, 1}},
      {{197, 3}, {3, 2}},
      {{198, 1}, {3, 3}},
      {{198, 3}, {3, 4}},
      {{193, 5}, {4, 1}},
      {{194, 5}, {4, 2}},
      {{195, 5}, {4, 3}},
      {{196, 5}, {4, 4}},
      {{197, 5}, {5, 1}},
      {{198, 5}, {5, 2}},
      {{199, 1}, {5, 3}},
      {{199, 3}, {5, 4}},
      {{201, 1}, {6, 1}},
      {{201, 3}, {6, 2}},
      {{202, 1}, {6, 3}},
      {{202, 3}, {6, 4}},
      {{203, 1}, {7, 1}},
      {{203, 3}, {7, 2}},
      {{192, 1}, {7, 3}},
      {{192, 3}, {7, 4}},
      {{201, 5}, {8, 1}},
      {{202, 5}, {8, 2}},
      {{203, 5}, {8, 3}},
      {{192, 5}, {8, 4}}};
#endif

  /**
      the mapping is a map where the key is the matrix index (U1, U2, U3, U4) and the value is
      a vector of the detector-oriented SiPM index on the matrix for a given
      electronics-oriented channel index

      electronics-oriented channel index is defined as following in the ALCOR-dual boards

      eoch = pixel + 8 * column + 32 * chip

      where pixel = [0, 7] and column = [0, 3] are withing an ALCOR chip
      whereas chip = [0, 1] defines the chip on the ALCOR-dual board (0 = left, 1 = right chip)

      the detector-oriented channel index on the SiPM matrix is defined as follong

      doch = row + column * 8

      where row = [0, 7] and column = [0, 7] are defined starting from the bottom-left corner
  **/

  /**
     maps the readout configuration to the detector
     readout_map[device][chip] --> {pdu, matrix}
  **/

  std::map<int, std::vector<int>> matrix_mapping = {

      {1, {3, 2, 1, 0, 8, 9, 10, 11, 17, 16, 12, 4, 18, 19, 5, 13, 25, 24, 21, 20, 26, 27, 28, 29, 30, 22, 14, 6, 7, 15, 23, 31, 35, 34, 33, 32, 36, 37, 38, 39, 43, 42, 41, 40, 44, 45, 46, 47, 51, 50, 49, 48, 52, 53, 54, 55, 59, 58, 57, 56, 60, 61, 62, 63}},

      {2, {59, 58, 57, 56, 60, 61, 62, 63, 51, 50, 49, 48, 52, 53, 54, 55, 43, 42, 41, 40, 44, 45, 46, 47, 35, 34, 33, 32, 36, 37, 38, 39, 0, 8, 16, 24, 25, 17, 26, 27, 29, 28, 1, 9, 30, 31, 18, 19, 21, 20, 2, 10, 22, 23, 11, 12, 3, 15, 14, 13, 4, 5, 6, 7}},

      {3, {4, 5, 6, 7, 3, 2, 1, 0, 12, 13, 14, 15, 11, 10, 9, 8, 20, 21, 22, 23, 19, 18, 17, 16, 28, 29, 30, 31, 27, 26, 25, 24, 46, 63, 55, 47, 54, 62, 39, 38, 34, 35, 36, 37, 33, 32, 45, 44, 42, 43, 61, 53, 41, 40, 52, 60, 48, 49, 50, 51, 59, 58, 57, 56}},

      {4, {60, 61, 62, 63, 55, 54, 53, 52, 46, 47, 51, 59, 45, 44, 58, 50, 38, 39, 42, 43, 37, 36, 35, 34, 33, 41, 49, 57, 56, 48, 40, 32, 28, 29, 30, 31, 27, 26, 25, 24, 20, 21, 22, 23, 19, 18, 17, 16, 12, 13, 14, 15, 11, 10, 9, 8, 4, 5, 6, 7, 3, 2, 1, 0}}};

#ifdef TESTBEAM2023
  bool rotateme[8] = {true, true, true, true, false, true, true, false};
#endif
#ifdef TESTBEAM2024
  bool rotateme[8] = {true, true, true, true, true, true, true, true};
#endif

#ifdef TESTBEAM2023
  std::map<int, int> placement = {
      {6, 1}, {4, 2}, {7, 3}, {2, 4}, {1, 6}, {8, 7}, {3, 8}, {5, 9}};
  std::map<int, std::array<float, 2>> placement_xy = {
      {6, {-82., 30.}}, {4, {-26., 35.}}, {7, {30., 30.}}, {2, {-82., -26.}}, {1, {30., -26.}}, {8, {-82., -82.}}, {3, {-26., -87.}}, {5, {30., -82.}}};
#endif
#ifdef TESTBEAM2024
  std::map<int, int> placement = {
      {1, 1}, {2, 2}, {3, 3}, {8, 4}, {4, 6}, {7, 7}, {6, 8}, {5, 9}};
  std::map<int, std::array<float, 2>> placement_xy = {
      {1, {-82., 30.}}, {2, {-26., 35.}}, {3, {30., 30.}}, {8, {-82., -26.}}, {4, {30., -26.}}, {7, {-82., -82.}}, {6, {-26., -87.}}, {5, {30., -82.}}};
#endif

  const std::array<float, 2> position_offset = {1.7, 1.7}; // the centre of the SiPM in the bottom-left corner (A1)
  const std::array<float, 2> position_pitch = {3.2, 3.2};  // the distance between the SiPM cetres

  int get_do_channel(int matrix, int eo_channel) { return matrix_mapping[matrix][eo_channel]; }

  /*******************************************************************************/

  std::array<int, 3>
  get_geo(int pdu, int matrix, int eo_channel)
  {
    auto do_channel = get_do_channel(matrix, eo_channel);
    std::array<int, 3> geo = {pdu, do_channel / 8, do_channel % 8};

    if (matrix == 2 || matrix == 4)
      geo[2] += 8;
    if (matrix == 3 || matrix == 4)
      geo[1] += 8;

    if (rotateme[pdu - 1])
    {
      geo[1] = 15 - geo[1];
      geo[2] = 15 - geo[2];
    }

    return geo;
  }

  std::array<int, 3>
  get_geo(lightdata cherenkov)
  {
    int device = cherenkov.device;
    auto chip = cherenkov.chip();
    auto eoch = cherenkov.eoch();
    auto pdu = pdu_matrix_map[{device, chip}][0];
    auto matrix = pdu_matrix_map[{device, chip}][1];
    return get_geo(pdu, matrix, eoch);
  }

  std::vector<std::array<int, 3>>
  get_fifo_geo(int device, int fifo)
  {
    std::vector<std::array<int, 3>> result;
    int target_device = device;
    int target_chip = fifo / 4;
    int target_chip_fifo = fifo % 4;
    auto target_pdu_matrix = pdu_matrix_map[{target_device, target_chip}];
    if (target_chip > 5 || target_pdu_matrix[0] == 0 || target_pdu_matrix[1] == 0)
      return {{-1, -1, -1}};

    for (auto i_eoch = 0; i_eoch < 8; i_eoch++)
      result.push_back(get_geo(target_pdu_matrix[0], target_pdu_matrix[1], (target_chip * 32 + target_chip_fifo * 8 + i_eoch) % 64));
    return result;
  }

  /*******************************************************************************/

  std::array<float, 2>
  get_position(std::array<int, 3> geo)
  {
    std::array<int, 3> null_geo = {{-1, -1, -1}};
    if (null_geo == geo)
      return {-999, -999};
    int pdu = geo[0];
    int col = geo[1];
    int row = geo[2];
    float x = 0.05 + 0.1 + 0.2 + 1.5 + 3.2 * col; // not clear why the 0.05, but it works to center the thing
    float y = 0.05 + 0.1 + 0.2 + 1.5 + 3.2 * row;
    if (col > 7)
      x += 0.3;
    if (row > 7)
      y += 0.3;
    x += placement_xy[pdu][0];
    y += placement_xy[pdu][1];

    return {x, y};
  }

  std::vector<std::array<float, 2>>
  get_fifo_position(int device, int fifo)
  {
    std::vector<std::array<float, 2>> result;
    auto fifo_geo = get_fifo_geo(device, fifo);
    for (auto current_geo : fifo_geo)
      result.push_back(get_position(current_geo));
    return result;
  }

  void decode_good_positions(lightio &io)
  {
    //  Active
    for (auto [device, fifo] : io.get_active_fifos())
      for (auto [i_fifo, active_fifo] : fifo)
        for (auto current_position : get_fifo_position(device, i_fifo))
          io.set_active_position(current_position, active_fifo);

    //  Dead
    for (auto [device, fifo] : io.get_dead_fifos())
      for (auto [i_fifo, active_fifo] : fifo)
        for (auto current_position : get_fifo_position(device, i_fifo))
          io.set_dead_position(current_position, active_fifo);

    //  Good
    for (auto [device, fifo] : io.get_good_fifos())
      for (auto [i_fifo, active_fifo] : fifo)
        for (auto current_position : get_fifo_position(device, i_fifo))
          io.set_good_position(current_position, active_fifo);
  }

}