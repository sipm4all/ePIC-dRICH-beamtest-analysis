#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/macros/recowriter.C"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/include/fit_functions.h"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/macro/test_bkg_investigation_step2.C"

/*
std::vector<std::string> no_aerogel_runs = {
    //"20240530-194756",
    "20240530-200252",
    "20240530-202344",
    "20240530-204413",
    "20240530-210507",
    "20240530-212533"};

//  --- Std first run
std::string run = "20240526-220929"; // 11.5 GeV
*/

template <typename T>
int find_index_of(const std::vector<T> &vec, const T &value)
{
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end())
    {
        return static_cast<int>(std::distance(vec.begin(), it));
    }
    return -1; // Not found
}

void baseline_dcr_subtraction(std::string current_run = run)
{
    //  Load lightdata from run
    std::string filename = Form("/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/Data/%s/lightdata.root", current_run.c_str());

    //  Read data
    sipm4eic::lightio io;
    io.read_from_tree(filename);

    //  Time difference, sensor by sensor
    std::map<std::array<float, 2>, TH1F *> h_delta_time_sensors;
    std::map<std::array<float, 2>, TH2F *> h_delta_time_sensors_fine;
    std::map<std::array<float, 2>, std::array<float, 2>> h_sensors_fine_parameters;
    std::vector<std::array<float, 2>> position_map_calibration_index;

    //  --- Overall sensors participation per spill
    TH2F *h_good_sensors_per_spill = new TH2F("h_good_sensors_per_spill", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);

    //  Loop over spills
    auto n_spills = 0;
    auto n_frames = 0;
    while (io.next_spill())
    {
        //  Spills counter
        n_spills++;

        //  Get positions and available positions for spill
        sipm4eic::decode_good_positions(io);
        auto good_positions = io.get_good_positions();
        for (auto [position, state] : good_positions)
        {
            if (!h_delta_time_sensors[position])
                h_delta_time_sensors[position] = new TH1F(Form("h_delta_time_sensors_%.2f_%.2f", position[0], position[1]), ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 5000, -250, 250);
            if (!h_delta_time_sensors_fine[position])
                h_delta_time_sensors_fine[position] = new TH2F(Form("h_delta_time_sensors_fine_%.2f_%.2f", position[0], position[1]), ";fine (a.u.);#Delta_{t_{hit}-t_{trg}}", 150, 0, 150, 5000, -250, 250);
            h_sensors_fine_parameters[position] = {1. / 60, 2. / 3.};
            if (find_index_of(position_map_calibration_index, position) == -1)
                position_map_calibration_index.push_back(position);
        }

        //  Loop over frames
        auto n_spill_frames = 0;
        while (io.next_frame())
        {
            //  Frames counters
            n_frames++;
            n_spill_frames++;

            // Trigger time reference
            auto trigger0_vector = io.get_trigger0_vector();
            auto reference_time = trigger0_vector[0].coarse;
            auto reference_time_fine = trigger0_vector[0].fine;

            // loop over cherenkov hits
            auto cherenkov_map = io.get_cherenkov_map();
            for (auto &[index, hits] : cherenkov_map)
            {
                //  Sort hits
                std::sort(hits.begin(), hits.end());

                //  Take first hit, syntax to take all hits
                for (auto i_hit = 0; i_hit < /*hits.size()*/ 1; i_hit++)
                {
                    //  Get hit
                    auto current_hit = hits[i_hit];

                    // Get position
                    int device = current_hit.device;
                    auto chip = current_hit.chip();
                    auto pdu = sipm4eic::pdu_matrix_map[{device, chip}][0] - 1;
                    auto matrix = sipm4eic::pdu_matrix_map[{device, chip}][1];
                    auto geo = sipm4eic::get_geo(current_hit);
                    auto pos = sipm4eic::get_position(geo);

                    // Get time information
                    auto coarse = current_hit.coarse;
                    auto fine = current_hit.fine;
                    auto time_difference = (coarse - (fine * h_sensors_fine_parameters[pos][0] + h_sensors_fine_parameters[pos][1])) - reference_time;

                    //  Fill reference plots
                    h_delta_time_sensors[pos]->Fill(time_difference);
                    h_delta_time_sensors_fine[pos]->Fill(fine, time_difference);
                }
            }
        }

        // Track each sensor availability throughout the spills
        for (auto [position, state] : good_positions)
        {
            auto x_scatter = position[0];
            auto y_scatter = position[1];
            if (x_scatter == -999 && y_scatter == -999)
                continue;
            h_good_sensors_per_spill->Fill(x_scatter, y_scatter, n_spill_frames);
        }
    }

    //  Guesstimate the fine correction
    for (auto [pos, hist] : h_delta_time_sensors_fine)
    {
        //  Exclude problematic positons
        if (pos[0] == -999 || pos[1] == -999)
            continue;
        auto min_val = hist->GetXaxis()->GetBinLowEdge(hist->FindFirstBinAbove(0.));
        auto max_val = hist->GetXaxis()->GetBinLowEdge(hist->FindLastBinAbove(0.));
    }

    //  Minimisation function
    auto calibration_sigma_chi2 = [&](const double *parameters)
    {
        //  Define result
        double current_chi2 = 0.;

        // re-loop over spills
        while (io.next_spill())
            while (io.next_frame())
            {
                // Trigger time reference
                auto trigger0_vector = io.get_trigger0_vector();
                auto reference_time = trigger0_vector[0].coarse;

                // loop over cherenkov hits
                auto cherenkov_map = io.get_cherenkov_map();
                for (auto &[index, hits] : cherenkov_map)
                {
                    //  Sort hits
                    std::sort(hits.begin(), hits.end());

                    //  Take first hit, syntax to take all hits
                    for (auto i_hit = 0; i_hit < hits.size(); i_hit++)
                    {
                        //  Get hit
                        auto current_hit = hits[i_hit];

                        // Get position
                        int device = current_hit.device;
                        auto chip = current_hit.chip();
                        auto pdu = sipm4eic::pdu_matrix_map[{device, chip}][0] - 1;
                        auto matrix = sipm4eic::pdu_matrix_map[{device, chip}][1];
                        auto geo = sipm4eic::get_geo(current_hit);
                        auto pos = sipm4eic::get_position(geo);
                        auto fine = current_hit.fine;
                        auto coarse = current_hit.coarse;
                        int universal_index = current_hit.cindex();

                        // Get time information
                        auto coarse = current_hit.coarse;
                        double current_correction = parameters[] + (double)current_hit.fine * parameters[];
                        auto time_difference = coarse - current_correction - reference_time;

                        current_chi2 += time_difference * time_difference;
                    }
                }
            }

        return current_chi2;
    }

    ROOT::Math::Functor minimisation_function(periodic_sigma_chi2, 9);
    ROOT::Fit::Fitter fitter;
    double pStart[9] = {0.5, 0.5, 0.5, 0.5, 0.0156, 0.0156, 0.0156, 0.0156, 320.};
    fitter.SetFCN(minimisation_function, pStart);
    fitter.Config().ParSettings(0).SetName("off_0");
    fitter.Config().ParSettings(1).SetName("off_1");
    fitter.Config().ParSettings(2).SetName("off_2");
    fitter.Config().ParSettings(3).SetName("off_3");
    fitter.Config().ParSettings(4).SetName("iif_0");
    fitter.Config().ParSettings(5).SetName("iif_1");
    fitter.Config().ParSettings(6).SetName("iif_2");
    fitter.Config().ParSettings(7).SetName("iif_3");
    */

        //  Save histograms
        auto fout = TFile::Open(Form("baseline_dcr_subtraction.%s.root", current_run.c_str()), "RECREATE");
    for (auto [pos, hist] : h_delta_time_sensors)
        hist->Write();
    for (auto [pos, hist] : h_delta_time_sensors_fine)
        hist->Write();
    fout->Close();
}

/*
 // Elaborate the time coincidences
 std::map<std::array<float, 2>, std::array<double, 2>> h_delta_time_sensors_signal;
 std::map<std::array<float, 2>, std::array<double, 2>> h_delta_time_sensors_background;
 for (auto &[position, current_histogram] : h_delta_time_sensors)
 {
     //  Exclude problematic positons
     if (position[0] == -999 || position[1] == -999)
         continue;

     //  Scale per frames
     current_histogram->Scale(1. / n_frames);

     double background_zone_1_err, background_signal_err, background_zone_2_err;
     auto background_zone_1 = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(-30.), current_histogram->GetXaxis()->FindBin(-10.), background_zone_1_err);
     auto background_signal = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(-10.), current_histogram->GetXaxis()->FindBin(+30.), background_signal_err);
     auto background_zone_2 = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(+30.), current_histogram->GetXaxis()->FindBin(+50.), background_zone_2_err);
     auto final_background = (background_zone_1 + background_zone_2) / 40;
     auto final_background_err = TMath::Sqrt(background_zone_1_err * background_zone_1_err + background_zone_2_err * background_zone_2_err) / 40;
     auto final_signal = background_signal / 40 - final_background;
     auto final_signal_err = background_signal_err / 40;
     auto final_signal_err2 = TMath::Sqrt(background_signal_err * background_signal_err + final_background_err * final_background_err);

     //  Fill histograms
     //  --- Fill signal histogram
     auto target_bin = h_signal_sensor->FindBin(position[0], position[1]);
     h_signal_sensor->SetBinContent(target_bin, final_signal);
     h_signal_sensor->SetBinError(target_bin, final_signal_err2);
     auto current_point = g_signal_sensor->GetN();
     g_signal_sensor->SetPoint(current_point, position[0], position[1], final_signal);
     g_signal_sensor->SetPointError(current_point, 0., 0., final_signal_err2);
     //  --- Fill background histogram
     target_bin = h_background_sensor->FindBin(position[0], position[1]);
     h_background_sensor->SetBinContent(target_bin, final_background);
     h_background_sensor->SetBinError(target_bin, final_background_err);
     current_point = g_background_sensor->GetN();
     g_background_sensor->SetPoint(current_point, position[0], position[1], final_background);
     g_background_sensor->SetPointError(current_point, 0., 0., final_background_err);
     //  --- Fill signal + background histogram
     target_bin = h_coincidences_per_sensor->FindBin(position[0], position[1]);
     h_coincidences_per_sensor->SetBinContent(target_bin, background_signal / 40);
     h_coincidences_per_sensor->SetBinError(target_bin, background_signal_err / 40);

     //  Fill signal container
     h_sensors_signal[position] = {final_signal, final_signal_err2};
     h_sensors_background[position] = {final_background, final_background_err};
 }
     */

void test_bkg_investigation____(std::string current_run = run) // Data/20231011-221500/ // Data/20240530-200252/
{
    std::string filename = Form("/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/Data/%s/lightdata.root", current_run.c_str());

    //  Read data
    sipm4eic::lightio io;
    io.read_from_tree(filename);

    //  Random generator
    TRandom *g_random = new TRandom();

    //  Vector of zone boundaries
    std::pair<std::array<float, 2>, std::vector<std::array<float, 2>>> zones_delimiter = {{0., 0.}, {{0., 60.}, {60., 100.}, {100., 999.}}};

    //  Output histograms
    //  --- Zones
    //  --- --- Participants counting
    TH1F *h_participants_per_zone = new TH1F("h_participants_per_zone", "", zones_delimiter.second.size() + 1, 0, zones_delimiter.second.size() + 1);
    //  --- Hitmaps
    //  --- --- Masks
    TH2F *h_good_sensors_per_spill = new TH2F("h_good_sensors_per_spill", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);
    TH2F *h_coincidences_per_sensor = new TH2F("h_coincidences_per_sensor", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);
    TH2F *h_background_per_sensor = new TH2F("h_background_per_sensor", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);
    //  --- --- Sensors
    TH2F *h_signal_sensor = new TH2F("h_signal_sensor", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);
    TGraph2DErrors *g_signal_sensor = new TGraph2DErrors();
    g_signal_sensor->SetName("g_signal_sensor");
    g_signal_sensor->SetTitle(";x (mm);y (mm)");
    TH2F *h_background_sensor = new TH2F("h_background_sensor", ";x (mm);y (mm)", 198, -99, 99, 198, -99, 99);
    TGraph2DErrors *g_background_sensor = new TGraph2DErrors();
    g_background_sensor->SetName("g_background_sensor");
    g_background_sensor->SetTitle(";x (mm);y (mm)");
    //  --- Reference time
    //  --- --- Zones
    std::map<std::array<float, 2>, TH1F *> h_delta_time_sensors;
    TProfile *h_signal_sensor_radius = new TProfile("h_signal_sensor_radius", ";radius (mm);", 150, 0, 150);
    TProfile *h_signal_sensor_radius_fit = new TProfile("h_signal_sensor_radius_fit", ";radius (mm);", 150, 0, 150);
    TProfile *h_signal_sensor0_radius = new TProfile("h_signal_sensor0_radius", ";radius (mm);", 150, 0, 150);
    TProfile *h_signal_sensor1_radius = new TProfile("h_signal_sensor1_radius", ";radius (mm);", 150, 0, 150);
    TProfile *h_background_sensor_radius = new TProfile("h_background_sensor_radius", ";radius (mm);", 150, 0, 150);
    TProfile *h_background_sensor0_radius = new TProfile("h_background_sensor0_radius", ";radius (mm);", 150, 0, 150);
    TProfile *h_background_sensor1_radius = new TProfile("h_background_sensor1_radius", ";radius (mm);", 150, 0, 150);

    TH1F *h_delta_time = new TH1F("h_delta_time", ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 500, -250, 250);
    TH1F *h_delta_time_zone1 = new TH1F("h_delta_time_zone1", ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 500, -250, 250);
    TH1F *h_delta_time_zone2 = new TH1F("h_delta_time_zone2", ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 500, -250, 250);
    TH1F *h_delta_time_zone3 = new TH1F("h_delta_time_zone3", ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 500, -250, 250);

    //  Loop over spills
    auto n_spills = 0;
    auto n_frames = 0;
    while (io.next_spill())
    {
        //  Spills counter
        n_spills++;

        //  Get positions and available positions for spill
        sipm4eic::decode_good_positions(io);
        auto good_positions = io.get_good_positions();
        for (auto [position, state] : good_positions)
        {
            if (!h_delta_time_sensors[position])
                h_delta_time_sensors[position] = new TH1F(Form("h_delta_time_sensors_%.2f_%.2f", position[0], position[1]), ";#Delta_{t_{hit}-t_{trg}};Hits per frame", 500, -250, 250);
        }

        for (auto [position, state] : good_positions)
        {
            // Track each sensor availability throughout the spills
            auto x_scatter = position[0];
            auto y_scatter = position[1];
            if (x_scatter == -999 && y_scatter == -999)
                continue;
            h_good_sensors_per_spill->Fill(x_scatter, y_scatter);
        }

        //  Loop over frames
        auto n_spill_frames = 0;
        while (io.next_frame())
        {
            //  Frames counters
            n_frames++;
            n_spill_frames++;

            // Trigger time reference
            auto trigger0_vector = io.get_trigger0_vector();
            auto reference_time = trigger0_vector[0].coarse;

            // loop over cherenkov hits
            auto cherenkov_map = io.get_cherenkov_map();
            for (auto &[index, hits] : cherenkov_map)
            {
                //  Sort hits
                std::sort(hits.begin(), hits.end());

                //  Take first hit, syntax to take all hits
                for (auto i_hit = 0; i_hit < hits.size(); i_hit++)
                {
                    //  Get hit
                    auto current_hit = hits[i_hit];

                    // Get position
                    int device = current_hit.device;
                    auto chip = current_hit.chip();
                    auto pdu = sipm4eic::pdu_matrix_map[{device, chip}][0] - 1;
                    auto matrix = sipm4eic::pdu_matrix_map[{device, chip}][1];
                    auto geo = sipm4eic::get_geo(current_hit);
                    auto pos = sipm4eic::get_position(geo);

                    // Get time information
                    auto coarse = current_hit.coarse;
                    auto time_difference = coarse - reference_time;

                    //  Avoid border effectsFill reference plots
                    h_delta_time->Fill(time_difference);
                    h_delta_time_sensors[pos]->Fill(time_difference);
                }
            }
        }
    }

    //  Process the delta time histograms
    std::map<std::array<float, 2>, std::array<double, 2>> h_sensors_signal;
    std::map<std::array<float, 2>, std::array<double, 2>> h_sensors_background;
    for (auto &[position, current_histogram] : h_delta_time_sensors)
    {
        //  Exclude problematic positons
        if (position[0] == -999 || position[1] == -999)
            continue;

        //  Scale per frames
        current_histogram->Scale(1. / n_frames);

        double background_zone_1_err, background_signal_err, background_zone_2_err;
        auto background_zone_1 = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(-30.), current_histogram->GetXaxis()->FindBin(-10.), background_zone_1_err);
        auto background_signal = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(-10.), current_histogram->GetXaxis()->FindBin(+30.), background_signal_err);
        auto background_zone_2 = 1.e3 * current_histogram->IntegralAndError(current_histogram->GetXaxis()->FindBin(+30.), current_histogram->GetXaxis()->FindBin(+50.), background_zone_2_err);
        auto final_background = (background_zone_1 + background_zone_2) / 40;
        auto final_background_err = TMath::Sqrt(background_zone_1_err * background_zone_1_err + background_zone_2_err * background_zone_2_err) / 40;
        auto final_signal = background_signal / 40 - final_background;
        auto final_signal_err = background_signal_err / 40;
        auto final_signal_err2 = TMath::Sqrt(background_signal_err * background_signal_err + final_background_err * final_background_err);

        //  Fill histograms
        //  --- Fill signal histogram
        auto target_bin = h_signal_sensor->FindBin(position[0], position[1]);
        h_signal_sensor->SetBinContent(target_bin, final_signal);
        h_signal_sensor->SetBinError(target_bin, final_signal_err2);
        auto current_point = g_signal_sensor->GetN();
        g_signal_sensor->SetPoint(current_point, position[0], position[1], final_signal);
        g_signal_sensor->SetPointError(current_point, 0., 0., final_signal_err2);
        //  --- Fill background histogram
        target_bin = h_background_sensor->FindBin(position[0], position[1]);
        h_background_sensor->SetBinContent(target_bin, final_background);
        h_background_sensor->SetBinError(target_bin, final_background_err);
        current_point = g_background_sensor->GetN();
        g_background_sensor->SetPoint(current_point, position[0], position[1], final_background);
        g_background_sensor->SetPointError(current_point, 0., 0., final_background_err);
        //  --- Fill signal + background histogram
        target_bin = h_coincidences_per_sensor->FindBin(position[0], position[1]);
        h_coincidences_per_sensor->SetBinContent(target_bin, background_signal / 40);
        h_coincidences_per_sensor->SetBinError(target_bin, background_signal_err / 40);

        //  Fill signal container
        h_sensors_signal[position] = {final_signal, final_signal_err2};
        h_sensors_background[position] = {final_background, final_background_err};
    }

    //  Fit signal histogram
    TF2 *h_find_center_single = new TF2("h_find_center", _h_find_center, -100, 100, -100, 100, 9);
    h_find_center_single->SetParNames("x0", "y0", "R", "sigma", "amplitude", "bkg_amplitude", "bkg_power", "bkg_drop", "bkg_drop_sigma");
    h_find_center_single->SetParLimits(0, -3, 3);    // x0
    h_find_center_single->SetParLimits(1, -3, 3);    // y0
    h_find_center_single->SetParLimits(2, 0, 100);   // R
    h_find_center_single->SetParLimits(3, 0, 5);     // sigma
    h_find_center_single->SetParLimits(4, 0, 20);    // amplitude
    h_find_center_single->SetParLimits(5, 0, 10);    // bkg_amplitude
    h_find_center_single->SetParLimits(6, -2, 0);    // bkg_power
    h_find_center_single->SetParLimits(7, 80, 120.); // bkg_drop
    h_find_center_single->SetParLimits(8, 5, 15.);   // bkg_drop_sigma
    h_find_center_single->SetParameters(0.0, 0.0, 75.0, 2.0, 10.0, 2.0, -1.1, 100.0, 10.0);
    TF2 *h_find_center_double = new TF2("h_find_center", _h_find_center_2, -100, 100, -100, 100, 11);
    h_find_center_double->SetParNames("x0", "y0", "R", "sigma", "amplitude", "bkg_amplitude", "bkg_power", "bkg_drop", "bkg_drop_sigma", "2nd_R", "2nd_amplitude");
    h_find_center_double->SetParLimits(0, -3, 3);    // x0
    h_find_center_double->SetParLimits(1, -3, 3);    // y0
    h_find_center_double->SetParLimits(2, 0, 100);   // R
    h_find_center_double->SetParLimits(3, 0, 5);     // sigma
    h_find_center_double->SetParLimits(4, 0, 20);    // amplitude
    h_find_center_double->SetParLimits(5, 0, 10);    // bkg_amplitude
    h_find_center_double->SetParLimits(6, -2, 0);    // bkg_power
    h_find_center_double->SetParLimits(7, 80, 120.); // bkg_drop
    h_find_center_double->SetParLimits(8, 5, 15.);   // bkg_drop_sigma
    h_find_center_double->SetParLimits(9, 0, 100);   // 2nd R
    h_find_center_double->SetParLimits(10, 0, 20);   // 2nd amplitude
    h_find_center_double->SetParameters(0.0, 0.0, 75.0, 2.0, 15.0, 2.0, -1.1, 100.0, 10.0, 65.0, 10.0);
    TF2 *h_find_center_triple = new TF2("h_find_center", _h_find_center_3, -100, 100, -100, 100, 13);
    h_find_center_triple->SetParNames("x0", "y0", "R", "sigma", "amplitude", "bkg_amplitude", "bkg_power", "bkg_drop", "bkg_drop_sigma", "2nd_R", "2nd_amplitude");
    h_find_center_double->SetParLimits(0, -3, 3);    // x0
    h_find_center_triple->SetParLimits(1, -3, 3);    // y0
    h_find_center_triple->SetParLimits(2, 0, 100);   // R
    h_find_center_triple->SetParLimits(3, 0, 5);     // sigma
    h_find_center_triple->SetParLimits(4, 0, 20);    // amplitude
    h_find_center_triple->SetParLimits(5, 0, 10);    // bkg_amplitude
    h_find_center_triple->SetParLimits(6, -2, 0);    // bkg_power
    h_find_center_triple->SetParLimits(7, 80, 120.); // bkg_drop
    h_find_center_triple->SetParLimits(8, 5, 15.);   // bkg_drop_sigma
    h_find_center_triple->SetParLimits(9, 0, 100);   // 2nd R
    h_find_center_triple->SetParLimits(10, 0, 20);   // 2nd amplitude
    h_find_center_triple->SetParLimits(11, 0, 100);  // 3rd R
    h_find_center_triple->SetParLimits(12, 0, 20);   // 3rd amplitude
    h_find_center_triple->SetParameters(0.0, 0.0, 75.0, 2.0, 15.0, 2.0, -1.1, 100.0, 10.0, 65.0, 10.0);
    h_find_center_triple->SetParameter(11, 60); // 3rd R
    h_find_center_triple->SetParameter(12, 2);  // 3rd amplitude

    TF2 *h_find_center = h_find_center_triple;

    TCanvas *c_show_fit = new TCanvas();
    g_signal_sensor->Draw("COLZ");
    g_signal_sensor->Fit(h_find_center);

    TCanvas *c_show_bkg = new TCanvas();
    h_background_per_sensor->Draw("COLZ");

    auto center_x = h_find_center->GetParameter(0);
    auto center_y = h_find_center->GetParameter(1);

    for (auto &[position, signal] : h_sensors_signal)
    {
        auto current_radius = std::sqrt((position[0] - center_x) * (position[0] - center_x) + (position[1] - center_y) * (position[1] - center_y));
        auto current_radius_bkg = std::sqrt(position[0] * position[0] + position[1] * position[1]);

        h_signal_sensor_radius->Fill(current_radius, signal[0]);
        h_background_sensor_radius->Fill(current_radius_bkg, h_sensors_background[position][0]);

        if (which_sensor(position) == 0)
        {
            h_signal_sensor0_radius->Fill(current_radius, signal[0]);
            h_background_sensor0_radius->Fill(current_radius_bkg, h_sensors_background[position][0]);
        }
        else
        {
            h_signal_sensor1_radius->Fill(current_radius, signal[0]);
            h_background_sensor1_radius->Fill(current_radius_bkg, h_sensors_background[position][0]);
        }
        h_signal_sensor_radius_fit->Fill(current_radius, h_find_center->Eval(position[0], position[1]));
    }

    TCanvas *c_show_fit2 = new TCanvas();
    h_signal_sensor_radius->Draw("");
    h_signal_sensor_radius->SetMarkerStyle(20);
    h_signal_sensor_radius_fit->Draw("HIST SAME");
    h_signal_sensor_radius_fit->SetLineColor(kRed);

    //  Save histograms
    auto fout = TFile::Open(Form("test_bkg_investigation.%s.root", current_run.c_str()), "RECREATE");
    g_signal_sensor->Write();
    h_signal_sensor->Write();
    h_signal_sensor0_radius->Write();
    h_signal_sensor1_radius->Write();
    h_signal_sensor_radius->Write();
    g_background_sensor->Write();
    h_background_sensor->Write();
    h_background_sensor0_radius->Write();
    h_background_sensor1_radius->Write();
    h_background_sensor_radius->Write();
    for (auto [pos, hist] : h_delta_time_sensors)
        hist->Write();

    fout->Close();
}