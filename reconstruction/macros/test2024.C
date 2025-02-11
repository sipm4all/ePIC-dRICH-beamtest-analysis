#include "../lib/alcor.h"

void test2024()
{
    std::string target_file = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/tracking_data.root";
    alcor::trkdata data(target_file);
    auto current_tree = data.get_tree();

    std::vector<std::array<std::array<unsigned char, 4>, 2>> eligible_pairs;
    std::vector<std::array<std::array<unsigned char, 4>, 2>> test_pairs;
    std::map<int, int> enough_entries;

    //  Front-back coincidences
    TH1F *uncorrected_coincidence = new TH1F("uncorrected_coincidence", ";t_{hit_{1}} - t_{hit_{2}} (clock cycle)", 20, -15, 5);
    TH1F *corrected_coincidence = new TH1F("corrected_coincidence", ";t_{hit_{1}} - t_{hit_{2}} (clock cycle, fine corr.)", 200, -10, 10);

    TH2F *fine_test = new TH2F("fine_test", ";fine;t_{hit} - t_{ref}", 100, 25, 125, 50, -10, 10);
    TProfile *fine_test_p = new TProfile("fine_test_p", ";fine", 100, 25, 125);

    TH2F *fine_test_2 = new TH2F("fine_test_2", ";fine;t_{hit} - t_{ref}", 100, 25, 125, 100, -20, 20);
    TH1F *check_channels = new TH1F("check_channels", "", 300, 0, 300);

    for (auto iev = 0; iev < current_tree->GetEntries(); iev++)
    {
        current_tree->GetEntry(iev);
        if (data.get_in_frame_hits() != 2)
            continue;
        for (auto ihit = 0; ihit < data.get_in_frame_hits(); ihit++)
            for (auto jhit = ihit + 1; jhit < data.get_in_frame_hits(); jhit++)
            {
                if (data.get_chip(ihit) != data.get_chip(jhit))
                {
                    if (fabs(6 + data.get_coarse(ihit) - data.get_coarse(jhit)) > 4)
                         continue;
                    if (data.get_index(ihit) != 101)
                        continue;
                    if (data.get_index(jhit) != 163)
                        continue;
                    if (fabs(data.get_fine(ihit) -65)< 5)
                        continue;
                    if (fabs(data.get_fine(jhit) -65)< 5)
                        continue;
                    auto i_index = data.get_index(ihit) < 140 ? data.get_index(ihit) % 32 : data.get_index(ihit) % 32 + 32;
                    auto j_index = data.get_index(jhit) < 140 ? data.get_index(jhit) % 32 : data.get_index(jhit) % 32 + 32;
                    test_pairs.push_back({{{data.get_tdc(ihit), data.get_index(ihit), data.get_coarse(ihit), data.get_fine(ihit)}, {data.get_tdc(jhit), data.get_index(jhit), data.get_coarse(jhit), data.get_fine(jhit)}}});
                    auto i_time = data.get_coarse(ihit) - (0.43 + 0.0143 * data.get_fine(ihit) - 1 * (data.get_fine(ihit) > 65));
                    auto j_time = data.get_coarse(jhit) - (0.43 + 0.0143 * data.get_fine(jhit) - 1 * (data.get_fine(jhit) > 65));
                    //if (fabs(i_time - j_time + 6) > 0.6)
                    //    continue;
                    enough_entries[i_index]++;
                    enough_entries[j_index]++;
                    eligible_pairs.push_back({{{data.get_tdc(ihit), data.get_index(ihit), data.get_coarse(ihit), data.get_fine(ihit)}, {data.get_tdc(jhit), data.get_index(jhit), data.get_coarse(jhit), data.get_fine(jhit)}}});
                }
            }
    }

    //  Minimisation function
    auto periodic_sigma_chi2 = [&](const double *parameters)
    {
        //  Define result
        double current_chi2 = 0.;

        //  current_pair (tdc, index, coarse, fine)
        for (auto current_pair : eligible_pairs)
        {
            //  fix "a manoni" to map everything from 1 to 513
            auto i_index = current_pair[0][1] < 140 ? current_pair[0][1] % 32 : current_pair[0][1] % 32 + 32;
            auto j_index = current_pair[1][1] < 140 ? current_pair[1][1] % 32 : current_pair[1][1] % 32 + 32;

            //  define a parameter index
            auto i_par_off_index = 1 + i_index * 4 * 3 + current_pair[0][0] * 3 + 0;
            auto i_par_slope_index = 1 + i_par_off_index;
            auto i_par_cut_index = 2 + i_par_off_index;
            auto j_par_off_index = 1 + j_index * 4 * 3 + current_pair[1][0] * 3 + 0;
            auto j_par_slope_index = 1 + j_par_off_index;
            auto j_par_cut_index = 2 + j_par_off_index;

            //  define the parameters involved
            //  define the parameters involved
            auto i_off = parameters[i_par_off_index];
            auto i_slope = parameters[i_par_slope_index];
            auto i_cut = parameters[i_par_cut_index];
            auto j_off = parameters[j_par_off_index];
            auto j_slope = parameters[j_par_slope_index];
            auto j_cut = parameters[j_par_cut_index];

            //  Compute time
            auto i_time = current_pair[0][2] - (i_off + i_slope * current_pair[0][3] - 1 * (current_pair[0][3] > i_cut));
            auto j_time = current_pair[1][2] - (j_off + j_slope * current_pair[1][3] - 1 * (current_pair[1][3] > j_cut));

            //  Compute the delta of times
            double current_delta = i_time - j_time + parameters[0];
            current_chi2 += current_delta * current_delta;
        }
        return current_chi2;
    };

    ROOT::Math::Functor minimisation_function(periodic_sigma_chi2, 769);
    ROOT::Fit::Fitter fitter;

    double pStart[769];
    pStart[0] = 5;
    for (auto iter = 1; iter < 769; iter += 3)
    {
        pStart[iter] = 0.5;
        pStart[iter + 1] = +0.015;
        pStart[iter + 2] = 65;
    }
    fitter.SetFCN(minimisation_function, pStart);
    fitter.Config().ParSettings(0).SetName("t_diff");
    for (auto iter = 1; iter < 769; iter += 3)
    {
        fitter.Config().ParSettings(iter).SetName(Form("index_%i_tdc_%i_off", iter / 12, (iter % 12) / 3));
        fitter.Config().ParSettings(iter + 1).SetName(Form("index_%i_tdc_%i_slope", iter / 12, (iter % 12) / 3));
        fitter.Config().ParSettings(iter + 2).SetName(Form("index_%i_tdc_%i_cut", iter / 12, (iter % 12) / 3));
        if (enough_entries[(iter - 1) / 12] < 1000000)
        {
            fitter.Config().ParSettings(iter).Fix();
            fitter.Config().ParSettings(iter + 1).Fix();
            fitter.Config().ParSettings(iter + 2).Fix();
        }
        else
        {
            cout << "iter: " << iter << " - name: " << Form("index_%i_tdc_%i_off", iter / 12, (iter % 12) / 2) << endl;
        }
    }
    bool ok = fitter.FitFCN();
    auto result = fitter.Result();
    result.Print(std::cout);

    //  current_pair (tdc, index, coarse, fine)
    for (auto current_pair : eligible_pairs)
    {
        //  fix "a manoni" to map everything from 1 to 769
        auto i_index = current_pair[0][1] < 140 ? current_pair[0][1] % 32 : current_pair[0][1] % 32 + 32;
        auto j_index = current_pair[1][1] < 140 ? current_pair[1][1] % 32 : current_pair[1][1] % 32 + 32;

        //  define a parameter index
        auto i_par_off_index = 1 + i_index * 4 * 3 + current_pair[0][0] * 3 + 0;
        auto i_par_slope_index = 1 + i_par_off_index;
        auto i_par_cut_index = 2 + i_par_off_index;
        auto j_par_off_index = 1 + j_index * 4 * 3 + current_pair[1][0] * 3 + 0;
        auto j_par_slope_index = 1 + j_par_off_index;
        auto j_par_cut_index = 2 + j_par_off_index;

        //  define the parameters involved
        auto i_off = result.Parameter(i_par_off_index);
        auto i_slope = result.Parameter(i_par_slope_index);
        auto i_cut = result.Parameter(i_par_cut_index);
        auto j_off = result.Parameter(j_par_off_index);
        auto j_slope = result.Parameter(j_par_slope_index);
        auto j_cut = result.Parameter(j_par_cut_index);

        //  Compute time
        auto i_time = current_pair[0][2] - (i_off + i_slope * current_pair[0][3] - 1 * (current_pair[0][3] > i_cut));
        auto j_time = current_pair[1][2] - (j_off + j_slope * current_pair[1][3] - 1 * (current_pair[1][3] > j_cut));

        //  Compute the delta of times
        uncorrected_coincidence->Fill(current_pair[0][2] - current_pair[1][2]);
        corrected_coincidence->Fill(i_time - j_time + result.Parameter(0));
        if ((current_pair[0][1] == 101) || (current_pair[0][0] == 0))
        {
            fine_test_p->Fill(current_pair[0][3], current_pair[0][2] - j_time + result.Parameter(0));
            fine_test->Fill(current_pair[0][3], i_time - j_time + result.Parameter(0));
        }
    }

    //  Plot
    TCanvas *c_coincidences = new TCanvas("c_coincidences", "", 1200, 600);
    gStyle->SetOptFit();
    c_coincidences->Divide(2, 1);
    c_coincidences->cd(1);
    uncorrected_coincidence->Draw("");
    uncorrected_coincidence->Fit("gaus");
    c_coincidences->cd(2);
    corrected_coincidence->Draw("");
    corrected_coincidence->Fit("gaus");

    TCanvas *c_fine_check = new TCanvas("c_fine_check", "", 1000, 1000);
    fine_test->Draw("COLZ");

    TCanvas *c_fine_check_2 = new TCanvas("c_fine_check_2", "", 1000, 1000);
    fine_test_p->Draw("");
    fine_test_p->Fit("pol1", "IRMES", "", 30, 97);
    TF1 *pol1_ = new TF1("pol1_", "[0]+[1]*x -1*(x>[2])", 0, 100);
    pol1_->SetLineColor(kGreen);
    pol1_->SetParameter(0, result.Parameter(61));
    pol1_->SetParameter(1, result.Parameter(62));
    pol1_->SetParameter(2, result.Parameter(63));
    pol1_->Draw("SAME");

    return;
}