#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/macros/recowriter.C"

void efficiency_calculation(std::string filename = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/lightdata.root") // Data/20231011-221500/ // Data/20240530-200252/
{
    //  Read data
    sipm4eic::lightio io;
    io.read_from_tree(filename);

    //  Plotting random walk
    TRandom *g_random = new TRandom();

    //  Output histograms
    //  --- Hitmaps
    //  --- --- Masks
    TH2F *h_deadmask = new TH2F("h_deadmask", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_partmask = new TH2F("h_partmask", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_goodpart = new TH2F("h_goodpart", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    //  --- --- Filtered hits
    TH2F *h_hitmap = new TH2F("h_hitmap", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_hitmap_hit_filter = new TH2F("h_hitmap_hit_filter", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_hitmap_hit_dcr = new TH2F("h_hitmap_hit_dcr", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_hitmap_zone_1 = new TH2F("h_hitmap_zone_1", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_hitmap_zone_2 = new TH2F("h_hitmap_zone_2", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_hitmap_zone_3 = new TH2F("h_hitmap_zone_3", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    //  --- Afterpulse
    TH1F *h_afterpulse = new TH1F("h_afterpulse", ";#Delta_{t}", 256, 0, 256);
    TH2F *h_afterpulse_ch = new TH2F("h_afterpulse_ch", ";channel;#Delta_{t}", 2500, 0, 2500, 256, 0, 256);
    TH1F *h_afterpulse_ch2 = new TH1F("h_afterpulse_ch2", ";channel;", 2500, 0, 2500);
    //  --- Reference time
    TH1F *h_delta_time = new TH1F("h_delta_time", ";channel;", 500, -250, 250);
    TH1F *h_delta_time_zone1 = new TH1F("h_delta_time_zone1", ";channel;", 500, -250, 250);
    TH1F *h_delta_time_zone2 = new TH1F("h_delta_time_zone2", ";channel;", 500, -250, 250);
    TH1F *h_delta_time_zone3 = new TH1F("h_delta_time_zone3", ";channel;", 500, -250, 250);

    auto n_frames_counter = 0;
    auto norm_zone_1 = 0;
    auto norm_zone_2 = 0;
    auto norm_zone_3 = 0;
    //  Loop over spills
    while (io.next_spill())
    {
        //  Get positions and available positions for spill
        sipm4eic::decode_good_positions(io);
        auto good_positions = io.get_good_positions();

        /*
// Get positions of participating SiPMs
        auto participating_map = io.get_good_fifos();
        for (auto [device, fifo] : participating_map)
            for (auto [i_fifo, active_fifo] : fifo)
                if (active_fifo)
                {
                    for (auto current_position : sipm4eic::get_fifo_position(device, i_fifo))
                    {
                        auto x_scatter = current_position[0] + 1.5 * g_random->Uniform(-1, 1);
                        auto y_scatter = current_position[1] + 1.5 * g_random->Uniform(-1, 1);
                        h_partmask->Fill(x_scatter, y_scatter);
                    }
                }
        participating_map = io.get_dead_fifos();
        for (auto [device, fifo] : participating_map)
            for (auto [i_fifo, active_fifo] : fifo)
                if (active_fifo)
                    for (auto current_position : sipm4eic::get_fifo_position(device, i_fifo))
                    {
                        auto x_scatter = current_position[0] + 1.5 * g_random->Uniform(-1, 1);
                        auto y_scatter = current_position[1] + 1.5 * g_random->Uniform(-1, 1);
                        h_deadmask->Fill(x_scatter, y_scatter);
                    }

        */

        //  Loop over frames
        while (io.next_frame())
        {
            //  Frames counter
            n_frames_counter++;

            // Trigger time reference
            auto trigger0_vector = io.get_trigger0_vector();
            auto reference_time = trigger0_vector[0].coarse;

            // loop over cherenkov hits
            auto cherenkov_map = io.get_cherenkov_map();
            for (auto &[index, hits] : cherenkov_map)
            {
                //  Sort hits
                std::sort(hits.begin(), hits.end());

                //  Loop on hits from the same channel
                if (hits.size() > 1)
                    for (auto i_hit = 0; i_hit < hits.size(); i_hit++)
                        for (auto j_hit = i_hit + 1; j_hit < hits.size(); j_hit++)
                        {
                            h_afterpulse->Fill(hits[j_hit].coarse - hits[i_hit].coarse);
                            h_afterpulse_ch->Fill(hits[i_hit].index + (hits[i_hit].device - 192) * 192, hits[j_hit].coarse - hits[i_hit].coarse);
                            if (hits[j_hit].coarse - hits[i_hit].coarse == 0)
                                h_afterpulse_ch2->Fill(hits[i_hit].index + (hits[i_hit].device - 192) * 192);
                        }
                //  Take first hit, syntax to take all hits
                for (auto i_hit = 0; i_hit < 1; /*hits.size();*/ i_hit++)
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
                    h_hitmap->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));

                    //  Fill expected hit hitmap
                    if (fabs(time_difference) < 10.)
                        h_hitmap_hit_filter->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));

                    if (fabs(fabs(time_difference) - 15) < 5.)
                        h_hitmap_hit_dcr->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));

                    //  Zones
                    //  - Zone 1: 0 < R < 60 in-time-bkg suspected
                    //  - Zone 2: 60 < R < 100 signal expected
                    //  - Zone 3: 100 < R in-time-bkg suppression suspected
                    auto current_hit_radius = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);

                    if (current_hit_radius < 60)
                        h_delta_time_zone1->Fill(time_difference);
                    if (current_hit_radius > 60 && current_hit_radius < 100)
                        h_delta_time_zone2->Fill(time_difference);
                    if (current_hit_radius > 100)
                        h_delta_time_zone3->Fill(time_difference);

                    if ((std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) >= 100) && (fabs(time_difference) < 10.))
                    {
                        h_hitmap_zone_1->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));
                        h_hitmap_zone_3->Fill(pos[0], pos[1], -1.);
                    }
                    if ((std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) >= 100) && (fabs(fabs(time_difference) - 15) < 5.))
                    {
                        h_hitmap_zone_2->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));
                        h_hitmap_zone_3->Fill(pos[0], pos[1], +1.);
                    }
                }
            }
        }
    }

    /*
    new TCanvas();
    h_partmask->Draw("colz");
    new TCanvas();
    h_deadmask->Draw("colz");
    */

    //  Poissonian + gaus?
    TF1 *fPoissonian = new TF1("fPoissonian", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 10);
    TF1 *fDoublePoissonian = new TF1("fDoublePoissonian", "fPoissonian(0)+fPoissonian(3)", 0, 10);
    fDoublePoissonian->SetParameters(1, 1, 1, 1, 1, 1);

    /*
    new TCanvas();
    h_afterpulse->Draw("");
    TF1 *fit_test = new TF1("fit_test", "pol3(0)+[4]*TMath::Gaus(x,[6],[5],1)+[7]*TMath::Gaus(x,[6]+[8],[5],1)+[9]*TMath::Gaus(x,[6]+2*[8],[5],1)+[10]*TMath::Gaus(x,[6]+3*[8],[5],1)");
    fit_test->SetParameters(300, 100, -3, 0.02, 75e3, 2, 10, 45e3, 10, 10e3, 1e3);
    h_afterpulse->Fit("fit_test", "IMEQ", "", 2, 45);
    new TCanvas();
    gStyle->SetOptStat(111111111);
    h_afterpulse_ch2->Draw("");

    new TCanvas();
    h_afterpulse_ch->Draw("COLZ");
*/

    TF1 *mixed_gaus_tail = new TF1("mixed_gaus_tail", "([0]+(2*(x<=0)-1)*[1]*x)+[4]*TMath::Gaus(x,[2],[3],1)*((x<[2]+[5])+(x>=[2]+[5])*TMath::Exp(-(x-[2]-[5])/[3]))+[8]*TMath::Gaus(x,[6],[7],1)", -1000, 1000); //+[4]*TMath::Gaus(x,[2],[3],kTrue)
    mixed_gaus_tail->SetParName(0, "offset");
    mixed_gaus_tail->SetParLimits(0, 0., 10.);
    mixed_gaus_tail->SetParName(1, "slope");
    mixed_gaus_tail->SetParLimits(1, 0., 2.);
    mixed_gaus_tail->SetParName(2, "gaus_mean");
    mixed_gaus_tail->SetParLimits(2, -2., 2.);
    mixed_gaus_tail->SetParName(3, "gaus_sigma");
    mixed_gaus_tail->SetParLimits(3, 0.3, 3.);
    mixed_gaus_tail->SetParName(4, "gaus_norm");
    mixed_gaus_tail->SetParLimits(4, 0., 50.);
    mixed_gaus_tail->SetParName(5, "gaus_exp_decay");
    mixed_gaus_tail->SetParLimits(5, 0., 10.);
    mixed_gaus_tail->SetParName(6, "gaus2_mean");
    mixed_gaus_tail->SetParLimits(6, 1., 6.);
    mixed_gaus_tail->SetParName(7, "gaus2_sigma");
    mixed_gaus_tail->SetParLimits(7, 0.3, 3.);
    mixed_gaus_tail->SetParName(8, "gaus2_norm");
    mixed_gaus_tail->SetParLimits(8, 0., 50.);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11111);

    new TCanvas();
    h_hitmap->Draw("COLZ");
    new TCanvas();
    h_hitmap_hit_filter->Draw("COLZ");
    new TCanvas();
    h_hitmap_hit_dcr->Draw("COLZ");

    new TCanvas();
    h_hitmap_zone_1->Draw("COLZ");
    new TCanvas();
    h_hitmap_zone_2->Draw("COLZ");
    new TCanvas();
    h_hitmap_zone_3->Draw("COLZ");

    TCanvas *c1 = new TCanvas("", "", 1000, 400);
    c1->Divide(3, 1);
    c1->cd(1);
    gPad->SetLogy();
    h_delta_time_zone1->Draw("");
    h_delta_time_zone1->Scale(1. / n_frames_counter);
    h_delta_time_zone1->GetXaxis()->SetRangeUser(-50., 50.);
    h_delta_time_zone1->Fit("mixed_gaus_tail", "");
    cout << "Zone 1: " << mixed_gaus_tail->GetParameter(2) << " +/- " << mixed_gaus_tail->GetParError(2) << endl;
    c1->cd(2);
    gPad->SetLogy();
    h_delta_time_zone2->Draw("");
    h_delta_time_zone2->Scale(1. / n_frames_counter);
    h_delta_time_zone2->GetXaxis()->SetRangeUser(-50., 50.);
    h_delta_time_zone2->Fit("mixed_gaus_tail", "");
    cout << "Zone 2: " << mixed_gaus_tail->GetParameter(2) << " +/- " << mixed_gaus_tail->GetParError(2) << endl;
    c1->cd(3);
    gPad->SetLogy();
    h_delta_time_zone3->Draw("");
    h_delta_time_zone3->Scale(1. / n_frames_counter);
    h_delta_time_zone3->GetXaxis()->SetRangeUser(-50., 50.);
    h_delta_time_zone3->Fit("mixed_gaus_tail", "");
    cout << "Zone 3: " << mixed_gaus_tail->GetParameter(2) << " +/- " << mixed_gaus_tail->GetParError(2) << endl;

    new TCanvas();
    h_delta_time->Draw("COLZ");

    TFile *fout = new TFile("fout_data.root", "RECREATE");
    h_afterpulse->Write();
    h_afterpulse_ch->Write();
    h_afterpulse_ch2->Write();
}