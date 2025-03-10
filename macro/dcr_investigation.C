#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/macros/recowriter.C"

void dcr_investigation(std::string filename = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/lightdatadcr.root")
{
    //  Read data
    sipm4eic::lightio io;
    io.read_from_tree(filename);

    //  Output histograms
    TH2F *h_nolight_check = new TH2F("h_nolight_check", ";x (mm);y (mm)", 1980, -99, 99, 1980, -99, 99);
    TH2F *h_highfrq_check = new TH2F("h_highfrq_check", ";x (mm);y (mm)", 1980, -99, 99, 1980, -99, 99);
    TGraphErrors *g_frequencies = new TGraphErrors();
    g_frequencies->SetName("g_frequencies");
    std::map<std::array<float, 2>, TH2F *> h_dcr_measurement;

    TRandom *g_random = new TRandom();

    //  Loop over spills
    while (io.next_spill())
    {
        cout << "[INFO] current spill: " << io.get_spill_current() << endl;

        //  Get positions and available positions for spill
        sipm4eic::decode_good_positions(io);
        auto good_positions = io.get_good_positions();

        //  Initialise histograms
        if (io.get_spill_current() == 1)
            for (auto [current_position, state] : good_positions)
                h_dcr_measurement[current_position] = new TH2F(Form("h_%.2f_%.2f", current_position[0], current_position[1]), Form("h_%.2f_%.2f;N_{#gamma};N_{spill}", current_position[0], current_position[1]), 10, 0, 10, 99, 1, 100);

        //  Loop over frames
        while (io.next_frame())
        {

            //  Initialise DCR counter
            std::map<std::array<float, 2>, int> dcr_counters;

            //  Get trigger vector
            auto trigger0_vector = io.get_trigger0_vector();
            auto ref = g_random->Uniform(0, io.get_frame_size());
            if (trigger0_vector.size() != 0)
                ref = trigger0_vector[0].coarse;

            //  Avoid border effects
            if (ref <= 25 || ref >= (io.get_frame_size() - 25))
                continue;

            //  Get cherenkov vector
            auto cherenkov_vector = io.get_cherenkov_vector();
            for (auto &cherenkov : cherenkov_vector)
            {
                auto delta = ref - cherenkov.coarse;
                if (fabs(delta) <= 20)
                    continue;
                auto pos = sipm4eic::get_position(sipm4eic::get_geo(cherenkov));
                dcr_counters[pos]++;
                h_nolight_check->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));
            }

            //  Fill dcr counter plots
            for (auto [pos, val] : dcr_counters)
                if (good_positions[pos])
                    h_dcr_measurement[pos]->Fill(val, io.get_spill_current());
        }

        //  Set all zero entries
        for (auto [pos, val] : h_dcr_measurement)
        {
            if (!good_positions[pos])
                continue;
            auto spill_projection = val->ProjectionX("tmp", io.get_spill_current(), io.get_spill_current());
            auto zero_entries = io.get_frame_current() - spill_projection->GetEntries();
            val->SetBinContent(1, io.get_spill_current(), zero_entries);
        }
    }

    //  Fit all DCR plots

    //  Fit
    std::map<std::array<float, 2>, TGraphErrors *> dcr_graphs;
    for (auto [pos, val] : h_dcr_measurement)
    {
        dcr_graphs[pos] = new TGraphErrors();
        dcr_graphs[pos]->SetMarkerStyle(20);
        dcr_graphs[pos]->SetMarkerColor(kBlack);

        TF1 *fPhotonFitFunction = new TF1("hPhotonFitFunction", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", -1000, 1000);
        fPhotonFitFunction->FixParameter(2, 1);

        auto avg_val = 0.;
        auto avg_err = 0.;
        auto eff_spills = 0;
        for (auto current_spill = 1; current_spill < 100; current_spill++)
        {
            fPhotonFitFunction->SetParameters(1, 1, 1);
            auto current_point = dcr_graphs[pos]->GetN();
            auto spill_projection = val->ProjectionX("tmp", current_spill, current_spill);
            if (spill_projection->GetEntries() == 0)
                continue;
            spill_projection->Fit(fPhotonFitFunction, "Q");
            auto norm_dcr = 1. / (io.get_frame_size() * sipm4eic::lightdata::coarse_to_ns * 1.e-9);
            dcr_graphs[pos]->SetPoint(current_point, current_spill, fPhotonFitFunction->GetParameter(1) * norm_dcr);
            dcr_graphs[pos]->SetPointError(current_point, 0, fPhotonFitFunction->GetParError(1) * norm_dcr);
            avg_val += fPhotonFitFunction->GetParameter(1) * norm_dcr;
            avg_err += fPhotonFitFunction->GetParError(1) * norm_dcr * fPhotonFitFunction->GetParError(1) * norm_dcr;
            eff_spills++;
        }
        if (eff_spills == 0)
            continue;
        avg_val /= eff_spills;
        avg_err = TMath::Sqrt(avg_err);
        avg_err /= eff_spills;
        auto current_point = g_frequencies->GetN();
        g_frequencies->SetPoint(current_point, pos[0], pos[1]);
        g_frequencies->SetPointError(current_point, avg_val, avg_err);

        for (auto iPnt = 0; iPnt < avg_val; iPnt++)
            h_highfrq_check->Fill(pos[0] + 1.5 * g_random->Uniform(-1, 1), pos[1] + 1.5 * g_random->Uniform(-1, 1));

    }

    //  Canvas
    new TCanvas();
    TH1F *frame = new TH1F("frame", ";spill;DCR (Hz)", 100, 0, 100);
    frame->SetMaximum(1000000000);
    frame->SetMinimum(1);
    frame->Draw();
    gPad->SetLogy();
    for (auto [pos, graph] : dcr_graphs)
        graph->Draw("SAME EP");

    new TCanvas();
    h_nolight_check->Draw("COLZ");

    new TCanvas();
    h_highfrq_check->Draw("COLZ");


    TFile *outfile = new TFile("dcr_test_dcr.root", "RECREATE");
    g_frequencies->Write();
    for (auto [key, val] : h_dcr_measurement)
        val->Write();
    for (auto [key, val] : dcr_graphs)
        val->Write();
    outfile->Close();
}
