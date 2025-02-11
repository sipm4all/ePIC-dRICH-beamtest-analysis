#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/macros/recowriter.C"

void efficiency_calculation(std::string filename = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/Data/20231011-221500/lightdata.root")
{
    //  Read data
    sipm4eic::lightio io;
    io.read_from_tree(filename);

    //
    TRandom *g_random = new TRandom();

    //  Output histograms
    TH2F *h_deadmask = new TH2F("h_deadmask", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_partmask = new TH2F("h_partmask", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);
    TH2F *h_goodpart = new TH2F("h_goodpart", ";x (mm);y (mm)", 594, -99, 99, 594, -99, 99);

    //  Loop over spills
    while (io.next_spill())
    {
        auto participating_map = io.get_active_fifos();
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
        //  Loop over frames
        /*
        while (io.next_frame())
        {
            continue;
        }
        */
    }

    new TCanvas();
    h_partmask->Draw("colz");
    new TCanvas();
    h_deadmask->Draw("colz");
}