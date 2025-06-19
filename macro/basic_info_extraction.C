#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"

void basic_info_extraction(std::string filename = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/recodata.filter_rollover.root") // fastmc_full
{
    //  Load data from run
    testbeam::data current_data(filename);
    auto current_tree = current_data.get_tree();
    current_data.set_up_global_variables();
    std::string path = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/";

    //  Event selection
    TH1F *hEventCounting = new TH1F("hEventCounting", "", 10, 0, 10);
    hEventCounting->GetXaxis()->SetBinLabel(1, "All evs.");
    hEventCounting->GetXaxis()->SetBinLabel(2, "Evs. w/ less than 4 hits");
    hEventCounting->GetXaxis()->SetBinLabel(3, "Evs. w/o clusters");
    hEventCounting->GetXaxis()->SetBinLabel(4, "Evs. w/ bigger cluster w/ less than 3 hits");
    hEventCounting->GetXaxis()->SetBinLabel(5, "Evs. w/ radius outside range");
    hEventCounting->GetXaxis()->SetBinLabel(6, "Evs. w/ bigger clusters w/ more than 20 hits");
    hEventCounting->GetXaxis()->SetBinLabel(7, "Evs. w/o at least 3 hits in radius range");
    hEventCounting->GetXaxis()->SetBinLabel(8, "Accepted events");
    TH2F *hEvntEffMap0 = new TH2F("hEvntEffMap0", "Evs. w/ less than 4 hits", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap1 = new TH2F("hEvntEffMap1", "Evs. w/o clusters", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap2 = new TH2F("hEvntEffMap2", "Evs. w/ bigger cluster w/ less than 3 hits", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap3 = new TH2F("hEvntEffMap3", "Evs. w/ radius outside range", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap4 = new TH2F("hEvntEffMap4", "Evs. w/ bigger clusters w/ more than 20 hits", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap5 = new TH2F("hEvntEffMap5", "Evs. w/o at least 3 hits in radius range", 100, -100, 100, 100, -100, 100);
    TH2F *hEvntEffMap6 = new TH2F("hEvntEffMap6", "Accepted events", 100, -100, 100, 100, -100, 100);

    //  Radius distribution
    TH1D *h_radius_distribution_all = new TH1D("h_radius_distribution_all", "", 80, -60, 60);
    TH1F *h_radius_distribution_sig = new TH1F("h_radius_distribution_sig", "", 80, -60, 60);
    TH1F *h_radius_distribution_bkg = new TH1F("h_radius_distribution_bkg", "", 80, -60, 60);
    TH1F *h_radius_distribution_dcr = new TH1F("h_radius_distribution_dcr", "", 80, -60, 60);
    TH1F *h_fit_radius = new TH1F("h_fit_radius", "", 1000, 0, 100);
    TH2F *h3_fit_radius = new TH2F("h3_fit_radius", "", 100, -10, 10, 40, 0, 40);
    TH2F *h2_fit_radius = new TH2F("h2_fit_radius", "", 1000, -50, 50, 40, 0, 40);

    //  Photon counting
    TH1F *h_photon_counting_DCR = new TH1F("h_photon_counting_DCR", ";N_{#gamma}", 100, 0, 100);
    TH1F *h_photon_counting_bkg = new TH1F("h_photon_counting_bkg", ";N_{#gamma}", 100, 0, 100);
    TH1F *h_photon_counting_all = new TH1F("h_photon_counting_all", ";N_{#gamma}", 100, 0, 100);
    TH2F *h_photon_counting_bkg_all = new TH2F("h_photon_counting_bkg_all", ";N_{#gamma} (bkg);N_{#gamma} (sig)", 100, 0, 100, 100, 0, 100);

    //  Correlations
    TH2F *h_delta_R_Phi = new TH2F("h_delta_R_Phi", ";#Delta_{#phi};#Delta_{R}", 45, -45, 135, 50, -100, 100);

    for (auto iev = 0; iev < current_tree->GetEntries(); iev++)
    {
        //  Load event
        current_tree->GetEntry(iev);

        //  Reject unpopulated frames (w/ less than overall 4 hits)
        if (current_data.get_n() < 4)
        {
            for (auto i = 0; i <= 1; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hEvntEffMap0->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        //  Reject frames w/o clusters
        auto cluster_points = current_data.get_r_clusters(5);
        if (!cluster_points.size())
        {
            for (auto i = 0; i <= 2; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hEvntEffMap1->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        //  Reject frames w/ bigger clusters w/ less than 3 hits
        if (cluster_points[0].size() < 3)
        {
            for (auto i = 0; i <= 3; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hEvntEffMap2->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        //  Reject frames w/ bigger clusters w/ more than 20 hits
        auto bkg_cluster_points = current_data.get_bkg_clusters();
        if (bkg_cluster_points.size())
            if (bkg_cluster_points[0].size() > 8)
            {
                for (auto i = 0; i <= 4; i++)
                    hEventCounting->Fill(i);
                for (auto i = 0; i < current_data.get_n(); i++)
                    hEvntEffMap3->Fill(current_data.get_x(i), current_data.get_y(i));
                continue;
            }

        //  Reject events w/o at least 3 hits in radius range
        auto good_points_array = current_data.select_points({current_data.get_common_center_x(), current_data.get_common_center_y(), current_data.get_common_radius(), current_data.get_timing_center()}, {999, 999, 1.8 * 3, 4 * current_data.get_timing_sigma()});
        auto good_points = good_points_array[0];
        if (good_points.size() < 3)
        {
            for (auto i = 0; i <= 5; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hEvntEffMap4->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        //  Fit selected points
        auto current_fit_results = current_data.fit_circle(good_points, {0., 0., 0.}, true);
        if (fabs(current_fit_results[2][0] - current_data.get_common_radius()) > 1.8 * 3)
        {
            for (auto i = 0; i <= 6; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hEvntEffMap5->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }
        for (auto i = 0; i <= 7; i++)
            hEventCounting->Fill(i);
        for (auto i = 0; i < current_data.get_n(); i++)
            hEvntEffMap6->Fill(current_data.get_x(i), current_data.get_y(i));

        //  Correlation test
        for (auto ihit = 0; ihit < current_data.get_n(); ihit++)
            for (auto jhit = ihit + 1; jhit < current_data.get_n(); jhit++)
            {
                if (ihit == jhit)
                    continue;
                h_delta_R_Phi->Fill(current_data.get_phi_delta(ihit, jhit), current_data.get_r_delta(ihit, jhit), 0.5*std::min(100.,1.*std::max(0.01,1.*fabs(current_data.get_t_delta(ihit, jhit)))));
            }

        //  Fill radius
        h_fit_radius->Fill(current_fit_results[2][0]);
        h2_fit_radius->Fill(current_fit_results[2][0] - current_data.get_common_radius(), good_points.size());
        auto diff_circle_vals = current_data.diff_circle(good_points);
        for (auto current_diff : diff_circle_vals)
            h3_fit_radius->Fill(current_diff, good_points.size());
        h_photon_counting_bkg->Fill(current_data.get_n() - good_points.size());
        h_photon_counting_all->Fill(good_points.size());
        auto count_dcr = 0;
        for (auto ihit = 0; ihit < current_data.get_n(); ihit++)
        {
            h_radius_distribution_all->Fill(current_data.get_radius(ihit) - current_data.get_common_radius());
            if (current_data.is_dcr(ihit))
            {
                h_radius_distribution_dcr->Fill(current_data.get_radius(ihit) - current_data.get_common_radius());
                count_dcr++;
            }
        }
        for (auto ihit : good_points)
            h_radius_distribution_sig->Fill(current_data.get_radius(ihit));
        h_photon_counting_DCR->Fill(count_dcr);
        h_photon_counting_bkg_all->Fill(current_data.get_n() - good_points.size(), good_points.size());
    }

    auto _graph = testbeam::measure_resolution_in_mult_photons(h2_fit_radius, 1000);
    TGraphErrors *graph = new TGraphErrors();
    for (auto current_point : _graph)
    {
        auto ipnt = graph->GetN();
        graph->SetPoint(ipnt, current_point[0][0], current_point[1][0]);
        graph->SetPointError(ipnt, current_point[0][1], current_point[1][1]);
    }

    TF1 *fitfunc = new TF1("fitfunc", "[0]+[1]/(TMath::Sqrt(x))");
    fitfunc->SetParameter(0, 0);
    fitfunc->FixParameter(0, 0);

    TF1 *fPhotonFitFunction = new TF1("hPhotonFitFunction", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", -1000, 1000);
    fPhotonFitFunction->SetParameters(1, 1, 1);
    fPhotonFitFunction->FixParameter(2, 1);

    new TCanvas();
    graph->Fit(fitfunc);
    graph->Draw("ALP");

    std::ofstream data_stream(path + "/out.txt");
    data_stream << "radius sigma " << fitfunc->GetParameter(1) << "\t" << fitfunc->GetParError(1) << "\n";

    //  norm dcr
    auto norm_all_n = 0, norm_bkg_n = 0;
    for (auto [key, value] : current_data.get_available_SiPMs())
    {
        auto current_x = key[0] - current_data.get_common_center_x();
        auto current_y = key[1] - current_data.get_common_center_y();
        auto current_r = current_data.get_common_radius() - TMath::Sqrt(current_x * current_x + current_y * current_y);
        if (fabs(current_r) < 3)
            norm_all_n++;
        else
            norm_bkg_n++;
    }

    auto norm_dcr = (1. / (9 * current_data.get_available_SiPMs().size())) * (1. / (1.e-9 * (80. - 4 * current_data.get_timing_sigma())));
    auto norm_all = (1. / (9 * norm_all_n)) * (1. / (1.e-9 * (4 * current_data.get_timing_sigma())));
    auto norm_bkg = (1. / (9 * norm_bkg_n)) * (1. / (1.e-9 * (4 * current_data.get_timing_sigma())));

    new TCanvas();
    h_photon_counting_DCR->Draw("SAME");
    h_photon_counting_DCR->Fit(fPhotonFitFunction, "Q");
    data_stream << "n_gamma dcr " << fPhotonFitFunction->GetParameter(1) * norm_dcr << "\t" << fPhotonFitFunction->GetParError(1) * norm_dcr << "\n";
    h_photon_counting_bkg->Draw("SAME");
    h_photon_counting_bkg->Fit(fPhotonFitFunction, "Q");
    data_stream << "n_gamma bkg " << fPhotonFitFunction->GetParameter(1) * norm_bkg << "\t" << fPhotonFitFunction->GetParError(1) * norm_bkg << "\n";
    h_photon_counting_all->Draw("SAME");
    h_photon_counting_all->Fit(fPhotonFitFunction, "Q");
    data_stream << "n_gamma all " << fPhotonFitFunction->GetParameter(1) * norm_all << "\t" << fPhotonFitFunction->GetParError(1) * norm_all << "\n";
    data_stream << "center_x " << current_data.get_common_center_x() << "\t" << "\n";
    data_stream << "center_y " << current_data.get_common_center_y() << "\t" << "\n";
    data_stream << "radius " << current_data.get_common_radius() << "\t" << "\n";

    new TCanvas();
    h_photon_counting_bkg_all->Draw("COLZ");
    new TCanvas();

    TGraphErrors *g_photon_counting_bkg = new TGraphErrors();
    g_photon_counting_bkg->SetName("g_photon_counting_bkg");
    g_photon_counting_bkg->SetTitle(";N_{#gamma} (sig);N_{#gamma} (bkg)");
    for (auto ybin = 1; ybin <= h_photon_counting_bkg_all->GetNbinsY(); ybin++)
    {
        //  Get current slice under test
        auto current_point = g_photon_counting_bkg->GetN();
        auto current_slice = h_photon_counting_bkg_all->ProjectionX(Form("proj_%i", ybin), ybin, ybin);
        if (current_slice->GetEntries() < 100)
            continue;

        auto y_low_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(ybin);
        auto y_high_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(ybin + 1);
        current_slice->Fit(fPhotonFitFunction, "I");

        g_photon_counting_bkg->SetPoint(current_point, 0.5 * (y_high_edge + y_low_edge) - 0.5, fPhotonFitFunction->GetParameter(1));
        g_photon_counting_bkg->SetPointError(current_point, 0.5 * (y_high_edge - y_low_edge), fPhotonFitFunction->GetParError(1));
    }

    TGraphErrors *g_photon_counting_sig = new TGraphErrors();
    g_photon_counting_sig->SetName("g_photon_counting_sig");
    g_photon_counting_sig->SetTitle(";N_{#gamma} (bkg);N_{#gamma} (sig)");
    for (auto xbin = 1; xbin <= h_photon_counting_bkg_all->GetNbinsX(); xbin++)
    {
        //  Get current slice under test
        auto current_point = g_photon_counting_sig->GetN();
        auto current_slice = h_photon_counting_bkg_all->ProjectionY(Form("proj_%i", xbin), xbin, xbin);
        if (current_slice->GetEntries() < 100)
            continue;

        auto x_low_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(xbin);
        auto x_high_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(xbin + 1);
        current_slice->Fit(fPhotonFitFunction, "I");

        g_photon_counting_sig->SetPoint(current_point, 0.5 * (x_high_edge + x_low_edge) - 0.5, fPhotonFitFunction->GetParameter(1));
        g_photon_counting_sig->SetPointError(current_point, 0.5 * (x_high_edge - x_low_edge), fPhotonFitFunction->GetParError(1));
    }

    new TCanvas();
    g_photon_counting_bkg->Draw("ALPE");
    new TCanvas();
    g_photon_counting_sig->Draw("ALPE");

    new TCanvas();
    h_delta_R_Phi->Draw("COLZ");

    new TCanvas();

    TFile *outfile = new TFile((path + "/preliminary_test.root").c_str(), "RECREATE");
    hEventCounting->Write();
    hEvntEffMap0->Write();
    hEvntEffMap1->Write();
    hEvntEffMap2->Write();
    hEvntEffMap3->Write();
    hEvntEffMap4->Write();
    hEvntEffMap5->Write();
    hEvntEffMap6->Write();

    h_photon_counting_bkg->Write();
    h_photon_counting_all->Write();
    h_photon_counting_DCR->Write();
    h_photon_counting_bkg_all->Write();

    g_photon_counting_sig->Write();
    g_photon_counting_bkg->Write();

    h_radius_distribution_all->Write();
    h_radius_distribution_dcr->Write();

    graph->Write();

    h_delta_R_Phi->Scale(1. / current_tree->GetEntries());
    h_delta_R_Phi->Write();

    outfile->Close();

    h_fit_radius->Draw();
    h_fit_radius->Fit("gaus");
    new TCanvas();
    h2_fit_radius->Draw();
    h2_fit_radius->Fit("gaus");
    new TCanvas();
    h3_fit_radius->Draw();
    h3_fit_radius->Fit("gaus");

    TCanvas *c_ev_selection = new TCanvas("c_ev_selection", "c_ev_selection", 1000, 500);
    c_ev_selection->Divide(2, 1);
    c_ev_selection->cd(1);
    hEventCounting->Draw("");
    c_ev_selection->cd(2);
    hEvntEffMap6->Draw("COLZ");

    TCanvas *c_ev_selection_map = new TCanvas("c_ev_selection_map", "c_ev_selection_map", 1000, 1000);
    c_ev_selection_map->Divide(3, 3);
    c_ev_selection_map->cd(1);
    hEvntEffMap0->Draw("COLZ");
    c_ev_selection_map->cd(2);
    hEvntEffMap1->Draw("COLZ");
    c_ev_selection_map->cd(3);
    hEvntEffMap2->Draw("COLZ");
    c_ev_selection_map->cd(4);
    hEvntEffMap3->Draw("COLZ");
    c_ev_selection_map->cd(5);
    hEvntEffMap4->Draw("COLZ");
    c_ev_selection_map->cd(6);
    hEvntEffMap5->Draw("COLZ");

    TF1 *f_radius_dist = new TF1("f_radius_dist", "[0]*TMath::Gaus(x,[1],[2],kTRUE)+[3]*TMath::Gaus(x,[4],[5],kTRUE)+[6]+[7]*x+[8]*x*x+[9]*x*x*x+[10]*x*x*x*x");
    f_radius_dist->SetParameters(1, 1, 1, 0.33, 10, 4, 0.1, -0.001, 0, 0);
    f_radius_dist->FixParameter(3, 0);
    f_radius_dist->FixParameter(4, 0);
    f_radius_dist->FixParameter(5, 0);
    f_radius_dist->SetParName(0, "Norm_gaus_1");
    f_radius_dist->SetParName(1, "Mean_gaus_1");
    f_radius_dist->SetParName(2, "Sigma_gaus_1");
    f_radius_dist->SetParName(3, "Norm_gaus_2");
    f_radius_dist->SetParName(4, "Mean_gaus_2");
    f_radius_dist->SetParName(5, "Sigma_gaus_2");
    f_radius_dist->SetParName(6, "pol4_par0");
    f_radius_dist->SetParName(7, "pol4_par1");
    f_radius_dist->SetParName(8, "pol4_par2");
    f_radius_dist->SetParName(9, "pol4_par3");
    f_radius_dist->SetParName(10, "pol4_par4");

    testbeam::fit_ring_peak(h_radius_distribution_all, -1);

    TCanvas *c_radius_resolution = new TCanvas("c_radius_resolution", "c_radius_resolution", 1000, 1000);
    gPad->SetLogy();
    gStyle->SetOptFit(1);
    h_radius_distribution_all->Scale(1. / hEventCounting->GetBinContent(8));
    h_radius_distribution_all->Scale(1., "width");
    h_radius_distribution_all->SetLineColor(kBlack);
    h_radius_distribution_all->Draw("SAME");
    h_radius_distribution_all->Fit("f_radius_dist");
}