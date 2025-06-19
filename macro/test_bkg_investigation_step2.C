#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/macro/test_bkg_investigation.C"
#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/Utility.h"

TF1 *pol2_bkg = new TF1("pol2_bkg", "pol2(0)", 0, 150);
//  Runlists
//  --- Std first run
//  std::string run = "20240526-220929"; // 11.5 GeV
//  --- No radiators
//  std::string run = "20240530-200252"; // 11 GeV
//  --- Momscan Pos Aerogel
//  std::string run = "20240528-052051"; // 2 GeV
//  std::string run = "20240528-011104"; // 11 GeV
//  --- Filter scan
std::vector<int> CreateWavelengthColorPalette()
{
    std::vector<int> palette;

    // Define RGB values for each wavelength interval (adjust as needed)
    std::vector<std::tuple<double, double, double>> colors = {
        {0.5, 0.0, 0.5}, // <300 nm — violet
        {0.4, 0.0, 0.8}, // 280–320 nm — near UV (bluish-purple)
        {0.0, 0.0, 1.0}, // 400–450 nm — blue
        {0.0, 1.0, 1.0}, // 450–500 nm — cyan
        {0.0, 1.0, 0.0}, // 500–550 nm — green
        {1.0, 1.0, 0.0}  // 550–600 nm — yellow
    };

    for (const auto &[r, g, b] : colors)
    {
        int colorIndex = TColor::GetFreeColorIndex();
        new TColor(colorIndex, r, g, b);
        palette.push_back(colorIndex);
    }

    return palette;
}

std::string run_bkg = "20240530-200252"; // 11 GeV

std::vector<std::pair<std::array<float, 2>, std::string>> run_list_filter_scan = {
    {{0, 300}, "20240528-181849"},   // <300nm
    {{280, 320}, "20240530-233229"}, // 280-320nm
    {{400, 450}, "20240529-125500"}, // 400-450nm
    {{450, 500}, "20240529-145727"}, // 450-500nm
    {{500, 550}, "20240529-173251"}, // 500-550nm
    {{550, 600}, "20240529-211829"}, // 550-600nm
    //{{600, 999}, "20240531-163657"} // >600nm ***
};

TF1 *gaus_pol2 = new TF1("gaus_pol2", "gaus(0)+pol2(3)", 0, 150);
gaus_pol2->SetParNames("Amplitude", "Mean", "Width", "p0", "p1", "p2");

void test_bkg_investigation_step2()
{
    //  Load sensors PDE
    TFile *fin = new TFile("Data/PDE.root");
    TGraph *pde_13_50 = (TGraph *)(fin->Get("13_50_PDE"));
    TGraph *pde_13_75 = (TGraph *)(fin->Get("13_75_PDE"));
    TGraph *pde_14_50 = (TGraph *)(fin->Get("14_50_PDE"));

    /*
    auto entries = 10000000;
    TH1F *production_probability_13_50 = new TH1F("production_probability_13_50", "", 500, 0, 1000);
    TH1F *production_probability_13_75 = new TH1F("production_probability_13_75", "", 500, 0, 1000);
    std::array<double, 2> cherenkov_emission_wavelength = {280, 850};
    for (auto i = 0; i < entries; i++)
    {
        double emission_wavelength = (cherenkov_emission_wavelength[0]) / (1 - (gRandom->Uniform(0, 1)) * (cherenkov_emission_wavelength[1] - cherenkov_emission_wavelength[0]) / (cherenkov_emission_wavelength[1]));
        auto current_PDE = eval_with_errors(pde_13_50, emission_wavelength);
        if (current_PDE[0] > gRandom->Uniform(0, 100))
            production_probability_13_50->Fill(emission_wavelength);
        current_PDE = eval_with_errors(pde_13_75, emission_wavelength);
        if (current_PDE[0] > gRandom->Uniform(0, 100))
            production_probability_13_75->Fill(emission_wavelength);
    }

    production_probability_13_50->Scale(1. / entries);
    production_probability_13_75->Scale(1. / entries);
    auto production_probability_ratio = new TH1F("production_probability_ratio", "", 500, 0, 1000);
    auto production_probability_diff = new TH1F("production_probability_diff", "", 500, 0, 1000);
    production_probability_ratio->Divide(production_probability_13_75, production_probability_13_50, 1., 1., "B");
    production_probability_diff->Add(production_probability_13_75, production_probability_13_50, 1., -1.);

    new TCanvas("c0", "c0", 800, 800);
    production_probability_13_75->Draw("SAME EP");
    production_probability_13_50->Draw("SAME EP");

    new TCanvas("c1", "c0", 800, 800);
    production_probability_ratio->Draw("SAME EP");

    new TCanvas("c2", "c0", 800, 800);
    production_probability_diff->Draw("SAME EP");

    cout << "Production probability (13_50): " << production_probability_13_50->Integral(production_probability_13_50->GetXaxis()->FindBin(280), production_probability_13_50->GetXaxis()->FindBin(320)) << endl;
    cout << "Production probability (13_75): " << production_probability_13_75->Integral(production_probability_13_75->GetXaxis()->FindBin(280), production_probability_13_75->GetXaxis()->FindBin(320)) << endl;
    cout << "Production probability (13_50): " << production_probability_13_50->Integral(production_probability_13_50->GetXaxis()->FindBin(0.), production_probability_13_50->GetXaxis()->FindBin(300)) << endl;
    cout << "Production probability (13_75): " << production_probability_13_75->Integral(production_probability_13_75->GetXaxis()->FindBin(0.), production_probability_13_75->GetXaxis()->FindBin(300)) << endl;

    return;
    */

    for (auto i_pnt = 0; i_pnt < pde_13_50->GetN(); i_pnt++)
    {
        auto _x = pde_13_50->GetPointX(i_pnt);
        auto _y = pde_13_50->GetPointY(i_pnt);
        pde_13_50->SetPoint(i_pnt, _x, (0.25 / 0.415) * _y / 100.);
    }
    for (auto i_pnt = 0; i_pnt < pde_13_75->GetN(); i_pnt++)
    {
        auto _x = pde_13_75->GetPointX(i_pnt);
        auto _y = pde_13_75->GetPointY(i_pnt);
        pde_13_75->SetPoint(i_pnt, _x, (0.25 / 0.415) * _y / 100.);
    }

    // Set Fit function
    gaus_pol2->SetParLimits(0, 0., 1.);
    gaus_pol2->SetParLimits(1, 70, 80.);
    gaus_pol2->SetParLimits(2, 0., 6.);

    auto palette = CreateWavelengthColorPalette();
    // test_bkg_investigation(run_bkg);
    TFile *bkg_file = new TFile(Form("/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/test_bkg_investigation.%s.root", run_bkg.c_str()));
    auto bkg_hist_sensor = (TH1F *)bkg_file->Get("h_signal_sensor_radius");
    auto bkg_hist_sensor0 = (TH1F *)bkg_file->Get("h_signal_sensor0_radius");
    auto bkg_hist_sensor1 = (TH1F *)bkg_file->Get("h_signal_sensor1_radius");
    TCanvas *c0 = new TCanvas("c0", "c0", 800, 800);
    bkg_hist_sensor->Draw("SAME EP");
    bkg_hist_sensor->SetMarkerStyle(20);
    bkg_hist_sensor0->Draw("SAME HIST");
    bkg_hist_sensor0->SetMarkerStyle(0);
    bkg_hist_sensor0->SetLineColor(kGreen - 2);
    bkg_hist_sensor1->Draw("SAME HIST");
    bkg_hist_sensor1->SetMarkerStyle(0);
    bkg_hist_sensor1->SetLineColor(kRed - 1);
    bkg->SetParameters(2., -1., 100., 12.);
    bkg_hist_sensor->Fit(bkg, "QIM", "", 30, 115);

    //  TGraphErrors
    //  --- Average radius
    TGraphErrors *g_average_radius = new TGraphErrors();
    g_average_radius->SetName("g_average_radius");
    g_average_radius->SetTitle(";Wavelength (nm);R (mm)");
    TGraphErrors *g_average_radius_sensor0 = new TGraphErrors();
    g_average_radius_sensor0->SetName("g_average_radius_sensor0");
    g_average_radius_sensor0->SetTitle(";Wavelength (nm);R (mm)");
    TGraphErrors *g_average_radius_sensor1 = new TGraphErrors();
    g_average_radius_sensor1->SetName("g_average_radius_sensor1");
    g_average_radius_sensor1->SetTitle(";Wavelength (nm);R (mm)");

    //  --- Average yield
    TGraphErrors *g_average_yield = new TGraphErrors();
    g_average_yield->SetName("g_average_yield");
    g_average_yield->SetTitle(";Wavelength (nm);Yield");
    TGraphErrors *g_average_yield_sensor0 = new TGraphErrors();
    g_average_yield_sensor0->SetName("g_average_yield_sensor0");
    g_average_yield_sensor0->SetTitle(";Wavelength (nm);Yield");
    TGraphErrors *g_average_yield_sensor1 = new TGraphErrors();
    g_average_yield_sensor1->SetName("g_average_yield_sensor1");
    g_average_yield_sensor1->SetTitle(";Wavelength (nm);Yield");

    //  --- Average resolution
    TGraphErrors *g_average_resolution = new TGraphErrors();
    g_average_resolution->SetName("g_average_resolution");
    g_average_resolution->SetTitle(";Wavelength (nm);#sigma_{R} (mm)");
    TGraphErrors *g_average_resolution_sensor0 = new TGraphErrors();
    g_average_resolution_sensor0->SetName("g_average_resolution_sensor0");
    g_average_resolution_sensor0->SetTitle(";Wavelength (nm);#sigma_{R} (mm)");
    TGraphErrors *g_average_resolution_sensor1 = new TGraphErrors();
    g_average_resolution_sensor1->SetName("g_average_resolution_sensor1");
    g_average_resolution_sensor1->SetTitle(";Wavelength (nm);#sigma_{R} (mm)");

    //  --- additional bkg
    TGraphErrors *g_average_yield_bkg = new TGraphErrors();
    g_average_yield_bkg->SetName("g_average_yield_bkg");
    g_average_yield_bkg->SetTitle(";Wavelength (nm);Bkg Yield");
    TGraphErrors *g_average_yield_bkg_sensor0 = new TGraphErrors();
    g_average_yield_bkg_sensor0->SetName("g_average_yield_bkg_sensor0");
    g_average_yield_bkg_sensor0->SetTitle(";Wavelength (nm);Bkg Yield");
    TGraphErrors *g_average_yield_bkg_sensor1 = new TGraphErrors();
    g_average_yield_bkg_sensor1->SetName("g_average_yield_bkg_sensor1");
    g_average_yield_bkg_sensor1->SetTitle(";Wavelength (nm);Bkg Yield");

    //  --- Average angle
    TGraphErrors *g_average_angle = new TGraphErrors();
    g_average_angle->SetName("g_average_angle");
    g_average_angle->SetTitle(";Wavelength (nm);Angle (rad)");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(3, 3);
    TCanvas *c_draw_out_of_time_bkg = new TCanvas("c_draw_out_of_time_bkg", "c_draw_out_of_time_bkg", 800, 800);
    gPad->DrawFrame(50, 0.0001, 100, 0.5);
    auto i_color = -1;
    for (auto [range, run] : run_list_filter_scan)
    {
        i_color++;
        // test_bkg_investigation(run);

        //  Get histograms
        TFile *current_file = new TFile(Form("/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/test_bkg_investigation.%s.root", run.c_str()));
        auto current_hist_sensors = (TH1F *)current_file->Get("h_signal_sensor_radius");
        auto current_hist_sensor_0 = (TH1F *)current_file->Get("h_signal_sensor0_radius");
        auto current_hist_sensor_1 = (TH1F *)current_file->Get("h_signal_sensor1_radius");
        auto current_hist_sensors_bkg = (TH1F *)current_file->Get("h_background_sensor_radius");
        auto current_hist_sensor_0_bkg = (TH1F *)current_file->Get("h_background_sensor0_radius");
        auto current_hist_sensor_1_bkg = (TH1F *)current_file->Get("h_background_sensor1_radius");

        //  Outoftime bkg investigation
        c_draw_out_of_time_bkg->cd();
        current_hist_sensors_bkg->Draw("SAME EP");
        current_hist_sensors_bkg->SetMarkerStyle(20);
        current_hist_sensors_bkg->SetMarkerColor(palette[i_color]);

        c1->cd(1 + i_color);
        gPad->SetLogy();
        current_hist_sensors->GetXaxis()->SetRangeUser(40, 100);
        current_hist_sensors->Draw("SAME EP");
        current_hist_sensors->SetMarkerStyle(20);
        current_hist_sensors->Add(bkg_hist_sensor, -1);
        current_hist_sensor_0->Draw("SAME HIST");
        current_hist_sensor_0->SetMarkerStyle(0);
        current_hist_sensor_0->SetLineColor(kGreen - 2);
        current_hist_sensor_0->Add(bkg_hist_sensor0, -1);
        current_hist_sensor_1->Draw("SAME HIST");
        current_hist_sensor_1->SetMarkerStyle(0);
        current_hist_sensor_1->SetLineColor(kRed - 1);
        current_hist_sensor_1->Add(bkg_hist_sensor1, -1);

        // Fit
        gaus_pol2->SetParameters(current_hist_sensors->GetMaximum(), 73, 0.5, 0.001, 0.001, 0.001);

        gaus_pol2->FixParameter(3, 0.);
        gaus_pol2->FixParameter(4, 0.);
        gaus_pol2->FixParameter(5, 0.);
        current_hist_sensors->Fit("gaus_pol2", "IM", "", 50, 90);
        gaus_pol2->ReleaseParameter(3);
        gaus_pol2->ReleaseParameter(4);
        gaus_pol2->ReleaseParameter(5);
        current_hist_sensors->Fit("gaus_pol2", "IM", "", 40, 100);
        g_average_radius->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(1));
        g_average_radius->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(1));
        g_average_yield->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(0));
        g_average_yield->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(0));
        g_average_resolution->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(2));
        g_average_resolution->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(2));
        // Calculate the bkg yield
        pol2_bkg->SetParameters(gaus_pol2->GetParameter(3), gaus_pol2->GetParameter(4), gaus_pol2->GetParameter(5));
        double bkg_yield = pol2_bkg->Integral(40, 100);
        double bkg_yield_err = pol2_bkg->IntegralError(40, 100);
        g_average_yield_bkg->SetPoint(i_color, (range[0] + range[1]) / 2., bkg_yield);
        g_average_yield_bkg->SetPointError(i_color, (range[1] - range[0]) / 2., bkg_yield_err);

        current_hist_sensor_0->Fit("gaus_pol2", "QIM", "", 40, 100);
        g_average_radius_sensor0->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(1));
        g_average_radius_sensor0->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(1));
        g_average_yield_sensor0->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(0));
        g_average_yield_sensor0->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(0));
        g_average_resolution_sensor0->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(2));
        g_average_resolution_sensor0->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(2));
        // Calculate the bkg yield
        pol2_bkg->SetParameters(gaus_pol2->GetParameter(3), gaus_pol2->GetParameter(4), gaus_pol2->GetParameter(5));
        bkg_yield = pol2_bkg->Integral(40, 100);
        bkg_yield_err = pol2_bkg->IntegralError(40, 100);
        g_average_yield_bkg_sensor0->SetPoint(i_color, (range[0] + range[1]) / 2., bkg_yield);
        g_average_yield_bkg_sensor0->SetPointError(i_color, (range[1] - range[0]) / 2., bkg_yield_err);

        current_hist_sensor_1->Fit("gaus_pol2", "QIM", "", 40, 100);
        g_average_radius_sensor1->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(1));
        g_average_radius_sensor1->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(1));
        g_average_yield_sensor1->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(0));
        g_average_yield_sensor1->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(0));
        g_average_resolution_sensor1->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(2));
        g_average_resolution_sensor1->SetPointError(i_color, (range[1] - range[0]) / 2., gaus_pol2->GetParError(2));
        // Calculate the bkg yield
        pol2_bkg->SetParameters(gaus_pol2->GetParameter(3), gaus_pol2->GetParameter(4), gaus_pol2->GetParameter(5));
        bkg_yield = pol2_bkg->Integral(40, 100);
        bkg_yield_err = pol2_bkg->IntegralError(40, 100);
        g_average_yield_bkg_sensor1->SetPoint(i_color, (range[0] + range[1]) / 2., bkg_yield);
        g_average_yield_bkg_sensor1->SetPointError(i_color, (range[1] - range[0]) / 2., bkg_yield_err);

        /*
        current_hist->Fit("gaus_pol2", "QIM", "", 40, 105);
        current_hist->Fit("gaus_pol2", "QIM", "", 40, 105);
        current_hist->Fit("gaus_pol2", "QIM", "", 40, 105);
        current_hist->Fit("gaus_pol2", "IM", "", 40, 105);
        g_average_radius->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(1));
        g_average_radius->SetPointError(i_color, 0., gaus_pol2->GetParError(1));
        g_average_yield->SetPoint(i_color, (range[0] + range[1]) / 2., gaus_pol2->GetParameter(0));
        g_average_yield->SetPointError(i_color, 0., gaus_pol2->GetParError(0));

        auto arm_length_aerogel = 373 + 75 - 75;
        auto current_theta = measure_ctheta(gaus_pol2->GetParameter(1), arm_length_aerogel);
        auto current_beta = get_beta(particle_mass["pion+-"], 11.5);
        auto ref_index = get_refindex(current_theta, current_beta);
        g_average_angle->SetPoint(i_color, (range[0] + range[1]) / 2., ref_index);

        */
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    auto hframe = gPad->DrawFrame(100, 73.5, 700, 77.5);
    hframe->SetTitle(";Wavelength (nm);R (mm)");
    g_average_radius->Draw("SAME LPE");
    g_average_radius->SetLineColor(kBlack);
    g_average_radius->SetMarkerColor(kBlack);
    g_average_radius->SetMarkerStyle(20);
    g_average_radius_sensor0->Draw("SAME LPE");
    g_average_radius_sensor0->SetLineColor(kGreen - 2);
    g_average_radius_sensor1->Draw("SAME LPE");
    g_average_radius_sensor1->SetLineColor(kRed - 1);

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
    hframe = gPad->DrawFrame(100, 0.01, 700, 0.4);
    hframe->SetTitle(";Wavelength (nm);Yield (a.u.)");
    g_average_yield->Draw("SAME LPE");
    g_average_yield->SetLineColor(kBlack);
    g_average_yield->SetMarkerColor(kBlack);
    g_average_yield->SetMarkerStyle(20);
    g_average_yield_sensor0->Draw("SAME LPE");
    g_average_yield_sensor0->SetLineColor(kGreen - 2);
    g_average_yield_sensor1->Draw("SAME LPE");
    g_average_yield_sensor1->SetLineColor(kRed - 1);
    TGraphErrors *g_average_yield_ratio = new TGraphErrors();
    g_average_yield_ratio->SetName("g_average_yield_ratio");
    g_average_yield_ratio->SetTitle(";Wavelength (nm);Yield Ratio");
    for (auto i_point = 0; i_point < g_average_yield->GetN(); i_point++)
    {
        auto x = g_average_yield_sensor1->GetPointX(i_point);
        auto ex = g_average_yield_sensor1->GetErrorX(i_point);
        auto y1 = g_average_yield_sensor1->GetPointY(i_point);
        auto y0 = g_average_yield_sensor0->GetPointY(i_point);
        auto ey1 = g_average_yield_sensor1->GetErrorY(i_point);
        auto ey0 = g_average_yield_sensor0->GetErrorY(i_point);
        auto y = y1 / y0;
        auto ey = y * sqrt((ey0 / y0) * (ey0 / y0) + (ey1 / y1) * (ey1 / y1));
        g_average_yield_ratio->SetPoint(i_point, x, y);
        g_average_yield_ratio->SetPointError(i_point, ex, ey);
    }

    TCanvas *c3_1 = new TCanvas("c3_1", "c3", 800, 800);
    g_average_yield_ratio->Draw("A LPE");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
    hframe = gPad->DrawFrame(100, 0.01, 700, 3.);
    hframe->SetTitle(";Wavelength (nm);Resolution (mm)");
    g_average_resolution->Draw("SAME LPE");
    g_average_resolution->SetLineColor(kBlack);
    g_average_resolution->SetMarkerColor(kBlack);
    g_average_resolution->SetMarkerStyle(20);
    g_average_resolution_sensor0->Draw("SAME LPE");
    g_average_resolution_sensor0->SetLineColor(kGreen - 2);
    g_average_resolution_sensor1->Draw("SAME LPE");
    g_average_resolution_sensor1->SetLineColor(kRed - 1);

    TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
    hframe = gPad->DrawFrame(100, 0.01, 700, 1.);
    hframe->SetTitle(";Wavelength (nm);Bkg Yield (a.u.)");
    g_average_yield_bkg_sensor0->Draw("SAME LPE");
    g_average_yield_bkg_sensor0->SetLineColor(kBlack);
    g_average_yield_bkg_sensor0->SetMarkerColor(kBlack);
    g_average_yield_bkg_sensor0->SetMarkerStyle(20);
    g_average_yield_bkg_sensor1->Draw("SAME LPE");
    g_average_yield_bkg_sensor1->SetLineColor(kRed - 1);

    /*
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
    g_average_angle->Draw("AP");
    g_average_angle->SetMarkerStyle(20);
    */

    // bkg_hist->SetMarkerStyle(20);
    // bkg_hist->SetMarkerColor(kBlack);
    // bkg_hist->SetLineColor(kBlack);
    // bkg_hist->Draw("SAME");
}