#include "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/reconstruction/lib/testbeam.h"

//  General physics utilities
const float alhpa_EM = 1. / 137.055;   // [n.u.]
const float kElectronMass = 0.0005110; // [GeV]
const float kMuonMass = 0.1057;        // [GeV]
const float kPionMass = 0.1396;        // [GeV]
const float kKaonMass = 0.4937;        // [GeV]
const float kProtonMass = 0.9383;      // [GeV]
const float kcLight = 299792458.;      // [m/s]

//  General functions utilities
// https://root-forum.cern.ch/t/fitting-a-poisson-distribution-to-a-histogram/12078/4
//  [0] = Normalizing parameter
//  [1] / [2] -> mean (mu)
//  x / [2] -> x
//  Gamma( x / [2] + 1 ) = factorial (x / [2])
TF1 *fPhotonFitFunction = new TF1("hPhotonFitFunction", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", -1000, 1000);
TF1 *fDoublePhotonFitFunction = new TF1("fDoublePhotonFitFunction", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)+[3]*TMath::Power(([4]/[5]),(x/[5]))*(TMath::Exp(-([4]/[5])))/TMath::Gamma((x/[5])+1)", -1000, 1000);

//  === Kinematics
inline float get_relativistic_beta(float gamma) { return TMath::Sqrt(1 - gamma * gamma); }
inline float get_relativistic_beta(float mass, float momentum) { return 1. / sqrt(1 + (mass * mass) / (momentum * momentum)); }
inline float get_relativistic_gamma(float beta) { return 1. / (TMath::Sqrt(1 - beta * beta)); }
inline float get_relativistic_gamma(float mass, float momentum) { return get_relativistic_gamma(get_relativistic_beta(mass, momentum)); }

//  === Cherenkov
//  === === A priori knowledge
inline float get_theta_cherenkov(float beta, float ref_index) { return TMath::ACos(1. / (beta * ref_index)); }
inline float get_theta_atan(float radius, float arm_length) { return TMath::ATan(radius / arm_length); }
inline float get_expected_radius(float mass, float momentum, float ref_index, float arm_length) { return TMath::Tan(get_theta_cherenkov(get_relativistic_beta(mass, momentum), ref_index)) * arm_length; }
inline float get_relativistic_beta_cherenkov_threshold(float ref_index) { return 1. / ref_index; }
inline float get_momentum_cherenkov_threshold(float ref_index, float mass) { return mass * get_relativistic_beta_cherenkov_threshold(ref_index) * get_relativistic_gamma(get_relativistic_beta_cherenkov_threshold(ref_index)); }
inline float get_mass_cherenkov_threshold(float ref_index, float momentum) { return momentum / (get_relativistic_beta_cherenkov_threshold(ref_index) * get_relativistic_gamma(get_relativistic_beta_cherenkov_threshold(ref_index))); }
inline float get_expected_photons(float ref_index, float mass, float momentum, int z_charge = 1, float lambda_low = 350.e-9, float lambda_high = 500.e-9)
{
    return ((2 * TMath::Pi() * z_charge * z_charge * alhpa_EM) * ((1. / lambda_low) - (1. / lambda_high)) * (1. - 1. / (get_relativistic_beta(mass, momentum) * get_relativistic_beta(mass, momentum) * ref_index * ref_index)));
}

//  === === Mirrors
inline float get_aerogel_initpoint(float mirror_position, int n_tiles = 2) { return 20 * (n_tiles) + 52 + 262.5 + mirror_position; }
inline float get_gas_initpoint(float mirror_position) { return (-10 + 1273 - mirror_position); }
//  === === Measurements
inline float measure_relativistic_beta(float ref_index, float radius, float arm_length) { return 1. / (ref_index * TMath::Cos(get_theta_atan(radius, arm_length))); }
inline float measure_arm_length(float ref_index, float mass, float momentum, float radius) { return radius / TMath::Tan(get_theta_cherenkov(get_relativistic_beta(mass, momentum), ref_index)); }
inline float measure_mass_hypothesis(float momentum, float ref_index, float radius, float arm_length) { return TMath::Sqrt((momentum * momentum) / (measure_relativistic_beta(ref_index, radius, arm_length) * measure_relativistic_beta(ref_index, radius, arm_length)) - (momentum * momentum)); }
//  === === Graphics
// inline void plot_prediction(float ref_index, float momentum, float mass, float distance, std::array<float, 2> xy = {0., 0.}) { return plot_circle({xy[0], xy[1], get_expected_radius(mass, momentum, ref_index, distance)}, kRed, kDashed, 2); }

std::array<float, 2>
eval_with_errors(TGraphErrors *gTarget, float _xtarget)
{
    auto iPnt = 0;
    auto current_x = 1.;
    auto previous_x = 1.;
    float interpolated_y = -1.;
    float interpolated_ey = -1.;
    for (iPnt = 1; iPnt < gTarget->GetN() - 1; iPnt++)
    {
        current_x = gTarget->GetPointX(iPnt);
        previous_x = gTarget->GetPointX(iPnt - 1);
        if ((current_x >= _xtarget) && (previous_x <= _xtarget))
        {
            break;
        }
    }
    if (iPnt == gTarget->GetN() - 1)
    {
        return {interpolated_y, interpolated_ey};
    }
    auto coeff_y1 = (_xtarget - previous_x) / (current_x - previous_x);
    auto coeff_y0 = (current_x - _xtarget) / (current_x - previous_x);
    interpolated_y = coeff_y1 * gTarget->GetPointY(iPnt) + coeff_y0 * gTarget->GetPointY(iPnt - 1);
    interpolated_ey = sqrt((coeff_y0 * gTarget->GetErrorY(iPnt)) * (coeff_y0 * gTarget->GetErrorY(iPnt)) + (coeff_y1 * gTarget->GetErrorY(iPnt - 1)) * (coeff_y1 * gTarget->GetErrorY(iPnt - 1)));
    return {interpolated_y, interpolated_ey};
}
std::array<float, 2>
eval_with_errors(TGraph *gTarget, float _xtarget)
{
    auto iPnt = 0;
    auto current_x = 1.;
    auto previous_x = 1.;
    float interpolated_y = -1.;
    float interpolated_ey = -1.;
    for (iPnt = 1; iPnt < gTarget->GetN() - 1; iPnt++)
    {
        current_x = gTarget->GetPointX(iPnt);
        previous_x = gTarget->GetPointX(iPnt - 1);
        if ((current_x >= _xtarget) && (previous_x <= _xtarget))
        {
            break;
        }
    }
    if (iPnt == gTarget->GetN() - 1)
    {
        return {interpolated_y, interpolated_ey};
    }
    auto coeff_y1 = (_xtarget - previous_x) / (current_x - previous_x);
    auto coeff_y0 = (current_x - _xtarget) / (current_x - previous_x);
    interpolated_y = coeff_y1 * gTarget->GetPointY(iPnt) + coeff_y0 * gTarget->GetPointY(iPnt - 1);
    interpolated_ey = 0.;
    return {interpolated_y, interpolated_ey};
}

void test()
{
    /*
    TFile *fout = new TFile("PDE.root", "RECREATE");
    auto test = new TGraph("13_50_PDE.csv", "%lg, %lg");
    test->SetName("13_50_PDE");
    for (auto i_pnt = 0; i_pnt < test->GetN(); i_pnt++)
    {
        auto _x = test->GetPointX(i_pnt);
        auto _y = test->GetPointY(i_pnt);
        test->SetPoint(i_pnt, _x, 60 - _y);
    }
    test->Sort();
    test->Write();
    test = new TGraph("13_75_PDE.csv", "%lg, %lg");
    test->SetName("13_75_PDE");
    for (auto i_pnt = 0; i_pnt < test->GetN(); i_pnt++)
    {
        auto _x = test->GetPointX(i_pnt);
        auto _y = test->GetPointY(i_pnt);
        test->SetPoint(i_pnt, _x, 60 - _y);
    }
    test->Sort();
    test->Write();
    test = new TGraph("14_50_PDE.csv", "%lg, %lg");
    test->SetName("14_50_PDE");
    for (auto i_pnt = 0; i_pnt < test->GetN(); i_pnt++)
    {
        auto _x = test->GetPointX(i_pnt);
        auto _y = test->GetPointY(i_pnt);
    }
    test->Sort();
    test->Write();
*/

    //  Loop over PDE of sensors
    TFile *fin = new TFile("PDE.root");
    TGraph *pde_13_50 = (TGraph *)(fin->Get("13_50_PDE"));
    TGraph *pde_13_75 = (TGraph *)(fin->Get("13_75_PDE"));
    TGraph *pde_14_50 = (TGraph *)(fin->Get("14_50_PDE"));

    //  Variables for Frank-Tamm
    double *x = new double[1];
    double *p = new double[2];
    x[0] = get_relativistic_beta(kElectronMass, 11.5);

    //  Graphs to check photon yield
    TGraphErrors *g13_50 = new TGraphErrors();
    TGraphErrors *g13_75 = new TGraphErrors();

    //  Loop from 200 to 900 nm
    auto nph_13_50 = 0.;
    auto nph_13_75 = 0.;
    for (auto i_lambda = 280; i_lambda <= 850; i_lambda++)
    {
        //  Calculate cherenkov photons emission
        p[0] = i_lambda * 1.e-9;
        p[1] = 1.0210;
        auto ch_photons_1 = 1.e-11 * testbeam::simp_frank_tamm(x, p);
        p[1] = 1.0207;
        auto ch_photons_2 = 1.e-11 * testbeam::simp_frank_tamm(x, p);

        //
        g13_50->SetPoint(i_lambda - 280, i_lambda, (ch_photons_1 + ch_photons_2) * 0.01 * eval_with_errors(pde_13_50, i_lambda)[0]);
        g13_75->SetPoint(i_lambda - 280, i_lambda, (ch_photons_1 + ch_photons_2) * 0.01 * eval_with_errors(pde_13_75, i_lambda)[0]);

        //  Integral
        nph_13_50 += (2 * ch_photons_1 + 2 * ch_photons_2) * 0.01 * eval_with_errors(pde_13_50, i_lambda)[0];
        nph_13_75 += (2 * ch_photons_1 + 2 * ch_photons_2) * 0.01 * eval_with_errors(pde_13_75, i_lambda)[0];
    }
    cout << nph_13_50 << endl;
    cout << nph_13_75 << endl;

    TCanvas *c1 = new TCanvas();
    auto h1 = gPad->DrawFrame(250, 0, 850, 0.16, ";photon wavelength (nm);N_{#gamma} (N / cm nm)");

    g13_75->SetLineColor(kRed);
    g13_75->SetMarkerColor(kRed);
    g13_75->SetMarkerStyle(kFullCircle);
    g13_75->Draw("SAME LPE");
    g13_75->SetName("13_75_PDE");

    g13_50->SetLineColor(kBlue);
    g13_50->SetMarkerColor(kBlue);
    g13_50->SetMarkerStyle(kFullCircle);
    g13_50->Draw("SAME LPE");
    g13_50->SetName("13_50_PDE");

    gPad->BuildLegend();
}

/*
void test()
{
    auto f_dcr = new TFile("dcr_test_dcr.root");
    auto graph_dcr = (TGraphErrors *)(f_dcr->Get("g_frequencies"));
    auto f_nodcr = new TFile("dcr_test_nodcr.root");
    auto graph_nodcr = (TGraphErrors *)(f_nodcr->Get("g_frequencies"));

    TH2F *h_dcr = new TH2F("h_dcr", ";x (mm);y (mm)", 1980, -99, 99, 1980, -99, 99);
    TH2F *h_nodcr = new TH2F("h_nodcr", ";x (mm);y (mm)", 1980, -99, 99, 1980, -99, 99);
    TH2F *h_ratdcr = new TH2F("h_ratdcr", ";x (mm);y (mm)", 1980, -99, 99, 1980, -99, 99);
    TH1F *h_rat_val = new TH1F("h_rat_val", "", 4000, 0, 4000);

    TGraphErrors *g_freq_corr = new TGraphErrors();
    g_freq_corr->GetXaxis()->SetTitle("DCR start of spill (Hz)");
    g_freq_corr->GetYaxis()->SetTitle("DCR out-of-time (Hz)");

    TRandom *g_random = new TRandom();

    for (auto ipnt = 0; ipnt < graph_dcr->GetN(); ipnt++)
    {
        auto pos_x = graph_dcr->GetPointX(ipnt);
        auto pos_y = graph_dcr->GetPointY(ipnt);
        auto dcr_v = graph_dcr->GetErrorX(ipnt);
        auto dcr_e = graph_dcr->GetErrorY(ipnt);
        auto nodcr_v = graph_nodcr->GetErrorX(ipnt);
        auto nodcr_e = graph_nodcr->GetErrorY(ipnt);

        g_freq_corr->SetPoint(ipnt, dcr_v, nodcr_v);
        // g_freq_corr->SetPointError(ipnt, dcr_e, nodcr_e);


                for (auto ipnt = 0; ipnt < dcr_v * 10; ipnt++)
                    h_dcr->Fill(pos_x + 1.5 * g_random->Uniform(-1, 1), pos_y + 1.5 * g_random->Uniform(-1, 1));
                for (auto ipnt = 0; ipnt < nodcr_v; ipnt++)
                    h_nodcr->Fill(pos_x + 1.5 * g_random->Uniform(-1, 1), pos_y + 1.5 * g_random->Uniform(-1, 1));
                for (auto ipnt = 0; ipnt < nodcr_v / dcr_v; ipnt++)
                    h_ratdcr->Fill(pos_x + 1.5 * g_random->Uniform(-1, 1), pos_y + 1.5 * g_random->Uniform(-1, 1));


        h_rat_val->Fill(nodcr_v / dcr_v);
    }

    TF1 *fPhotonFitFunction = new TF1("hPhotonFitFunction", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", -10000, 10000);
    fPhotonFitFunction->SetParameters(30, 60, 20);

    new TCanvas();
    h_dcr->Draw("COLZ");
    new TCanvas();
    h_nodcr->Draw("COLZ");
    new TCanvas();
    h_ratdcr->Draw("COLZ");
    new TCanvas();
    h_rat_val->Draw("");
    // h_rat_val->Fit(fPhotonFitFunction);

    new TCanvas();
    TH1F *hframe = new TH1F("", ";DCR start of spill (Hz);DCR out-of-time (Hz)", 1.e5, 1.e2, 1.e7);
    hframe->SetMinimum(1.e2);
    hframe->SetMaximum(1.e7);
    hframe->Draw();
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    gPad->SetLogy();
    g_freq_corr->SetMarkerStyle(20);
    g_freq_corr->Draw("SAME PE ");
}

void process_runs(std::string path = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/Data/20231010-163636/");

void test()
{
    std::string start_path = "/Users/nrubini/Analysis/ePIC/_production_repositories/ePIC-dRICH-beamtest-analysis/Data/";
    std::vector<std::string> list_of_runs = {
        //"20231011-174928", // 49.0V
        //"20231011-214110", // 49.5V
        //"20231011-181206", // 50.0V
        "20231011-223457", // 50.5V
        "20231011-182852", // 51.0V
        "20231011-221500", // 51.5V
        "20231011-191435", // 52.0V
        "20231011-225350", // 52.5V
        "20231010-144804", // 53.0V
        "20231010-152047", // 54.0V
        "20231010-154910"  // 55.0V
    };
    std::vector<float> list_of_vbias = {
        // 49.,
        // 50.,
        50.5,
        51.0,
        51.5,
        52.0,
        52.5,
        53.0,
        54.0,
        55.0};

    std::string label1, label2;
    double value1, value2;
    std::vector<TGraphErrors *> final_graph;
    for (auto i = 0; i < 4; i++)
        final_graph.push_back(new TGraphErrors());

    auto jter = -1;
    for (auto current_run : list_of_runs)
    {
        jter++;
        auto iter = 0;
        recowriter(start_path + current_run + "/lightdata.root", start_path + current_run + "/recodata.root");
        process_runs(start_path + current_run);
        std::ifstream file(start_path + current_run + "/out.txt");
        if (!file)
            continue;
        while (file >> label1 >> label2 >> value1 >> value2)
        {
            std::cout << label1 << " " << label2 << " " << value1 << " " << value2 << std::endl;
            final_graph[iter]->SetNameTitle((label1 + " " + label2).c_str(), (label1 + " " + label2).c_str());
            final_graph[iter]->SetPoint(jter, list_of_vbias[jter], value1);
            final_graph[iter]->SetPointError(jter, 0, value2);
            iter++;
        }
    }

    for (auto current_graph : final_graph)
    {
        new TCanvas();
        current_graph->Draw("ALPE");
    }
}

void process_runs(std::string path)
{
    //  Load data from run
    testbeam::data current_data(path + "/recodata.root");
    auto current_tree = current_data.get_tree();

    //  Event selection
    TH1F *hEventCounting = new TH1F("hEventCounting", "", 10, 0, 10);
    TH2F *hMapTest0 = new TH2F("hMapTest0", "Evs. w/ less than 4 hits", 100, -100, 100, 100, -100, 100);
    TH2F *hMapTest1 = new TH2F("hMapTest1", "Evs. w/ no clusters", 100, -100, 100, 100, -100, 100);
    TH2F *hMapTest2 = new TH2F("hMapTest2", "Evs. w/ bigger cluster w/ less than 3 hits", 100, -100, 100, 100, -100, 100);
    TH2F *hMapTest3 = new TH2F("hMapTest3", "Evs. w/ radius outsied range", 100, -100, 100, 100, -100, 100);
    TH2F *hMapTest4 = new TH2F("hMapTest4", "Accepted events", 100, -100, 100, 100, -100, 100);

    //
    TH1F *h_fit_radius = new TH1F("h_fit_radius", "", 1000, 0, 100);
    TH2F *h3_fit_radius = new TH2F("h3_fit_radius", "", 100, -10, 10, 40, 0, 40);
    TH2F *h2_fit_radius = new TH2F("h2_fit_radius", "", 1000, 0, 100, 40, 0, 40);

    //  Photon counting
    TH1F *h_photon_counting_DCR = new TH1F("h_photon_counting_DCR", ";N_{#gamma}", 100, 0, 100);
    TH1F *h_photon_counting_bkg = new TH1F("h_photon_counting_bkg", ";N_{#gamma}", 100, 0, 100);
    TH1F *h_photon_counting_all = new TH1F("h_photon_counting_all", ";N_{#gamma}", 100, 0, 100);
    TH2F *h_photon_counting_bkg_all = new TH2F("h_photon_counting_bkg_all", ";N_{#gamma} (bkg);N_{#gamma} (sig)", 100, 0, 100, 100, 0, 100);

    for (auto iev = 0; iev < current_tree->GetEntries(); iev++)
    {
        current_tree->GetEntry(iev);
        if (current_data.get_n() < 4)
        {
            for (auto i = 0; i < 1; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hMapTest0->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        auto cluster_points = current_data.get_r_clusters(5);
        if (!cluster_points.size())
        {
            for (auto i = 0; i < 2; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hMapTest1->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }
        if (cluster_points[0].size() < 3)
        {
            for (auto i = 0; i < 3; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hMapTest2->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }

        auto good_points_array = current_data.select_points({current_data.get_common_center_x(), current_data.get_common_center_y(), current_data.get_common_radius(), current_data.get_timing_center()}, {999, 999, 3, 4 * current_data.get_timing_sigma()});
        auto good_points = good_points_array[0];
        if (good_points.size() < 3)
        {
            for (auto i = 0; i < 3; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hMapTest2->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }
        //  Fit selected points
        auto current_fit_results = current_data.fit_circle(good_points, {0., 0., 0.}, true);
        if (fabs(current_fit_results[2][0] - current_data.get_common_radius()) > 3)
        {
            for (auto i = 0; i < 4; i++)
                hEventCounting->Fill(i);
            for (auto i = 0; i < current_data.get_n(); i++)
                hMapTest3->Fill(current_data.get_x(i), current_data.get_y(i));
            continue;
        }
        for (auto i = 0; i < 5; i++)
            hEventCounting->Fill(i);
        for (auto i = 0; i < current_data.get_n(); i++)
            hMapTest4->Fill(current_data.get_x(i), current_data.get_y(i));

        //  Fill radius
        h_fit_radius->Fill(current_fit_results[2][0]);
        h2_fit_radius->Fill(current_fit_results[2][0], good_points.size());
        auto diff_circle_vals = current_data.diff_circle(good_points);
        for (auto current_diff : diff_circle_vals)
            h3_fit_radius->Fill(current_diff, good_points.size());
        h_photon_counting_bkg->Fill(current_data.get_n() - good_points.size());
        h_photon_counting_all->Fill(good_points.size());
        auto count_dcr = 0;
        for (auto ihit = 0; ihit < current_data.get_n(); ihit++)
            if (current_data.is_dcr(ihit))
                count_dcr++;
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

    new TCanvas();
    h_photon_counting_bkg_all->Draw("COLZ");

    new TCanvas();
    TGraphErrors *test = new TGraphErrors();
    for (auto ybin = 1; ybin <= h_photon_counting_bkg_all->GetNbinsY(); ybin++)
    {
        //  Get current slice under test
        auto current_point = test->GetN();
        auto current_slice = h_photon_counting_bkg_all->ProjectionX(Form("proj_%i", ybin), ybin, ybin);
        if (current_slice->GetEntries() < 100)
            continue;

        auto y_low_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(ybin);
        auto y_high_edge = h_photon_counting_bkg_all->GetYaxis()->GetBinLowEdge(ybin + 1);
        current_slice->Fit(fPhotonFitFunction, "I");

        test->SetPoint(current_point, 0.5 * (y_high_edge + y_low_edge) - 0.5, fPhotonFitFunction->GetParameter(1));
        test->SetPointError(current_point, 0.5 * (y_high_edge - y_low_edge), fPhotonFitFunction->GetParError(1));
    }

    test->Draw("ALPE");

    TFile *outfile = new TFile((path + "/preliminary_test.root").c_str(), "RECREATE");
    hEventCounting->Write();
    hMapTest0->Write();
    hMapTest1->Write();
    hMapTest2->Write();
    hMapTest3->Write();
    hMapTest4->Write();

    h_photon_counting_DCR->Write();
    h_photon_counting_bkg->Write();
    h_photon_counting_all->Write();

    h_photon_counting_bkg->Write();
    h_photon_counting_all->Write();
    h_photon_counting_DCR->Write();
    h_photon_counting_bkg_all->Write();

    test->Write();
    outfile->Close();

    return;

    h_fit_radius->Draw();
    h_fit_radius->Fit("gaus");
    new TCanvas();
    h2_fit_radius->Draw();
    h2_fit_radius->Fit("gaus");
    new TCanvas();
    h3_fit_radius->Draw();
    h3_fit_radius->Fit("gaus");
    new TCanvas();
    hEventCounting->Draw();
    new TCanvas();
    hMapTest0->Draw("COLZ");
    new TCanvas();
    hMapTest1->Draw("COLZ");
    new TCanvas();
    hMapTest2->Draw("COLZ");
    new TCanvas();
    hMapTest3->Draw("COLZ");
    new TCanvas();
    hMapTest4->Draw("COLZ");

}
*/