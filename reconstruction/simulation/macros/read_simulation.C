#include "../lib/simulation.h"

void read_simulation(std::string input_filename = "./simulation/data/run_2041.root", int max_ev = INT_MAX, std::string output_filename = "out3.root")
{
  //  Output
  std::map<std::string, TProfile *> _TProfile;
  std::map<std::string, TH1F *> _TH1F;
  std::map<std::string, TH2F *> _TH2F;
  std::map<std::string, TH3F *> _TH3F;
  std::map<std::string, TH1D *> _TH1D;
  std::map<std::string, TH2D *> _TH2D;
  std::map<std::string, TH3D *> _TH3D;
  std::map<std::string, TGraphErrors *> _TGraphErrors;
  _TH1F["hWavelength"] = new TH1F("hAerogelWavelength", ";photon wavelength (nm);Entries", 100, 0, 1000);
  _TH1F["hWavelengthPDEcorrected"] = new TH1F("hWavelengthPDEcorrected", ";photon wavelength (nm);Entries", 100, 0, 1000);
  _TH1F["hPhotonPerEvent"] = new TH1F("hPhotonPerEvent", ";N_{#gamma}", 30, 0, 30);
  _TH2F["hXYMap"] = new TH2F("hXYMap", ";x (mm); y (mm)", 300, -300, 300, 300, -300, 300);
  _TH2F["hXYMapPDEcorrected"] = new TH2F("hXYMapPDEcorrected", ";x (mm); y (mm)", 300, -300, 300, 300, -300, 300);
  _TH2F["hEfficiencyMap"] = new TH2F("hEfficiencyMap", ";#phi (rad); R (mm)", 315, -TMath::Pi(), TMath::Pi(), 1000, 0, 1000);

  // TRandom utility
  TRandom *util_rnd = new TRandom();

  // Read tree
  simdata target_data;
  auto read_data = read_tree(input_filename, target_data);
  auto input_tree = read_data.first;
  auto input_file = read_data.second;

  //  Loop on tree
  auto max_entries = min(1. * input_tree->GetEntries(), 1. * max_ev);
  for (auto iEv = 0; iEv < max_entries; iEv++)
  {
    //  Recover event
    input_tree->GetEntry(iEv);

    //  Count interacting photons
    int n_photons = 0;

    //  Loop over entries
    for (auto iEntry = 0; iEntry < target_data.nhits; iEntry++)
    {
      //  Select photons emitted in the aerogel arriving on the detection plane
      // if (!is_aerogel_photon_on_detection_plane(target_data, iEntry))
      //  Select photons emitted in the gas arriving on the detection plane
      if (!is_aerogel_photon_on_detection_plane(target_data, iEntry))
        continue;

      // Record photon wavelength
      _TH1F["hWavelength"]->Fill(get_wavelenght(target_data.E[iEntry]));

      // Record photon wavelength
      _TH2F["hXYMap"]->Fill(target_data.hit_x[iEntry], target_data.hit_y[iEntry]);
      if (util_rnd->Uniform() < get_approximate_PDE("HPK S13360-3050VS", get_wavelenght(target_data.E[iEntry])))
      {
        _TH2F["hXYMapPDEcorrected"]->Fill(target_data.hit_x[iEntry], target_data.hit_y[iEntry]);
        _TH1F["hWavelengthPDEcorrected"]->Fill(get_wavelenght(target_data.E[iEntry]));
        n_photons++;
      }

      // First attempt at efficiency calculation
      auto current_x = target_data.hit_x[iEntry] + 10;
      auto current_y = target_data.hit_y[iEntry] - 3;
      auto current_phi = TMath::ATan2(current_y, current_x);
      auto current_radius = TMath::Sqrt(current_x * current_x + current_y * current_y);
      if (fabs(current_phi - TMath::Pi() * 0.5) < TMath::Pi() * 2 * 10. / 360.)
        continue;
      _TH2F["hEfficiencyMap"]->Fill(current_phi, current_radius);
    }

    //  Save n photons
    _TH1F["hPhotonPerEvent"]->Fill(n_photons);
  }

  TF1 *fit_center = new TF1("fit_center", "[0]+[1]*TMath::Sin(x-[2])");
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  fit_center->SetParameter(0, 45);
  fit_center->SetParameter(1, 10);
  fit_center->SetParameter(2, 1);
  auto test_prf = _TH2F["hEfficiencyMap"]->ProfileX("TPrf");
  test_prf->Draw();
  test_prf->Fit(fit_center);
  auto pred_R = fit_center->GetParameter(0);
  auto pred_x = -fit_center->GetParameter(1) * TMath::Sin(fit_center->GetParameter(2));
  auto pred_y = fit_center->GetParameter(1) * TMath::Cos(fit_center->GetParameter(2));
  cout << "Prediction: X: " << pred_x << " Y: " << pred_y << endl;

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
  _TH2F["hXYMap"]->Draw("COLZ");
  auto result = new TEllipse(pred_x, pred_y, pred_R);
  result->SetFillStyle(0);
  result->SetLineColor(kRed);
  result->DrawEllipse(pred_x, pred_y, pred_R, 0, 0, 360, 0, "same");

  //  Normalise photons
  _TH1F["hPhotonPerEvent"]->Scale(1. / _TH1F["hPhotonPerEvent"]->GetEntries());

  //  Save to file
  TFile *output_file = new TFile(output_filename.c_str(), "RECREATE");
  for (auto [label, object] : _TProfile)
    object->Write(label.c_str());
  for (auto [label, object] : _TH1F)
    object->Write(label.c_str());
  for (auto [label, object] : _TH2F)
    object->Write(label.c_str());
  for (auto [label, object] : _TH3F)
    object->Write(label.c_str());
  for (auto [label, object] : _TGraphErrors)
    object->Write(label.c_str());
  output_file->Close();
}