
struct recodata
{
  // MC
  unsigned short N = 0;
  float X[1024];
  float Y[1024];
  float R[1024];
  unsigned short id[1024];
  // Data-like
  unsigned short n = 0;
  float x[1024];
  float y[1024];
  float t[1024];
};

void load_data(TTree *reco_data_tree, recodata &target_data_struct)
{
  reco_data_tree->SetBranchAddress("N", &target_data_struct.N);
  reco_data_tree->SetBranchAddress("X0", &target_data_struct.X);
  reco_data_tree->SetBranchAddress("Y0", &target_data_struct.Y);
  reco_data_tree->SetBranchAddress("R", &target_data_struct.R);
  reco_data_tree->SetBranchAddress("n", &target_data_struct.n);
  reco_data_tree->SetBranchAddress("x", &target_data_struct.x);
  reco_data_tree->SetBranchAddress("y", &target_data_struct.y);
  reco_data_tree->SetBranchAddress("t", &target_data_struct.t);
  reco_data_tree->SetBranchAddress("id", &target_data_struct.id);
}

void read_fastmc(std::string input_filename)
{
  TFile *input_file = new TFile(input_filename.c_str());
  TTree *input_tree = (TTree *)(input_file->Get("recodata"));
  recodata input_data;
  load_data(input_tree, input_data);

  TH2F *SigMap = new TH2F("SigMap", "SigMap", 99, -99, 99, 99, -99, 99);
  TH2F *BkgMap = new TH2F("BkgMap", "BkgMap", 99, -99, 99, 99, -99, 99);

  for (int iEv = 0; iEv < input_tree->GetEntries(); iEv++)
  {
    input_tree->GetEntry(iEv);
    for (int iHit = 0; iHit < input_data.n; iHit++)
    {
      if (input_data.id[iHit] == 0)
      {
        SigMap->Fill(input_data.x[iHit], input_data.y[iHit]);
      }
      else
      {
        BkgMap->Fill(input_data.x[iHit], input_data.y[iHit]);
      }
    }
  }

  TCanvas *c1 = new TCanvas("", "", 1500, 750);
  c1->Divide(2, 1);
  c1->cd(1);
  SigMap->Draw("COLZ");
  c1->cd(2);
  BkgMap->Draw("COLZ");
}