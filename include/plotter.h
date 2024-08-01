#pragma once

class plotter
{

public:
  plotter()
  {
    // constructor
  } // const

  ~plotter(
      //  rad_ delete;
      // time delete;
      // nhits delete;
  )
  {

    // destructor

  } // dest

  // void setBinWidhth(int bw =1){
  //    bin_width= bw ;
  // }
  // plot the circle
  void globalBeauty()
  {
    // gStyle->SetPalette(kSolar);
    //  TColor *color = gROOT->GetColor(20);
    // color->SetRGB(0.4, 0.6, 0.8); // Soft Blue
    gStyle->SetOptStat(0);
    // gPad->SetGrid();
    gStyle->SetLineWidth(3);

    gStyle->SetLegendTextSize(0.03);
    gStyle->SetLegendBorderSize(0);
  }

  void setCanvasPixel(vector<float> canshape = {1200, 1200})
  {
    canvas_x = canshape[0];
    canvas_y = canshape[1];
  }
  void plot_circle(std::array<float, 3> parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1)
  {
    auto result = new TEllipse(parameters[0], parameters[1], parameters[2]);
    result->SetFillStyle(0);
    result->SetLineColor(line_color);
    result->SetLineStyle(line_style);
    result->SetLineWidth(line_width);
    result->DrawEllipse(parameters[0], parameters[1], parameters[2], 0, 0, 360, 0, "same");
  }

  virtual void beautify(TH1D *object, int col, int style)
  {
    object->SetLineWidth(3);
    object->SetLineColor(col);
    object->SetMarkerColor(col);
    object->SetMarkerSize(2);
    // object->SetOptStat(0);
    object->SetMarkerStyle(style);
  }

  virtual void beautify(TLine *object, int col, int style)
  {
    object->SetLineWidth(3);
    object->SetLineColor(col);
    // object->SetMarkerColor(col);
    //   object->SetMarkerSize(3);
    // object->SetOptStat(0);

    object->SetLineStyle(style);
  }

  virtual void beautify(TArc *object, int col, int style)
  {
    object->SetLineWidth(3);
    object->SetLineColor(col);
    // object->SetMarkerColor(col);
    //  object->SetMarkerSize(3);
    // object->SetOptStat(0);
    object->SetLineStyle(style);
    // object->SetMarkerStyle(style);
  }

  virtual void SetBinWidthX(float bwidth)
  {
    bin_width_x = bwidth;
  }

  virtual void SetBinWidthY(float bwidth)
  {
    bin_width_y = bwidth;
  }

  // to divide the canvas

  virtual void Set_Divide_Canvas(bool dodivide, int nrow, int ncol)
  {
    divide_canvas = dodivide;
    can_row = nrow;
    can_column = ncol;
  }

  TH2D *plot_histo2d(const char *histname, vector<float> x, vector<float> y)
  {
    // TCanvas * can = new TCanvas("can", "", 1200,1200);

    if (x.empty() || y.empty())
    {
      TH2D *empty2 = new TH2D("empty2", "", 100, -10, 10, 100, -10, 10);
      cout << "size of x  is " << x.size() << endl;
      cout << "size of y  is " << y.size() << endl;

      return empty2;
    }
    // Use std::max_element to find the maximum element
    auto xmax_ = std::max_element(x.begin(), x.end());
    auto xmin_ = std::min_element(x.begin(), x.end());

    auto ymax_ = std::max_element(y.begin(), y.end());
    auto ymin_ = std::min_element(y.begin(), y.end());

    float xmax = *xmax_;
    float xmin = *xmin_;
    float ymax = *ymax_;
    float ymin = *ymin_;
    int nbin_x = (int)TMath::Abs(xmax - xmin) / bin_width_x;
    int nbin_y = (int)TMath::Abs(ymax - ymin) / bin_width_y;

    TH2D *hist_2 = new TH2D(histname, "", nbin_x, xmin, xmax, nbin_y, ymin, ymax);

    for (int i = 0; i < x.size(); i++)
    {

      hist_2->Fill(x[i], y[i]);
    } // end loop
    // can->cd();
    // hist_2->Draw(option);
    // can->SaveAs(Form("%s.png",histname));
    return hist_2;

  } //

  TH2D *plot_histo2d(const char *histname, vector<pair<float, float>> p)
  {
    // TCanvas * can = new TCanvas("can", "", 1200,1200);
    if (p.empty())
    {
      TH2D *empty3 = new TH2D("empty3", "", 100, -10, 10, 100, -10, 10);
      cout << "size of x  is " << p.size() << endl;

      return empty3;
    }

    vector<float> x;
    vector<float> y;
    for (auto pt : p)
    {

      x.push_back(pt.first);
      y.push_back(pt.second);
    }

    // Use std::max_element to find the maximum element
    auto xmax_ = std::max_element(x.begin(), x.end());
    auto xmin_ = std::min_element(x.begin(), x.end());

    auto ymax_ = std::max_element(y.begin(), y.end());
    auto ymin_ = std::min_element(y.begin(), y.end());

    float xmax = *xmax_;
    float xmin = *xmin_;
    float ymax = *ymax_;
    float ymin = *ymin_;
    int nbin_x = (int)TMath::Abs(xmax - xmin) / bin_width_x;
    int nbin_y = (int)TMath::Abs(ymax - ymin) / bin_width_y;

    TH2D *hist_2 = new TH2D(histname, "", nbin_x, xmin, xmax, nbin_y, ymin, ymax);

    for (int i = 0; i < x.size(); i++)
    {

      hist_2->Fill(x[i], y[i]);
    } // end loop
    // can->cd();
    // hist_2->Draw(option);
    // can->SaveAs(Form("%s.png",histname));
    return hist_2;

  } //

  TH1D *plot_histo1d(const char *histname, vector<float> x)
  {

    if (x.empty())
    {
      TH1D *empty4 = new TH1D("empty4", "", 100, -10, 10);
      cout << "size of x  is " << x.size() << endl;

      return empty4;
    }
    // Use std::max_element to find the maximum element
    auto xmax_ = std::max_element(x.begin(), x.end());
    auto xmin_ = std::min_element(x.begin(), x.end());

    float xmax = *xmax_;
    float xmin = *xmin_;
    int nbin_x = (int)TMath::Abs(xmax - xmin) / bin_width_x;

    // cout << "bin_width " <<  bin_width_x << "   " <<nbin_x << endl;
    TH1D *hist_1 = new TH1D(histname, "", nbin_x, xmin, xmax);

    for (int i = 0; i < x.size(); i++)
    {

      hist_1->Fill(x[i]);
    } // end loop

    return hist_1;
  }

  TH1D *plot_histo1d(const char *histname, vector<pair<float, float>> x)
  {

    if (x.empty())
    {
      TH1D *empty4 = new TH1D("empty4", "", 100, -10, 10);
      cout << "size of x  is " << x.size() << endl;

      return empty4;
    }
    // Use std::max_element to find the maximum element
    //  auto xmax_ = std::max_element(x.begin().first, x.end().first);
    // auto  xmin_ = std::min_element(x.begin().first, x.end().first);

    float xmax = 100; //*xmax_ ;
    float xmin = 0;   //*xmin_;
    int nbin_x = (int)TMath::Abs(xmax - xmin) / bin_width_x;

    // cout << "bin_width " <<  bin_width_x << "   " <<nbin_x << endl;
    TH1D *hist_1 = new TH1D(histname, "", nbin_x, xmin, xmax);

    for (int i = 0; i < x.size(); i++)
    {
      cout << "radius_Scan, " << x[i].first << "   " << x[i].second << endl;
      hist_1->Fill(x[i].first, x[i].second);
    } // end loop

    return hist_1;
  }

  TCanvas *BuildCanvas(const char *name)
  {

    TCanvas *can = new TCanvas(name, "", canvas_x, canvas_y);
    if (divide_canvas)
    {
      can->Divide(can_row, can_column);
    }

    return can;
  }

private:
  int canvas_x = 1200;
  int canvas_y = 1200;
  bool divide_canvas = false;
  int can_row = 2;
  int can_column = 2;
  float bin_width_x = 3.0;
  float bin_width_y = 3.0;
};
