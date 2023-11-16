#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "TMath.h"
#include "Riostream.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTreeReader.h"

#include <gsl/gsl_statistics.h>

using namespace std;

const int NparX = 5;
const int NparY = 5;
const int Nparam = (NparX + 1) * (NparY + 1);
float X2_tree, Y2_tree;

int main(int argc, char **argv)
{
  char *fRootName = argv[1];

  // Open ROOT file
  TFile *root_file = new TFile(fRootName, "UPDATE");
  TTreeReader Reader("treeMCP", root_file);
  TTreeReaderValue<double> X0_tree(Reader, "X0");
  TTreeReaderValue<double> Y0_tree(Reader, "Y0");
  TTree *tree_MCP_corr = new TTree("treeMCP_corr", "treeMCP_corr");
  tree_MCP_corr->Branch("X2", &X2_tree, "X2/F");
  tree_MCP_corr->Branch("Y2", &Y2_tree, "Y2/F");
  TH2F *h_Image_corr = new TH2F("h_Image_corr", "h_Image_corr", 500, -10., 10., 500, -10., 10.);

  // fit function
    Double_t fit_points(Double_t *, Double_t *);

  if (argc > 2)
  {
    // Open data file
    int N = 0;
    string line;
    ifstream testfile("coordinate_file_back.txt");
    while (std::getline(testfile, line)){N++;}
    ifstream infile("coordinate_file_back.txt");
    

    // Inititalise varaibales
    double X0, Y0, X1, Y1, X2, Y2;
    double tab_X0[N];
    double tab_Y0[N];
    double tab_X1[N];
    double tab_Y1[N];

    double diff_r;
    double r;

    TF2 *f_X1 = new TF2("fit_points_fX1", fit_points, -1., 1., -1., 1., Nparam);
    TF2 *f_Y1 = new TF2("fit_points_fY1", fit_points, -1., 1., -1., 1., Nparam);

    for (int i = 0; i < Nparam; i++)
    {
      f_X1->SetParameter(i, 0.);
      f_Y1->SetParameter(i, 0.);
    }

    TGraph2D *g_X1 = new TGraph2D();
    TGraph2D *g_Y1 = new TGraph2D();

    TCanvas *c1 = new TCanvas("Residus", "Residus", 800, 800);
    TGraph *Residus = new TGraph();

    for (int i = 0; i < N; i++)
    {

      infile >> X1 >> Y1 >> X0 >> Y0;

      tab_X0[i] = X0;
      tab_Y0[i] = Y0;
      tab_X1[i] = X1;
      tab_Y1[i] = Y1;

      g_X1->SetPoint(i, X0, Y0, X1);
      g_Y1->SetPoint(i, X0, Y0, Y1);
    }

    g_X1->Fit("fit_points_fX1");
    g_Y1->Fit("fit_points_fY1");

    for (int i = 0; i < N; i++)
    {
      X2 = f_X1->Eval(tab_X0[i], tab_Y0[i]);
      Y2 = f_Y1->Eval(tab_X0[i], tab_Y0[i]);

      diff_r = sqrt(pow(X2 - tab_X1[i], 2) + pow(Y2 - tab_Y1[i], 2));
      r = sqrt(pow(tab_X1[i], 2) + pow(tab_Y1[i], 2));
      Residus->SetPoint(i, r, diff_r);
    }

    size_t size = sizeof(tab_X0) / sizeof(tab_X0[0]);

    while (Reader.Next())
    {
      if (*X0_tree < gsl_stats_max(tab_X0, 1, size) && *X0_tree > gsl_stats_min(tab_X0, 1, size) && *Y0_tree < gsl_stats_max(tab_Y0, 1, size) && *Y0_tree > gsl_stats_min(tab_Y0, 1, size))
      {
        X2_tree = f_X1->Eval(*X0_tree, *Y0_tree);
        Y2_tree = f_Y1->Eval(*X0_tree, *Y0_tree);
        tree_MCP_corr->Fill();
        h_Image_corr->Fill(X2_tree, Y2_tree);
      }
    }

    ofstream out_resultats("fit_params.txt");

    out_resultats << gsl_stats_max(tab_X0, 1, size) << "\t" << gsl_stats_min(tab_X0, 1, size) << "\t" << gsl_stats_max(tab_Y0, 1, size) << "\t" << gsl_stats_min(tab_Y0, 1, size) << endl;
    for (int i = 0; i < Nparam; i++)
      out_resultats << f_X1->GetParameter(i) << "	";
    out_resultats << endl;
    for (int i = 0; i < Nparam; i++)
      out_resultats << f_Y1->GetParameter(i) << "	";

    Residus->Draw("AP");
    c1->Write();
    root_file->Write();
    root_file->Close();

    return (0);
  }
  else
  {
    // Open data file
    ifstream infile("fit_params.txt");
    double xmax, xmin, ymax, ymin;
    infile >> xmax >> xmin >> ymax >> ymin;

    std::string line;
    std::vector<std::vector<double>> data;
    while (std::getline(infile, line))
    {
      std::vector<double> values;
      std::istringstream iss(line);
      double value;
      while (iss >> value)
      {
        values.push_back(value);
      }
      data.push_back(values);
    }

    infile.close();

    double *param_X = new double[data[1].size()];
    std::copy(data[1].begin(), data[1].end(), param_X);

    double *param_Y = new double[data[2].size()];
    std::copy(data[2].begin(), data[2].end(), param_Y);

    while (Reader.Next())
    {
      if (*X0_tree < xmax && *X0_tree > xmin && *Y0_tree < ymax && *Y0_tree > ymin)
      {
        double_t data[2] = {*X0_tree, *Y0_tree};
        X2_tree = fit_points(data, param_X);
        Y2_tree = fit_points(data, param_Y);
        tree_MCP_corr->Fill();
        h_Image_corr->Fill(X2_tree, Y2_tree);
      }
    }
    root_file->Write();
    root_file->Close();
    return (0);
  }
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
Double_t fit_points(Double_t *x, Double_t *par)
{
  int i, j, ij;
  double res = 0;
  for (int i = 0; i <= NparX; i++)
    for (int j = 0; j <= NparY; j++)
    {
      ij = i * NparY + j;
      res = res + par[ij] * pow(x[0], i) * pow(x[1], j);
    }

  return res;
}
