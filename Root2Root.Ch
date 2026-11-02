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
#include "TF3.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTreeReader.h"

#include <gsl/gsl_statistics.h>

using namespace std;

int Root2Root()
{
    TFile* fileRaw = new TFile("../../../../../../../../mnt/hgfs/shared-2/2024_DATA/MCP_DATA/RAW/MCP_010_4T.root", "READ");
    TFile *file = new TFile("MCP_010_4T.root", "RECREATE");

    TH2D* h2 = (TH2D*)fileRaw->Get("anode_image_func");

    TTree* tree = new TTree("treeMCP", "treeMCP");
    double X0_tree, Y0_tree, nothing;
    tree->Branch("X0", &X0_tree, "X0/D");
    tree->Branch("Y0", &Y0_tree, "Y0/D");
    tree->Branch("MCP_HG_q1", &nothing, "MCP_HG_q1/D");
    tree->Branch("MCP_BG_q1", &nothing, "MCP_BG_q1/D");
    tree->Branch("MCP_HD_q1", &nothing, "MCP_HD_q1/D");
    tree->Branch("MCP_BD_q1", &nothing, "MCP_BD_q1/D");

    nothing = 0;
    for (int i = 0; i < 1000000; i++)
    {
        h2->GetRandom2(X0_tree, Y0_tree);
        tree->Fill();
    }

    h2->SetTitle("h_Image");
    h2->SetName("h_Image"); 
    h2->Write();

    tree->Write();

    file->Close();
    fileRaw->Close();

    return 0;
}