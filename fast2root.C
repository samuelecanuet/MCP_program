/* Convert .fast files from MCP TEST data in .root files (Tree & histos)
 * Requires fasterac & Root libraries (should be on acquisition computer)
 * XF 13/12/2021
 *
 * Usage: make, and then command ./group2tree_XF_MCP_ADC_cor
 */
 

//  std includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//  root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Riostream.h"
#include "TString.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraph.h"

//  fasterac includes
#include "fasterac/fasterac.h"
#include "fasterac/fast_data.h"
#include "fasterac/qdc.h"
#include "fasterac/spectro.h"
#include "fasterac/group.h"


using namespace std;


#define LABEL_MCP_HG   3
#define LABEL_MCP_HD   2
#define LABEL_MCP_BG   4
#define LABEL_MCP_BD   1
#define LABEL_MCP_SIG  5

#define NO_DATA 0.

const int N=60;
const int NparX=5;
const int NparY=5;
const int Nparam = (NparX+1)*(NparY+1);

const double MCP_q1_min=4000.;

#include <TFile.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TCanvas.h>

int main (int argc, char** argv) 
{

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////lecture data fast/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

  /************/
  /*  FASTER  */
  /************/
  //  file reader
  faster_file_reader_p   reader;
  //  data
  faster_data_p          data;
  unsigned char          alias;
  unsigned short         label;
  double 		             clock_ns;
  double 		             clock_sec;
  
  
  //  group data
  faster_buffer_reader_p group_reader;
  char                   group_buffer [1500];
  unsigned short         lsize;
  faster_data_p          group_data; //struct for group
  //data
  qdc_x2                 qdc2; //struct for channels
  crrc4_spectro			 spectro_data;
  
  
  /**********/
  /*  ROOT  */
  /**********/
  
    
  Double_t                  MCP_HG_q1;
  Double_t                  MCP_HD_q1;
  Double_t                  MCP_BG_q1;
  Double_t                  MCP_BD_q1;
  Double_t                  MCP_SIG_q1;
  
  Double_t					X0;
  Double_t					Y0;
   
    
  //root & data files
  TString        fRootName;
  TString		 fDataName;
  TString		 RootFolder = "./";
  TString		 DataFolder = "./";
  TString    FileName;

  TFile*         root_file;

  FileName = argv[1];
  fRootName = RootFolder += FileName + ".root";
  fDataName = DataFolder += FileName + ".fast";

  //////////////////////////////////////////////////////////////////////
  //  output root file
  root_file = new TFile (fRootName.Data (), "recreate");
  
  // output root trees & branches
  //goup_tree
  TTree *treeMCP = new TTree ("treeMCP", "treeMCP");
  
  treeMCP->Branch ("MCP_HG_q1", &MCP_HG_q1,  "MCP_HG_q1/D");
  treeMCP->Branch ("MCP_HD_q1", &MCP_HD_q1,  "MCP_HD_q1/D");
  treeMCP->Branch ("MCP_BG_q1", &MCP_BG_q1,  "MCP_BG_q1/D");
  treeMCP->Branch ("MCP_BD_q1", &MCP_BD_q1,  "MCP_BD_q1/D");
  treeMCP->Branch ("MCP_SIG_q1", &MCP_SIG_q1, "MCP_SIG_q1/D");
  treeMCP->Branch ("X0", &X0, "X0/D");
  treeMCP->Branch ("Y0", &Y0, "Y0/D");
   
   
  //  output histograms
  TH2F *h_Image = new TH2F("h_Image","h_Image",500,-1.,1.,500,-1.,1.);
  TH1F *h_MCP_HG_q1 = new TH1F("h_MCP_HG_q1","h_MCP_HG_q1",1000,0.,100000);
  TH1F *h_MCP_HD_q1 = new TH1F("h_MCP_HD_q1","h_MCP_HD_q1",1000,0.,100000);
  TH1F *h_MCP_BG_q1 = new TH1F("h_MCP_BG_q1","h_MCP_BG_q1",1000,0.,100000);
  TH1F *h_MCP_BD_q1 = new TH1F("h_MCP_BD_q1","h_MCP_BD_q1",1000,0.,100000);
  TH1F *h_MCP_SIG_q1 = new TH1F("h_MCP_SIG_q1","h_MCP_SIG_q1",1000,0.,100000);
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//  loop on runs to read
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


 //////////////////////////////////////////////////////////////////////
 //open faster file reader
  reader = faster_file_reader_open (fDataName.Data());
  if (reader == NULL) {
    printf ("error opening file %s\n", fDataName.Data());
    return EXIT_FAILURE;}
  
  printf ("     - read faster file '%s'\n", fDataName.Data());
  

  //  print infos
  printf ("\n");
  printf ("  group2tree_XF_MCP :\n");
  printf ("     - get data labels  (grouped or not)\n");
  printf ("     - output data to root file '%s'\n", fRootName.Data());
  printf ("         ->  special value '%f' means no data\n", NO_DATA);
  printf ("         ->  dont forget to cut this value when histogramming\n");
  printf ("\n");

int i_group=0;
int i_no_group=0;
///////////////////////////////////////////////////////////////////////
/////////////////////////////main loop//////////////////////////////////
///////////////////////////////////////////////////////////////////////
  // main loop : for each data of the file
  while ((data = faster_file_reader_next (reader)) != NULL)
  {
    //  reset leaves
    
     
    MCP_HG_q1  = NO_DATA;
    MCP_HD_q1  = NO_DATA;
    MCP_BG_q1  = NO_DATA;
    MCP_BD_q1  = NO_DATA;
	  MCP_SIG_q1 = NO_DATA;
	  X0		     = NO_DATA;
	  Y0		     = NO_DATA;
        
    //  get type and label of the data
    alias = faster_data_type_alias (data);
    label = faster_data_label      (data);
    clock_ns = faster_data_clock_ns (data);
    clock_sec = faster_data_clock_sec (data);
    
///////////////////////////////////////////////////////////////////////
/////////////////////////////if within group//////////////////////////////////
///////////////////////////////////////////////////////////////////////
    
    //  test its type
    if (alias == GROUP_TYPE_ALIAS)
    { //  it's a group of data (ie : coincidence)
      //  get the specific GROUP part of the data
      lsize = faster_data_load (data, group_buffer);
      //  open a reader for that group
      group_reader = faster_buffer_reader_open (group_buffer, lsize);
      //  loop on each data of the group
	  i_group++;
	  if (i_group%100000==0) cout<<"group "<<i_group<<endl;
	  
      while ((group_data = faster_buffer_reader_next (group_reader)) != NULL)
      {
        //  get type & label & clock of the grouped data
        alias = faster_data_type_alias (group_data);
        label = faster_data_label      (group_data);
		
        	
	  if (alias == QDC_X2_TYPE_ALIAS)
    { //  it's a QDC x 2
	    //  get the QDC specific part of the data
      faster_data_load (group_data, &qdc2);
      //  test its label to assign its charge to the corresponding leaf
            
			
			if (label == LABEL_MCP_SIG)
			{
			   MCP_SIG_q1 = qdc2.q1;
        if ( MCP_SIG_q1 > MCP_q1_min)
			    h_MCP_SIG_q1->Fill(MCP_SIG_q1);
			}
    }
	
	  if (alias == CRRC4_SPECTRO_TYPE_ALIAS)
    { //  it's a CRRC4
      //  get the QDC specific part of the data
      faster_data_load (group_data, &spectro_data);
      //  test its label to assign its charge to the corresponding leaf
      if (label == LABEL_MCP_HG)
			{
			MCP_HG_q1 = spectro_data.measure;
			h_MCP_HG_q1->Fill(MCP_HG_q1);
			}
			
			if (label == LABEL_MCP_HD)
			{
			MCP_HD_q1 = spectro_data.measure;
			h_MCP_HD_q1->Fill(MCP_HD_q1);
			}
			
			if (label == LABEL_MCP_BG)
			{
			MCP_BG_q1 = spectro_data.measure;
			h_MCP_BG_q1->Fill(MCP_BG_q1);
			}
			
			if (label == LABEL_MCP_BD)
			{
			MCP_BD_q1 = spectro_data.measure;
			h_MCP_BD_q1->Fill(MCP_BD_q1);
			}

    }
	}// end of loop on data while in groups
	
	
	X0 = (MCP_HD_q1 + MCP_BD_q1 - MCP_HG_q1 - MCP_BG_q1)/(MCP_HD_q1 + MCP_BD_q1 + MCP_HG_q1 + MCP_BG_q1);
	Y0 = (MCP_HD_q1 + MCP_HG_q1 - MCP_BG_q1 - MCP_BD_q1)/(MCP_HD_q1 + MCP_BD_q1 + MCP_HG_q1 + MCP_BG_q1);
	
	if (MCP_SIG_q1 > MCP_q1_min && MCP_HD_q1 != NO_DATA && MCP_HG_q1 != NO_DATA && MCP_BD_q1 != NO_DATA && MCP_BG_q1 != NO_DATA) {h_Image->Fill(X0,Y0); treeMCP->Fill();}
	
	
	// close the group reader
      faster_buffer_reader_close (group_reader);     
    }// end (if group)
	
    else
		{
		i_no_group++;
	    if (i_no_group%100000==0) cout<<endl<<"no group "<<i_no_group<<endl;
	  }
      
    }//end (main loop)


  //  close the faster file and quit
  faster_file_reader_close (reader);
  root_file->Write ();
  root_file->Close ();

  

  return EXIT_SUCCESS;

}




