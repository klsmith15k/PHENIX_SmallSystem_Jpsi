
#include <TF1.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <TVirtualFitter.h>
#include <TMatrixDSym.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLine.h>
#include <TRandom1.h>
#include <TPolyLine.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TNtuple.h>
#include <TPaveText.h>

using namespace std;

void Draw_RAB_pT_cent()
{
  char name1[500];
  char name2[500];
  char name3[500];

  TFile *file_pAu;	  
  TFile *file_pAl;	  
  TFile *file_HeAu;	  
  
 // RAB arrays for TGraphs
  double RpAu_05_array[2][19] = {0};
  double RpAu_510_array[2][19] = {0};
  double RpAu_1020_array[2][19] = {0};
  double RpAu_020_array[2][19] = {0};
  double RpAu_2040_array[2][19] = {0};
  double RpAu_4060_array[2][19] = {0};
  double RpAu_6084_array[2][19] = {0}; 
  
  double RpAl_020_array[2][16] = {0};
  double RpAl_2040_array[2][16] = {0};
  double RpAl_4072_array[2][16] = {0};
  
  double RHeAu_020_array[2][13] = {0};
  double RHeAu_2040_array[2][13] = {0};
  double RHeAu_4088_array[2][13] = {0};

  // Type A error arrays for TGraphs
  double x_errors[19] = {0};
  // Width of type B Errors
  double width[19] = {0};

  double RpAu_05_ERR_array[2][19] = {0};
  double RpAu_510_ERR_array[2][19] = {0};
  double RpAu_1020_ERR_array[2][19] = {0};
  double RpAu_020_ERR_array[2][19] = {0};
  double RpAu_2040_ERR_array[2][19] = {0};
  double RpAu_4060_ERR_array[2][19] = {0};
  double RpAu_6084_ERR_array[2][19] = {0};
  
  double RpAl_020_ERR_array[2][16] = {0};
  double RpAl_2040_ERR_array[2][16] = {0};
  double RpAl_4072_ERR_array[2][16] = {0};
  
  double RHeAu_020_ERR_array[2][13] = {0};
  double RHeAu_2040_ERR_array[2][13] = {0};
  double RHeAu_4088_ERR_array[2][13] = {0};

 // Type B absolute error arrays for TGraphs
  double RpAu_05_sys_array[2][19] = {0};
  double RpAu_510_sys_array[2][19] = {0};
  double RpAu_1020_sys_array[2][19] = {0};
  double RpAu_020_sys_array[2][19] = {0};
  double RpAu_2040_sys_array[2][19] = {0};
  double RpAu_4060_sys_array[2][19] = {0};
  double RpAu_6084_sys_array[2][19] = {0};
  
  double RpAl_020_sys_array[2][16] = {0};
  double RpAl_2040_sys_array[2][16] = {0};
  double RpAl_4072_sys_array[2][16] = {0};
  
  double RHeAu_020_sys_array[2][13] = {0};
  double RHeAu_2040_sys_array[2][13] = {0};
  double RHeAu_4088_sys_array[2][13] = {0};

  // RAB
  TH1F *RpAu_N_05 = new TH1F("pAu_05_N","pAu_05_N",19, 0, 7);
  TH1F *RpAu_N_510 = new TH1F("pAu_510_N","pAu_510_N",19, 0, 7);
  TH1F *RpAu_N_1020 = new TH1F("pAu_1020_N","pAu_1020_N",19, 0, 7); 
  TH1F *RpAu_N_020 = new TH1F("pAu_020_N","pAu_020_N",19, 0, 7);
  TH1F *RpAu_N_2040 =  new TH1F("pAu_2040_N","pAu_2040_N",19, 0, 7);
  TH1F *RpAu_N_4060 = new TH1F("pAu_4060_N","pAu_4060_N",19, 0, 7);
  TH1F *RpAu_N_6084 = new TH1F("pAu_6084_N","pAu_6084_N",19, 0, 7);

  TH1F *RpAu_S_05 = new TH1F("pAu_05_S","pAu_05_S",19, 0, 7);
  TH1F *RpAu_S_510 = new TH1F("pAu_510_S","pAu_510_S",19, 0, 7);
  TH1F *RpAu_S_1020 = new TH1F("pAu_1020_S","pAu_1020_S",19, 0, 7); 
  TH1F *RpAu_S_020 = new TH1F("pAu_020_S","pAu_020_S",19, 0, 7);
  TH1F *RpAu_S_2040 =  new TH1F("pAu_2040_S","pAu_2040_S",19, 0, 7);
  TH1F *RpAu_S_4060 = new TH1F("pAu_4060_S","pAu_4060_S",19, 0, 7);
  TH1F *RpAu_S_6084 = new TH1F("pAu_6084_S","pAu_6084_S",19, 0, 7);

  TH1F *RpAl_N_020 = new TH1F("RpAl_N_020","RpAl_N_020",16, 0, 4);
  TH1F *RpAl_N_2040 = new TH1F("pAl_2040_N","pAl_2040_N",16, 0, 4);
  TH1F *RpAl_N_4072 = new TH1F("pAl_4072_N","pAl_4072_N",16, 0, 4);
 
  TH1F *RpAl_S_020 = new TH1F("RpAl_S_020","RpAl_S_020",16, 0, 4);
  TH1F *RpAl_S_2040 = new TH1F("pAl_2040_S","pAl_2040_S",16, 0, 4);
  TH1F *RpAl_S_4072 = new TH1F("pAl_4072_S","pAl_4072_S",16, 0, 4);
    
  TH1F *RHeAu_N_020 = new TH1F("HeAu_020_N","HeAu_020_N",16, 0, 4);
  TH1F *RHeAu_N_2040 = new TH1F("HeAu_2040_N","HeAu_2040_N",16, 0, 4);
  TH1F *RHeAu_N_4088 = new TH1F("HeAu_4072_N","HeAu_4072_N",16, 0, 4);

  TH1F *RHeAu_S_020 = new TH1F("HeAu_020_S","HeAu_020_S",16, 0, 4);
  TH1F *RHeAu_S_2040 = new TH1F("HeAu_2040_S","HeAu_2040_S",16, 0, 4);
  TH1F *RHeAu_S_4088 = new TH1F("HeAu_4072_S","HeAu_4072_S",16, 0, 4);

 // RAB Fractional Systematic Error

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *RpAu_05_N_SYS_FRAC_ERR = new TH1F("pAu_05_N_sys","pAu_05_N_sys",19, 0, 7);
  TH1F *RpAu_510_N_SYS_FRAC_ERR = new TH1F("pAu_510_N_sys","pAu_510_N_sys",19, 0, 7);
  TH1F *RpAu_1020_N_SYS_FRAC_ERR = new TH1F("pAu_1020_N_sys","pAu_1020_N_sys",19, 0, 7);
  TH1F *RpAu_020_N_SYS_FRAC_ERR = new TH1F("pAu_020_N_sys","pAu_020_N_sys",19, 0, 7);
  TH1F *RpAu_2040_N_SYS_FRAC_ERR = new TH1F("pAu_2040_N_sys","pAu_2040_N_sys",19, 0, 7);
  TH1F *RpAu_4060_N_SYS_FRAC_ERR = new TH1F("pAu_4060_N_sys","pAu_4060_N_sys",19, 0, 7);
  TH1F *RpAu_6084_N_SYS_FRAC_ERR = new TH1F("pAu_6084_N_sys","pAu_6084_N_sys",19, 0, 7);

  TH1F *RpAu_05_S_SYS_FRAC_ERR = new TH1F("pAu_05_S_sys","pAu_05_S_sys",19, 0, 7);
  TH1F *RpAu_510_S_SYS_FRAC_ERR = new TH1F("pAu_510_S_sys","pAu_510_S_sys",19, 0, 7);
  TH1F *RpAu_1020_S_SYS_FRAC_ERR = new TH1F("pAu_1020_S_sys","pAu_1020_S_sys",19, 0, 7);
  TH1F *RpAu_020_S_SYS_FRAC_ERR = new TH1F("pAu_020_S_sys","pAu_020_S_sys",19, 0, 7);
  TH1F *RpAu_2040_S_SYS_FRAC_ERR = new TH1F("pAu_2040_S_sys","pAu_2040_S_sys",19, 0, 7);
  TH1F *RpAu_4060_S_SYS_FRAC_ERR = new TH1F("pAu_4060_S_sys","pAu_4060_S_sys",19, 0, 7);
  TH1F *RpAu_6084_S_SYS_FRAC_ERR = new TH1F("pAu_6084_S_sys","pAu_6084_S_sys",19, 0, 7);

  TH1F *RpAl_020_N_SYS_FRAC_ERR = new TH1F("pAl_020_N_sys","pAl_020_N_sys",16, 0, 4);
  TH1F *RpAl_2040_N_SYS_FRAC_ERR = new TH1F("pAl_2040_N_sys","pAl_2040_N_sys",16, 0, 4);
  TH1F *RpAl_4072_N_SYS_FRAC_ERR = new TH1F("pAl_4072_N_sys","pAl_4072_N_sys",16, 0, 4);

  TH1F *RpAl_020_S_SYS_FRAC_ERR = new TH1F("pAl_020_S_sys","pAl_020_S_sys",16, 0, 4);
  TH1F *RpAl_2040_S_SYS_FRAC_ERR = new TH1F("pAl_2040_S_sys","pAl_2040_S_sys",16, 0, 4);
  TH1F *RpAl_4072_S_SYS_FRAC_ERR = new TH1F("pAl_4072_S_sys","pAl_4072_S_sys",16, 0, 4);
  
  TH1F *RHeAu_020_N_SYS_FRAC_ERR = new TH1F("HeAu_020_N_sys","HeAu_020_N_sys",13, 0, 4);
  TH1F *RHeAu_2040_N_SYS_FRAC_ERR = new TH1F("HeAu_2040_N_sys","HeAu_2040_N_sys",13, 0, 4);
  TH1F *RHeAu_4088_N_SYS_FRAC_ERR = new TH1F("HeAu_4072_N_sys","HeAu_4072_N_sys",13, 0, 4);

  TH1F *RHeAu_020_S_SYS_FRAC_ERR = new TH1F("HeAu_020_S_sys","HeAu_020_S_sys",13, 0, 4);
  TH1F *RHeAu_2040_S_SYS_FRAC_ERR = new TH1F("HeAu_2040_S_sys","HeAu_2040_S_sys",13, 0, 4);
  TH1F *RHeAu_4088_S_SYS_FRAC_ERR = new TH1F("HeAu_4072_S_sys","HeAu_4072_S_sys",13, 0, 4);

  file_pAu = TFile::Open("RpAu_YpAu_vs_pT_cent.root");
  
  // South Arm 
  // RpAu
  file_pAu->GetObject("RpAu_vs_pT_S_05",RpAu_S_05);
  file_pAu->GetObject("RpAu_vs_pT_S_510",RpAu_S_510);
  file_pAu->GetObject("RpAu_vs_pT_S_1020",RpAu_S_1020);
  file_pAu->GetObject("RpAu_vs_pT_S_020",RpAu_S_020);
  file_pAu->GetObject("RpAu_vs_pT_S_2040",RpAu_S_2040);
  file_pAu->GetObject("RpAu_vs_pT_S_4060",RpAu_S_4060);
  file_pAu->GetObject("RpAu_vs_pT_S_6084",RpAu_S_6084);
  
 // North Arm 
  // RpAu
  file_pAu->GetObject("RpAu_vs_pT_N_05",RpAu_N_05);
  file_pAu->GetObject("RpAu_vs_pT_N_510",RpAu_N_510);
  file_pAu->GetObject("RpAu_vs_pT_N_1020",RpAu_N_1020);
  file_pAu->GetObject("RpAu_vs_pT_N_020",RpAu_N_020);
  file_pAu->GetObject("RpAu_vs_pT_N_2040",RpAu_N_2040);
  file_pAu->GetObject("RpAu_vs_pT_N_4060",RpAu_N_4060);
  file_pAu->GetObject("RpAu_vs_pT_N_6084",RpAu_N_6084);
 
 // South Arm 
  // RpAu sys err
  file_pAu->GetObject("RpAu_vs_pT_S_05_sys_frac",RpAu_05_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_510_sys_frac",RpAu_510_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_1020_sys_frac",RpAu_1020_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_020_sys_frac",RpAu_020_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_2040_sys_frac",RpAu_2040_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_4060_sys_frac",RpAu_4060_S_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_S_6084_sys_frac",RpAu_6084_S_SYS_FRAC_ERR);

// North Arm 
  // RpAu sys err
  file_pAu->GetObject("RpAu_vs_pT_N_05_sys_frac",RpAu_05_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_510_sys_frac",RpAu_510_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_1020_sys_frac",RpAu_1020_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_020_sys_frac",RpAu_020_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_2040_sys_frac",RpAu_2040_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_4060_sys_frac",RpAu_4060_N_SYS_FRAC_ERR);
  file_pAu->GetObject("RpAu_vs_pT_N_6084_sys_frac",RpAu_6084_N_SYS_FRAC_ERR);

  // file_pAu->Close();

//////////////////////////////////////////////////////////////////////////////

  file_pAl = TFile::Open("RpAl_YpAl_vs_pT_cent.root");

  if(file_pAl->IsOpen())
    cout << "Hello world" << endl;
  else
    cout << "Not open " << endl;

  // // South arm RpAl
  file_pAl->GetObject("RpAl_vs_pT_S_020",RpAl_S_020);
  file_pAl->GetObject("RpAl_vs_pT_S_2040",RpAl_S_2040);
  file_pAl->GetObject("RpAl_vs_pT_S_4072",RpAl_S_4072);
 
  // North arm 
  // RpAl
  file_pAl->GetObject("RpAl_vs_pT_N_020",RpAl_N_020);
  file_pAl->GetObject("RpAl_vs_pT_N_2040",RpAl_N_2040);;
  file_pAl->GetObject("RpAl_vs_pT_N_4072",RpAl_N_4072);

  // South arm 
  // RpAl sys err
  file_pAl->GetObject("RpAl_vs_pT_S_020_sys_frac",RpAl_020_S_SYS_FRAC_ERR);
  file_pAl->GetObject("RpAl_vs_pT_S_2040_sys_frac",RpAl_2040_S_SYS_FRAC_ERR);
  file_pAl->GetObject("RpAl_vs_pT_S_4072_sys_frac",RpAl_4072_S_SYS_FRAC_ERR);

  // North arm 
  // RpAl sys err
  file_pAl->GetObject("RpAl_vs_pT_N_020_sys_frac",RpAl_020_N_SYS_FRAC_ERR);
  file_pAl->GetObject("RpAl_vs_pT_N_2040_sys_frac",RpAl_2040_N_SYS_FRAC_ERR);
  file_pAl->GetObject("RpAl_vs_pT_N_4072_sys_frac",RpAl_4072_N_SYS_FRAC_ERR);

  // file_pAl->Close();
 
  ////////////////////////////////////////////////////////////////////////////////
 
  file_HeAu = TFile::Open("RHeAu_YHeAu_vs_pT_cent.root");
 
  //South arm y
  // RHeAu
  file_HeAu->GetObject("RHeAu_vs_pT_S_020",RHeAu_S_020);
  file_HeAu->GetObject("RHeAu_vs_pT_S_2040",RHeAu_S_2040);
  file_HeAu->GetObject("RHeAu_vs_pT_S_4088",RHeAu_S_4088);

  //North arm 
  // RHeAu
  file_HeAu->GetObject("RHeAu_vs_pT_N_020",RHeAu_N_020);
  file_HeAu->GetObject("RHeAu_vs_pT_N_2040",RHeAu_N_2040);
  file_HeAu->GetObject("RHeAu_vs_pT_N_4088",RHeAu_N_4088);

  //South arm 
  // RHeAu sys err
  file_HeAu->GetObject("RHeAu_vs_pT_S_020_sys_frac",RHeAu_020_S_SYS_FRAC_ERR);
  file_HeAu->GetObject("RHeAu_vs_pT_S_2040_sys_frac",RHeAu_2040_S_SYS_FRAC_ERR);
  file_HeAu->GetObject("RHeAu_vs_pT_S_4088_sys_frac",RHeAu_4088_S_SYS_FRAC_ERR);

  //North arm 
  // RHeAu sys err
  file_HeAu->GetObject("RHeAu_vs_pT_N_020_sys_frac",RHeAu_020_N_SYS_FRAC_ERR);
  file_HeAu->GetObject("RHeAu_vs_pT_N_2040_sys_frac",RHeAu_2040_N_SYS_FRAC_ERR);
  file_HeAu->GetObject("RHeAu_vs_pT_N_4088_sys_frac",RHeAu_4088_N_SYS_FRAC_ERR);	      
  
  // file_HeAu->Close();
  
  // Type C error 
  double err_pAu_c[8] = {11.82, 11.69, 12.12, 11.77, 12, 12.14, 13.97}; // 0-5,5-10,10-20,0-20,20-40,40-60,60-84,0-84
  double err_HeAu_c[4] = {12.62, 12.50, 13.28, 12.11}; // 0-20, 20-40, 40-88, 0-88
  double err_dAu_c[4] = {9.0,8.5,9.1,10.6}; // 0-20, 20-40, 40-60, 60-88
  double err_pAl_c[4] = {12.51, 12.26, 12.62, 11.35}; // 0-20, 20-40, 40-72, 0-72
  
  double temp1[2][19] = {0};
  double temp2[2][19] = {0};
  double temp3[2][19] = {0};
  double temp4[2][19] = {0};
  double temp5[2][19] = {0};
  double temp6[2][19] = {0};
  double temp7[2][19] = {0};
  double temp8[2][19] = {0};
  double temp9[2][19] = {0};
  double temp10[2][19] = {0};
  double temp11[2][19] = {0};
  double temp12[2][19] = {0};
  double temp13[2][19] = {0};
  double temp14[2][19] = {0};

  for(int arm = 0; arm < 2; arm++)
    {  
      for(int i = 0; i < 16; i++)
	{
	  cout << RpAl_S_020->GetBinContent(i+1) << endl;

	  if(arm == 0)
	    {
	      //Fill RAB
	      RpAl_020_array[arm][i] = RpAl_S_020->GetBinContent(i+1);
	      RpAl_2040_array[arm][i] = RpAl_S_2040->GetBinContent(i+1);
	      RpAl_4072_array[arm][i] = RpAl_S_4072->GetBinContent(i+1);
	      // FIll Type A
	      RpAl_020_ERR_array[arm][i] = RpAl_S_020->GetBinError(i+1);
	      RpAl_2040_ERR_array[arm][i] = RpAl_S_2040->GetBinError(i+1);
	      RpAl_4072_ERR_array[arm][i] = RpAl_S_4072->GetBinError(i+1);
	      // Fill type B
	      temp1[arm][i] = RpAl_020_S_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp2[arm][i] = RpAl_2040_S_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp3[arm][i] = RpAl_4072_S_SYS_FRAC_ERR ->GetBinError(i+1);

	      RpAl_020_sys_array[arm][i] = temp1[arm][i]*RpAl_020_array[arm][i];
	      RpAl_2040_sys_array[arm][i] = temp2[arm][i]*RpAl_2040_array[arm][i];
	      RpAl_4072_sys_array[arm][i] = temp3[arm][i]*RpAl_4072_array[arm][i];
	    }
	  if(arm == 1)
	    {
	      // Fill RAB
	      RpAl_020_array[arm][i] = RpAl_N_020->GetBinContent(i+1);
	      RpAl_2040_array[arm][i] = RpAl_N_2040->GetBinContent(i+1);
	      RpAl_4072_array[arm][i] = RpAl_N_4072->GetBinContent(i+1);
	      // FIll Type A
	      RpAl_020_ERR_array[arm][i] = RpAl_N_020->GetBinError(i+1);
	      RpAl_2040_ERR_array[arm][i] = RpAl_N_2040->GetBinError(i+1);
	      RpAl_4072_ERR_array[arm][i] = RpAl_N_4072->GetBinError(i+1);
	      // Fill type B
	      temp4[arm][i] = RpAl_020_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp5[arm][i] = RpAl_2040_N_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp6[arm][i] = RpAl_4072_N_SYS_FRAC_ERR ->GetBinError(i+1);

	      RpAl_020_sys_array[arm][i] = temp4[arm][i]*RpAl_020_array[arm][i];
	      RpAl_2040_sys_array[arm][i] = temp5[arm][i]*RpAl_2040_array[arm][i];
	      RpAl_4072_sys_array[arm][i] = temp6[arm][i]*RpAl_4072_array[arm][i];
	    }
	}
      
    }

 // for(int arm = 0; arm < 2; arm++)
 //    {  
 //      for(int i = 0; i < 16; i++)
 // 	{
 // 	  cout << RpAl_S_020->GetBinContent(i+1) << endl;

 // 	  if(arm == 0)
 // 	    {
 // 	      RpAl_020_sys_array[arm][i] = temp1[arm][i]*RpAl_020_array[arm][i];
 // 	      RpAl_2040_sys_array[arm][i] = temp2[arm][i]*RpAl_2040_array[arm][i];
 // 	      RpAl_4072_sys_array[arm][i] = temp3[arm][i]*RpAl_4072_array[arm][i];
 // 	    }
 // 	  if(arm == 1)
 // 	    {
 // 	      RpAl_020_sys_array[arm][i] = temp4[arm][i]*RpAl_020_array[arm][i];
 // 	      RpAl_2040_sys_array[arm][i] = temp5[arm][i]*RpAl_2040_array[arm][i];
 // 	      RpAl_4072_sys_array[arm][i] = temp6[arm][i]*RpAl_4072_array[arm][i];
 // 	    }
 // 	}
      
 //    }





  
  temp1[2][19] = {0};
  temp2[2][19] = {0};
  temp3[2][19] = {0};
  temp4[2][19] = {0};
  temp5[2][19] = {0};
  temp6[2][19] = {0};
  
  // extract all bincontents and store in arrays to plot any system and centrality
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < 13; i++)
	{
	  if(arm == 0)
	    {
	      //Fill RAB
	      RHeAu_020_array[arm][i] = RHeAu_S_020->GetBinContent(i+1);
	      RHeAu_2040_array[arm][i] = RHeAu_S_2040->GetBinContent(i+1);
	      RHeAu_4088_array[arm][i] = RHeAu_S_4088->GetBinContent(i+1);
	      // FIll Type A
	      RHeAu_020_ERR_array[arm][i] = RHeAu_S_020->GetBinError(i+1);
	      RHeAu_2040_ERR_array[arm][i] = RHeAu_S_2040->GetBinError(i+1);
	      RHeAu_4088_ERR_array[arm][i] = RHeAu_S_4088->GetBinError(i+1);
	      // Fill type B
	      temp1[arm][i] = RHeAu_020_S_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp2[arm][i] = RHeAu_2040_S_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp3[arm][i] = RHeAu_4088_S_SYS_FRAC_ERR ->GetBinError(i+1);

	      RHeAu_020_sys_array[arm][i] = temp1[arm][i]*RHeAu_020_array[arm][i];
	      RHeAu_2040_sys_array[arm][i] = temp2[arm][i]*RHeAu_2040_array[arm][i];
	      RHeAu_4088_sys_array[arm][i] = temp3[arm][i]*RHeAu_4088_array[arm][i];
	    }
	  if(arm == 1)
	    {
	      // Fill RAB
	      RHeAu_020_array[arm][i] = RHeAu_N_020->GetBinContent(i+1);
	      RHeAu_2040_array[arm][i] = RHeAu_N_2040->GetBinContent(i+1);
	      RHeAu_4088_array[arm][i] = RHeAu_N_4088->GetBinContent(i+1);
	      // FIll Type A
	      RHeAu_020_ERR_array[arm][i] = RHeAu_N_020->GetBinError(i+1);
	      RHeAu_2040_ERR_array[arm][i] = RHeAu_N_2040->GetBinError(i+1);
	      RHeAu_4088_ERR_array[arm][i] = RHeAu_N_4088->GetBinError(i+1);
	      // Fill type B
	      temp4[arm][i] = RHeAu_020_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp5[arm][i] = RHeAu_2040_N_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp6[arm][i] = RHeAu_4088_N_SYS_FRAC_ERR ->GetBinError(i+1);

	      RHeAu_020_sys_array[arm][i] = temp4[arm][i]*RHeAu_020_array[arm][i];
	      RHeAu_2040_sys_array[arm][i] = temp5[arm][i]*RHeAu_2040_array[arm][i];
	      RHeAu_4088_sys_array[arm][i] = temp6[arm][i]*RHeAu_4088_array[arm][i];
	    }
	}
    }

  temp1[2][19] = {0};
  temp2[2][19] = {0};
  temp3[2][19] = {0};
  temp4[2][19] = {0};
  temp5[2][19] = {0};
  temp6[2][19] = {0};  

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < 19; i++)
	{
	  if(arm == 0)
	    {
	      //Fill RAB
	      RpAu_05_array[arm][i] = RpAu_S_05->GetBinContent(i+1);
	      RpAu_510_array[arm][i] = RpAu_S_510->GetBinContent(i+1);
	      RpAu_1020_array[arm][i] = RpAu_S_1020->GetBinContent(i+1);
	      RpAu_020_array[arm][i] = RpAu_S_020->GetBinContent(i+1);
	      RpAu_2040_array[arm][i] = RpAu_S_2040->GetBinContent(i+1);
	      RpAu_4060_array[arm][i] = RpAu_S_4060->GetBinContent(i+1);
	      RpAu_6084_array[arm][i] = RpAu_S_6084->GetBinContent(i+1);
	      // FIll Type A
	      RpAu_05_ERR_array[arm][i] = RpAu_S_05->GetBinError(i+1);
	      RpAu_510_ERR_array[arm][i] = RpAu_S_510->GetBinError(i+1);
	      RpAu_1020_ERR_array[arm][i] = RpAu_S_1020->GetBinError(i+1);
	      RpAu_020_ERR_array[arm][i] = RpAu_S_020->GetBinError(i+1);
	      RpAu_2040_ERR_array[arm][i] = RpAu_S_2040->GetBinError(i+1);
	      RpAu_4060_ERR_array[arm][i] = RpAu_S_4060->GetBinError(i+1);
	      RpAu_6084_ERR_array[arm][i] = RpAu_S_6084->GetBinError(i+1);
	      // Fill type B
	      temp1[arm][i] = RpAu_05_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp2[arm][i] = RpAu_510_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp3[arm][i] = RpAu_1020_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp4[arm][i] = RpAu_020_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp5[arm][i] = RpAu_2040_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp6[arm][i] = RpAu_4060_S_SYS_FRAC_ERR->GetBinError(i+1);
	      temp7[arm][i] = RpAu_6084_S_SYS_FRAC_ERR->GetBinError(i+1);

	      RpAu_05_sys_array[arm][i] = temp1[arm][i]*RpAu_05_array[arm][i];
	      RpAu_510_sys_array[arm][i] = temp2[arm][i]*RpAu_510_array[arm][i];
	      RpAu_1020_sys_array[arm][i] = temp3[arm][i]*RpAu_1020_array[arm][i];
	      RpAu_020_sys_array[arm][i] = temp4[arm][i]*RpAu_020_array[arm][i];
	      RpAu_2040_sys_array[arm][i] = temp5[arm][i]*RpAu_2040_array[arm][i];
	      RpAu_4060_sys_array[arm][i] = temp6[arm][i]*RpAu_4060_array[arm][i];
	      RpAu_6084_sys_array[arm][i] = temp7[arm][i]*RpAu_6084_array[arm][i];
	    }
	   if(arm == 1)
	    {
	        // Fill RAB
	      RpAu_05_array[arm][i] = RpAu_N_05->GetBinContent(i+1);
	      RpAu_510_array[arm][i] = RpAu_N_510->GetBinContent(i+1);
	      RpAu_1020_array[arm][i] = RpAu_N_1020->GetBinContent(i+1);
	      RpAu_020_array[arm][i] = RpAu_N_020->GetBinContent(i+1);
	      RpAu_2040_array[arm][i] = RpAu_N_2040->GetBinContent(i+1);
	      RpAu_4060_array[arm][i] = RpAu_N_4060->GetBinContent(i+1);
	      RpAu_6084_array[arm][i] = RpAu_N_6084->GetBinContent(i+1);
	      // FIll Type A
	      RpAu_05_ERR_array[arm][i] = RpAu_N_05->GetBinError(i+1);
	      RpAu_510_ERR_array[arm][i] = RpAu_N_510->GetBinError(i+1);
	      RpAu_1020_ERR_array[arm][i] = RpAu_N_1020->GetBinError(i+1);
	      RpAu_020_ERR_array[arm][i] = RpAu_N_020->GetBinError(i+1);
	      RpAu_2040_ERR_array[arm][i] = RpAu_N_2040->GetBinError(i+1);
	      RpAu_4060_ERR_array[arm][i] = RpAu_N_4060->GetBinError(i+1);
	      RpAu_6084_ERR_array[arm][i] = RpAu_N_6084->GetBinError(i+1);
	      // Fill type B
	      temp8[arm][i] = RpAu_05_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp9[arm][i] = RpAu_510_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp10[arm][i] = RpAu_1020_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp11[arm][i] = RpAu_020_N_SYS_FRAC_ERR->GetBinError(i+1);
	      temp12[arm][i] = RpAu_2040_N_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp13[arm][i] = RpAu_4060_N_SYS_FRAC_ERR ->GetBinError(i+1);
	      temp14[arm][i] = RpAu_6084_N_SYS_FRAC_ERR ->GetBinError(i+1);

	      RpAu_05_sys_array[arm][i] = temp8[arm][i]*RpAu_05_array[arm][i]; 
	      RpAu_510_sys_array[arm][i] = temp9[arm][i]*RpAu_510_array[arm][i]; 
	      RpAu_1020_sys_array[arm][i] = temp10[arm][i]*RpAu_1020_array[arm][i]; 
	      RpAu_020_sys_array[arm][i] = temp11[arm][i]*RpAu_020_array[arm][i]; 
	      RpAu_2040_sys_array[arm][i] = temp12[arm][i]*RpAu_2040_array[arm][i]; 
	      RpAu_4060_sys_array[arm][i] = temp13[arm][i]*RpAu_4060_array[arm][i]; 
	      RpAu_6084_sys_array[arm][i] = temp14[arm][i]*RpAu_6084_array[arm][i]; 
	    }
	   cout << "RpAu " << RpAu_05_array[arm][i] << endl;
	   width[i] = 0.125;
	}
    }

  // syst uncertainty boxes
  double box_pt = 0.1;
  int box_color = 0; 
  
  //plot arrays
  double pt_array_pAu[19] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0};
  double pt_array_pAl[16] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875};
  double pt_array_HeAu[13] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75};
      
  // define unique TGraphs
  // RAB
  TGraphErrors *gr_pAu_05_N = new TGraphErrors(19,pt_array_pAu, RpAu_05_array[1],x_errors, RpAu_05_ERR_array[1]);
  TGraphErrors *gr_pAu_510_N = new TGraphErrors(19,pt_array_pAu, RpAu_510_array[1],x_errors, RpAu_510_ERR_array[1]);
  TGraphErrors *gr_pAu_1020_N = new TGraphErrors(19,pt_array_pAu, RpAu_1020_array[1],x_errors, RpAu_1020_ERR_array[1]);
  TGraphErrors *gr_pAu_020_N = new TGraphErrors(19,pt_array_pAu, RpAu_020_array[1],x_errors, RpAu_020_ERR_array[1]);
  TGraphErrors *gr_pAu_2040_N = new TGraphErrors(19,pt_array_pAu, RpAu_2040_array[1],x_errors, RpAu_2040_ERR_array[1]);
  TGraphErrors *gr_pAu_4060_N = new TGraphErrors(19,pt_array_pAu, RpAu_4060_array[1],x_errors, RpAu_4060_ERR_array[1]);
  TGraphErrors *gr_pAu_6084_N = new TGraphErrors(19,pt_array_pAu, RpAu_6084_array[1],x_errors, RpAu_6084_ERR_array[1]);

  TGraphErrors *gr_pAu_05_S = new TGraphErrors(19,pt_array_pAu, RpAu_05_array[0],x_errors, RpAu_05_ERR_array[0]);
  TGraphErrors *gr_pAu_510_S = new TGraphErrors(19,pt_array_pAu, RpAu_510_array[0],x_errors, RpAu_510_ERR_array[0]);
  TGraphErrors *gr_pAu_1020_S = new TGraphErrors(19,pt_array_pAu, RpAu_1020_array[0],x_errors, RpAu_1020_ERR_array[0]);
  TGraphErrors *gr_pAu_020_S = new TGraphErrors(19,pt_array_pAu, RpAu_020_array[0],x_errors, RpAu_020_ERR_array[0]);
  TGraphErrors *gr_pAu_2040_S = new TGraphErrors(19,pt_array_pAu, RpAu_2040_array[0],x_errors, RpAu_2040_ERR_array[0]);
  TGraphErrors *gr_pAu_4060_S = new TGraphErrors(19,pt_array_pAu, RpAu_4060_array[0],x_errors, RpAu_4060_ERR_array[0]);
  TGraphErrors *gr_pAu_6084_S = new TGraphErrors(19,pt_array_pAu, RpAu_6084_array[0],x_errors, RpAu_6084_ERR_array[0]);

  TGraphErrors *gr_pAl_020_N = new TGraphErrors(16,pt_array_pAl, RpAl_020_array[1],x_errors,RpAl_020_ERR_array[1]);
  TGraphErrors *gr_pAl_2040_N = new TGraphErrors(16,pt_array_pAl, RpAl_2040_array[1],x_errors,RpAl_2040_ERR_array[1]);
  TGraphErrors *gr_pAl_4072_N = new TGraphErrors(16,pt_array_pAl, RpAl_4072_array[1],x_errors,RpAl_4072_ERR_array[1]);
      
  TGraphErrors *gr_pAl_020_S = new TGraphErrors(16,pt_array_pAl, RpAl_020_array[0],x_errors,RpAl_020_ERR_array[0]);
  TGraphErrors *gr_pAl_2040_S = new TGraphErrors(16,pt_array_pAl, RpAl_2040_array[0],x_errors,RpAl_2040_ERR_array[0]);
  TGraphErrors *gr_pAl_4072_S = new TGraphErrors(16,pt_array_pAl, RpAl_4072_array[0],x_errors,RpAl_4072_ERR_array[0]);
 
  TGraphErrors *gr_HeAu_020_N = new TGraphErrors(13,pt_array_HeAu, RHeAu_020_array[1], x_errors,RHeAu_020_ERR_array[1]);
  TGraphErrors *gr_HeAu_2040_N = new TGraphErrors(13,pt_array_HeAu, RHeAu_2040_array[1], x_errors,RHeAu_2040_ERR_array[1]);
  TGraphErrors *gr_HeAu_4088_N = new TGraphErrors(13,pt_array_HeAu, RHeAu_4088_array[1], x_errors,RHeAu_4088_ERR_array[1]);

  TGraphErrors *gr_HeAu_020_S = new TGraphErrors(13,pt_array_HeAu, RHeAu_020_array[0], x_errors,RHeAu_020_ERR_array[0]);
  TGraphErrors *gr_HeAu_2040_S = new TGraphErrors(13,pt_array_HeAu, RHeAu_2040_array[0], x_errors,RHeAu_2040_ERR_array[0]);
  TGraphErrors *gr_HeAu_4088_S = new TGraphErrors(13,pt_array_HeAu, RHeAu_4088_array[0], x_errors,RHeAu_4088_ERR_array[0]);
 
  // Type B

  TGraphErrors *gr_pAu_05_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_05_array[1],width,RpAu_05_sys_array[1]);
  TGraphErrors *gr_pAu_510_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_510_array[1],width,RpAu_510_sys_array[1]);
  TGraphErrors *gr_pAu_1020_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_1020_array[1],width, RpAu_1020_sys_array[1]);
  TGraphErrors *gr_pAu_020_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_020_array[1],width, RpAu_020_sys_array[1]);
  TGraphErrors *gr_pAu_2040_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_2040_array[1],width, RpAu_2040_sys_array[1]);
  TGraphErrors *gr_pAu_4060_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_4060_array[1],width, RpAu_4060_sys_array[1]);
  TGraphErrors *gr_pAu_6084_N_sys = new TGraphErrors(19,pt_array_pAu, RpAu_6084_array[1],width, RpAu_6084_sys_array[1]);

  TGraphErrors *gr_pAu_05_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_05_array[0],width, RpAu_05_sys_array[1]);
  TGraphErrors *gr_pAu_510_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_510_array[0],width, RpAu_510_sys_array[1]);
  TGraphErrors *gr_pAu_1020_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_1020_array[0],width, RpAu_1020_sys_array[1]);
  TGraphErrors *gr_pAu_020_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_020_array[0],width, RpAu_020_sys_array[1]);
  TGraphErrors *gr_pAu_2040_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_2040_array[0],width, RpAu_2040_sys_array[1]);
  TGraphErrors *gr_pAu_4060_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_4060_array[0],width, RpAu_4060_sys_array[1]);
  TGraphErrors *gr_pAu_6084_S_sys = new TGraphErrors(19,pt_array_pAu, RpAu_6084_array[0],width, RpAu_6084_sys_array[1]);

  TGraphErrors *gr_pAl_020_N_sys = new TGraphErrors(16,pt_array_pAl, RpAl_020_array[1],width,RpAl_020_sys_array[1]);
  TGraphErrors *gr_pAl_2040_N_sys = new TGraphErrors(16,pt_array_pAl, RpAl_2040_array[1],width,RpAl_2040_sys_array[1]);
  TGraphErrors *gr_pAl_4072_N_sys = new TGraphErrors(16,pt_array_pAl, RpAl_4072_array[1],width,RpAl_4072_sys_array[1]);
      
  TGraphErrors *gr_pAl_020_S_sys = new TGraphErrors(16,pt_array_pAl, RpAl_020_array[0],width,RpAl_020_sys_array[1]);
  TGraphErrors *gr_pAl_2040_S_sys = new TGraphErrors(16,pt_array_pAl, RpAl_2040_array[0],width,RpAl_2040_sys_array[1]);
  TGraphErrors *gr_pAl_4072_S_sys = new TGraphErrors(16,pt_array_pAl, RpAl_4072_array[0],width,RpAl_4072_sys_array[1]);
 
  TGraphErrors *gr_HeAu_020_N_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_020_array[1], width,RHeAu_020_sys_array[1]);
  TGraphErrors *gr_HeAu_2040_N_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_2040_array[1], width,RHeAu_2040_sys_array[1]);
  TGraphErrors *gr_HeAu_4088_N_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_4088_array[1], width,RHeAu_4088_sys_array[1]);

  TGraphErrors *gr_HeAu_020_S_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_020_array[0], width,RHeAu_020_sys_array[1]);
  TGraphErrors *gr_HeAu_2040_S_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_2040_array[0], width,RHeAu_2040_sys_array[1]);
  TGraphErrors *gr_HeAu_4088_S_sys = new TGraphErrors(13,pt_array_HeAu, RHeAu_4088_array[0], width,RHeAu_4088_sys_array[1]);
   
  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
  c1->SetGrid();
     
  // set marker style RAB
  gr_pAu_05_N->SetMarkerStyle(20);  
  gr_pAu_510_N->SetMarkerStyle(20);  
  gr_pAu_1020_N->SetMarkerStyle(20);  
  gr_pAu_020_N->SetMarkerStyle(20);  
  gr_pAu_2040_N->SetMarkerStyle(20);  
  gr_pAu_4060_N->SetMarkerStyle(20);  
  gr_pAu_6084_N->SetMarkerStyle(20);  

  gr_pAu_05_S->SetMarkerStyle(20);  
  gr_pAu_510_S->SetMarkerStyle(20);  
  gr_pAu_1020_S->SetMarkerStyle(20);  
  gr_pAu_020_S->SetMarkerStyle(20);  
  gr_pAu_2040_S->SetMarkerStyle(20);  
  gr_pAu_4060_S->SetMarkerStyle(20);  
  gr_pAu_6084_S->SetMarkerStyle(20);  

  gr_pAl_020_N->SetMarkerStyle(20); 
  gr_pAl_2040_N->SetMarkerStyle(20);  
  gr_pAl_4072_N->SetMarkerStyle(20);  
      
  gr_pAl_020_S->SetMarkerStyle(20);   
  gr_pAl_2040_S->SetMarkerStyle(20);  
  gr_pAl_4072_S->SetMarkerStyle(20);  
 
  gr_HeAu_020_N->SetMarkerStyle(20); 
  gr_HeAu_2040_N->SetMarkerStyle(20);  
  gr_HeAu_4088_N ->SetMarkerStyle(20); 

  gr_HeAu_020_S->SetMarkerStyle(20); 
  gr_HeAu_2040_S->SetMarkerStyle(20);  
  gr_HeAu_4088_S->SetMarkerStyle(20);  
  ////////////////////////  set marker style Type B
  gr_pAu_05_N_sys->SetMarkerStyle(20);  
  gr_pAu_510_N_sys->SetMarkerStyle(20);  
  gr_pAu_1020_N_sys->SetMarkerStyle(20);  
  gr_pAu_020_N_sys->SetMarkerStyle(20);  
  gr_pAu_2040_N_sys->SetMarkerStyle(20);  
  gr_pAu_4060_N_sys->SetMarkerStyle(20);  
  gr_pAu_6084_N_sys->SetMarkerStyle(20);  

  gr_pAu_05_S_sys->SetMarkerStyle(20);  
  gr_pAu_510_S_sys->SetMarkerStyle(20);  
  gr_pAu_1020_S_sys->SetMarkerStyle(20);  
  gr_pAu_020_S_sys->SetMarkerStyle(20);  
  gr_pAu_2040_S_sys->SetMarkerStyle(20);  
  gr_pAu_4060_S_sys->SetMarkerStyle(20);  
  gr_pAu_6084_S_sys->SetMarkerStyle(20);  

  gr_pAl_020_N_sys->SetMarkerStyle(20); 
  gr_pAl_2040_N_sys->SetMarkerStyle(20);  
  gr_pAl_4072_N_sys->SetMarkerStyle(20);  
      
  gr_pAl_020_S_sys->SetMarkerStyle(20);   
  gr_pAl_2040_S_sys->SetMarkerStyle(20);  
  gr_pAl_4072_S_sys->SetMarkerStyle(20);  
 
  gr_HeAu_020_N_sys->SetMarkerStyle(20); 
  gr_HeAu_2040_N_sys->SetMarkerStyle(20);  
  gr_HeAu_4088_N_sys ->SetMarkerStyle(20); 

  gr_HeAu_020_S_sys->SetMarkerStyle(20); 
  gr_HeAu_2040_S_sys->SetMarkerStyle(20);  
  gr_HeAu_4088_S_sys->SetMarkerStyle(20);  

  // set marker size RAB
  gr_pAu_05_N->SetMarkerSize(1.0);  
  gr_pAu_510_N->SetMarkerSize(1.0);  
  gr_pAu_1020_N->SetMarkerSize(1.0);  
  gr_pAu_020_N->SetMarkerSize(1.0);  
  gr_pAu_2040_N->SetMarkerSize(1.0);  
  gr_pAu_4060_N->SetMarkerSize(1.0);  
  gr_pAu_6084_N->SetMarkerSize(1.0);  

  gr_pAu_05_S->SetMarkerSize(1.0);  
  gr_pAu_510_S->SetMarkerSize(1.0);  
  gr_pAu_1020_S->SetMarkerSize(1.0);  
  gr_pAu_020_S->SetMarkerSize(1.0);  
  gr_pAu_2040_S->SetMarkerSize(1.0);  
  gr_pAu_4060_S->SetMarkerSize(1.0);  
  gr_pAu_6084_S->SetMarkerSize(1.0);  

  gr_pAl_020_N->SetMarkerSize(1.0); 
  gr_pAl_2040_N->SetMarkerSize(1.0);  
  gr_pAl_4072_N->SetMarkerSize(1.0);  
      
  gr_pAl_020_S->SetMarkerSize(1.0);   
  gr_pAl_2040_S->SetMarkerSize(1.0);  
  gr_pAl_4072_S->SetMarkerSize(1.0);  
 
  gr_HeAu_020_N->SetMarkerSize(1.0); 
  gr_HeAu_2040_N->SetMarkerSize(1.0);  
  gr_HeAu_4088_N ->SetMarkerSize(1.0); 

  gr_HeAu_020_S->SetMarkerSize(1.0); 
  gr_HeAu_2040_S->SetMarkerSize(1.0);  
  gr_HeAu_4088_S->SetMarkerSize(1.0); 
  /////////////////// Set Marker Size Type B
  gr_pAu_05_N_sys->SetMarkerSize(1.0);  
  gr_pAu_510_N_sys->SetMarkerSize(1.0);  
  gr_pAu_1020_N_sys->SetMarkerSize(1.0);  
  gr_pAu_020_N_sys->SetMarkerSize(1.0);  
  gr_pAu_2040_N_sys->SetMarkerSize(1.0);  
  gr_pAu_4060_N_sys->SetMarkerSize(1.0);  
  gr_pAu_6084_N_sys->SetMarkerSize(1.0);  

  gr_pAu_05_S_sys->SetMarkerSize(1.0);  
  gr_pAu_510_S_sys->SetMarkerSize(1.0);  
  gr_pAu_1020_S_sys->SetMarkerSize(1.0);  
  gr_pAu_020_S_sys->SetMarkerSize(1.0);  
  gr_pAu_2040_S_sys->SetMarkerSize(1.0);  
  gr_pAu_4060_S_sys->SetMarkerSize(1.0);  
  gr_pAu_6084_S_sys->SetMarkerSize(1.0);  

  gr_pAl_020_N_sys->SetMarkerSize(1.0); 
  gr_pAl_2040_N_sys->SetMarkerSize(1.0);  
  gr_pAl_4072_N_sys->SetMarkerSize(1.0);  
      
  gr_pAl_020_S_sys->SetMarkerSize(1.0);   
  gr_pAl_2040_S_sys->SetMarkerSize(1.0);  
  gr_pAl_4072_S_sys->SetMarkerSize(1.0);  
 
  gr_HeAu_020_N_sys->SetMarkerSize(1.0); 
  gr_HeAu_2040_N_sys->SetMarkerSize(1.0);  
  gr_HeAu_4088_N_sys ->SetMarkerSize(1.0); 

  gr_HeAu_020_S_sys->SetMarkerSize(1.0); 
  gr_HeAu_2040_S_sys->SetMarkerSize(1.0);  
  gr_HeAu_4088_S_sys->SetMarkerSize(1.0); 

  // Set Fill Style type B
  gr_pAu_05_N_sys->SetFillStyle(0);  
  gr_pAu_510_N_sys->SetFillStyle(0);  
  gr_pAu_1020_N_sys->SetFillStyle(0);  
  gr_pAu_020_N_sys->SetFillStyle(0);  
  gr_pAu_2040_N_sys->SetFillStyle(0);  
  gr_pAu_4060_N_sys->SetFillStyle(0);  
  gr_pAu_6084_N_sys->SetFillStyle(0);  

  gr_pAu_05_S_sys->SetFillStyle(0);  
  gr_pAu_510_S_sys->SetFillStyle(0);  
  gr_pAu_1020_S_sys->SetFillStyle(0);  
  gr_pAu_020_S_sys->SetFillStyle(0);  
  gr_pAu_2040_S_sys->SetFillStyle(0);  
  gr_pAu_4060_S_sys->SetFillStyle(0);  
  gr_pAu_6084_S_sys->SetFillStyle(0);  

  gr_pAl_020_N_sys->SetFillStyle(0); 
  gr_pAl_2040_N_sys->SetFillStyle(0);  
  gr_pAl_4072_N_sys->SetFillStyle(0);  
      
  gr_pAl_020_S_sys->SetFillStyle(0);   
  gr_pAl_2040_S_sys->SetFillStyle(0);  
  gr_pAl_4072_S_sys->SetFillStyle(0);  
 
  gr_HeAu_020_N_sys->SetFillStyle(0); 
  gr_HeAu_2040_N_sys->SetFillStyle(0);  
  gr_HeAu_4088_N_sys ->SetFillStyle(0); 

  gr_HeAu_020_S_sys->SetFillStyle(0); 
  gr_HeAu_2040_S_sys->SetFillStyle(0);  
  gr_HeAu_4088_S_sys->SetFillStyle(0); 
 /////// Set Line Width Type B
  gr_pAu_05_N_sys->SetLineWidth(2);  
  gr_pAu_510_N_sys->SetLineWidth(2);  
  gr_pAu_1020_N_sys->SetLineWidth(2);  
  gr_pAu_020_N_sys->SetLineWidth(2);  
  gr_pAu_2040_N_sys->SetLineWidth(2);  
  gr_pAu_4060_N_sys->SetLineWidth(2);  
  gr_pAu_6084_N_sys->SetLineWidth(2);  

  gr_pAu_05_S_sys->SetLineWidth(2);  
  gr_pAu_510_S_sys->SetLineWidth(2);  
  gr_pAu_1020_S_sys->SetLineWidth(2);  
  gr_pAu_020_S_sys->SetLineWidth(2);  
  gr_pAu_2040_S_sys->SetLineWidth(2);  
  gr_pAu_4060_S_sys->SetLineWidth(2);  
  gr_pAu_6084_S_sys->SetLineWidth(2);  

  gr_pAl_020_N_sys->SetLineWidth(2); 
  gr_pAl_2040_N_sys->SetLineWidth(2);  
  gr_pAl_4072_N_sys->SetLineWidth(2);  
      
  gr_pAl_020_S_sys->SetLineWidth(2);   
  gr_pAl_2040_S_sys->SetLineWidth(2);  
  gr_pAl_4072_S_sys->SetLineWidth(2);  
 
  gr_HeAu_020_N_sys->SetLineWidth(2); 
  gr_HeAu_2040_N_sys->SetLineWidth(2);  
  gr_HeAu_4088_N_sys ->SetLineWidth(2); 

  gr_HeAu_020_S_sys->SetLineWidth(2); 
  gr_HeAu_2040_S_sys->SetLineWidth(2);  
  gr_HeAu_4088_S_sys->SetLineWidth(2); 
  // set x axis limits
  gr_pAu_05_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_510_N->GetXaxis()->SetLimits(0,7);  
  gr_pAu_1020_N->GetXaxis()->SetLimits(0,7);  
  gr_pAu_020_N->GetXaxis()->SetLimits(0,7);  
  gr_pAu_2040_N->GetXaxis()->SetLimits(0,7);  
  gr_pAu_4060_N->GetXaxis()->SetLimits(0,7);  
  gr_pAu_6084_N->GetXaxis()->SetLimits(0,7);  

  gr_pAu_05_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_510_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_1020_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_020_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_2040_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_4060_S->GetXaxis()->SetLimits(0,7);  
  gr_pAu_6084_S->GetXaxis()->SetLimits(0,7);  

  gr_pAl_020_N->GetXaxis()->SetLimits(0,7); 
  gr_pAl_2040_N->GetXaxis()->SetLimits(0,7);  
  gr_pAl_4072_N->GetXaxis()->SetLimits(0,7);  
      
  gr_pAl_020_S->GetXaxis()->SetLimits(0,7);   
  gr_pAl_2040_S->GetXaxis()->SetLimits(0,7);  
  gr_pAl_4072_S->GetXaxis()->SetLimits(0,7);  
 
  gr_HeAu_020_N->GetXaxis()->SetLimits(0,7); 
  gr_HeAu_2040_N->GetXaxis()->SetLimits(0,7);  
  gr_HeAu_4088_N ->GetXaxis()->SetLimits(0,7); 

  gr_HeAu_020_S->GetXaxis()->SetLimits(0,7); 
  gr_HeAu_2040_S->GetXaxis()->SetLimits(0,7);  
  gr_HeAu_4088_S->GetXaxis()->SetLimits(0,7); 

  // set y aaxis limits
  gr_pAu_05_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_510_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_1020_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_020_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_2040_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_4060_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_6084_N->GetYaxis()->SetRangeUser(0,2);  

  gr_pAu_05_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_510_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_1020_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_020_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_2040_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_4060_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAu_6084_S->GetYaxis()->SetRangeUser(0,2);  

  gr_pAl_020_N->GetYaxis()->SetRangeUser(0,2); 
  gr_pAl_2040_N->GetYaxis()->SetRangeUser(0,2);  
  gr_pAl_4072_N->GetYaxis()->SetRangeUser(0,2);  
      
  gr_pAl_020_S->GetYaxis()->SetRangeUser(0,2);   
  gr_pAl_2040_S->GetYaxis()->SetRangeUser(0,2);  
  gr_pAl_4072_S->GetYaxis()->SetRangeUser(0,2);  
 
  gr_HeAu_020_N->GetYaxis()->SetRangeUser(0,2); 
  gr_HeAu_2040_N->GetYaxis()->SetRangeUser(0,2);  
  gr_HeAu_4088_N ->GetYaxis()->SetRangeUser(0,2); 

  gr_HeAu_020_S->GetYaxis()->SetRangeUser(0,2); 
  gr_HeAu_2040_S->GetYaxis()->SetRangeUser(0,2);  
  gr_HeAu_4088_S->GetYaxis()->SetRangeUser(0,2); 

  // set x axis title
  gr_pAu_05_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_510_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_1020_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_020_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_2040_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_4060_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_6084_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  

  gr_pAu_05_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_510_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_1020_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_020_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_2040_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_4060_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAu_6084_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  

  gr_pAl_020_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)"); 
  gr_pAl_2040_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAl_4072_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
      
  gr_pAl_020_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");   
  gr_pAl_2040_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_pAl_4072_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
 
  gr_HeAu_020_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)"); 
  gr_HeAu_2040_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_HeAu_4088_N ->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)"); 

  gr_HeAu_020_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)"); 
  gr_HeAu_2040_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");  
  gr_HeAu_4088_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)"); 

 
  
  ///////////////////////////////////////////////////////// RAB set marker color (comment out everything else)
  //  gr_pAu_05_N->SetMarkerColor(kBlue);
  // gr_pAu_510_N->SetMarkerColor(kBlack); 
  // gr_pAu_1020_N->SetMarkerColor(kViolet); 
  // gr_pAu_020_N->SetMarkerColor(kRed); 
  // gr_pAu_2040_N->SetMarkerColor(kCyan); 
  // gr_pAu_4060_N->SetMarkerColor(kGreen); 
  // gr_pAu_6084_N->SetMarkerColor(kMagenta); 

  // gr_pAu_05_S->SetMarkerColor(kMagenta); 
  //  gr_pAu_510_S->SetMarkerColor(kBlue); 
  //  gr_pAu_1020_S->SetMarkerColor(kBlack); 
  gr_pAu_020_S->SetMarkerColor(kRed); 
  // gr_pAu_2040_S->SetMarkerColor(kViolet); 
  // gr_pAu_4060_S->SetMarkerColor(kGreen); 
  // gr_pAu_6084_S->SetMarkerColor(kCyan); 

  // gr_pAl_020_N->SetMarkerColor(kBlue);
  // gr_pAl_2040_N->SetMarkerColor(kBlue); 
  // gr_pAl_4072_N->SetMarkerColor(kBlue); 
      
  // gr_pAl_020_S->SetMarkerColor(kBlue);  
  // gr_pAl_2040_S->SetMarkerColor(kBlue); 
  // gr_pAl_4072_S->SetMarkerColor(kBlue); 
 
  // gr_HeAu_020_N->SetMarkerColor(kBlue);
  // gr_HeAu_020_N->SetMarkerColor(kBlue); 
  // gr_HeAu_020_N->SetMarkerColor(kBlue); 
  // gr_HeAu_020_S->SetMarkerColor(kBlue);
  // gr_HeAu_020_S->SetMarkerColor(kBlue); 
  // gr_HeAu_020_S->SetMarkerColor(kBlue); 
   ///////////////////////////////////////////////////////////////  Type B set marker color

  //  gr_pAu_05_N_sys->SetMarkerColor(kBlue);  
  //  gr_pAu_510_N_sys->SetMarkerColor(kBlack);  
  // gr_pAu_1020_N_sys->SetMarkerColor(kViolet);  
  //  gr_pAu_020_N_sys->SetMarkerColor(kRed);  
  //  gr_pAu_2040_N_sys->SetMarkerColor(kCyan);  
  //  gr_pAu_4060_N_sys->SetMarkerColor(kGreen);  
  //  gr_pAu_6084_N_sys->SetMarkerColor(kMagenta);  

  // gr_pAu_05_S_sys->SetMarkerColor(kMagenta);  
  // gr_pAu_510_S_sys->SetMarkerColor(kBlue);  
  // gr_pAu_1020_S_sys->SetMarkerColor(kBlack);  
  gr_pAu_020_S_sys->SetMarkerColor(kRed);  
  // gr_pAu_2040_S_sys->SetMarkerColor(kViolet);  
  // gr_pAu_4060_S_sys->SetMarkerColor(kGreen);  
  // gr_pAu_6084_S_sys->SetMarkerColor(kCyan);  

  // gr_pAl_020_N_sys->SetMarkerColor(kBlue); 
  // gr_pAl_2040_N_sys->SetMarkerColor(kBlue);  
  // gr_pAl_4072_N_sys->SetMarkerColor(kBlue);  
      
  // gr_pAl_020_S_sys->SetMarkerColor(kBlue);   
  // gr_pAl_2040_S_sys->SetMarkerColor(kBlue);  
  // gr_pAl_4072_S_sys->SetMarkerColor(kBlue);  
 
  // gr_HeAu_020_N_sys->SetMarkerColor(kBlue); 
  // gr_HeAu_2040_N_sys->SetMarkerColor(kBlue);  
  // gr_HeAu_4088_N_sys ->SetMarkerColor(kBlue); 

  // gr_HeAu_020_S_sys->SetMarkerColor(kBlue); 
  // gr_HeAu_2040_S_sys->SetMarkerColor(kBlue);  
  // gr_HeAu_4088_S_sys->SetMarkerColor(kBlue)

   ///////////////////////////////////////////////// Set Line Color Type B box

  //  gr_pAu_05_N_sys->SetLineColor(kBlue);  
  // gr_pAu_510_N_sys->SetLineColor(kBlack);  
  // gr_pAu_1020_N_sys->SetLineColor(kViolet);  
  // gr_pAu_020_N_sys->SetLineColor(kRed);  
  // gr_pAu_2040_N_sys->SetLineColor(kCyan);  
  // gr_pAu_4060_N_sys->SetLineColor(kGreen);  
  // gr_pAu_6084_N_sys->SetLineColor(kMagenta);  

  // gr_pAu_05_S_sys->SetLineColor(kMagenta);  
  // gr_pAu_510_S_sys->SetLineColor(kBlue);  
  // gr_pAu_1020_S_sys->SetLineColor(kBlack);  
  gr_pAu_020_S_sys->SetLineColor(kRed);  
  // gr_pAu_2040_S_sys->SetLineColor(kViolet);  
  // gr_pAu_4060_S_sys->SetLineColor(kGreen);  
  // gr_pAu_6084_S_sys->SetLineColor(kCyan);  

  // gr_pAl_020_N_sys->SetLineColor(kBlue); 
  // gr_pAl_2040_N_sys->SetLineColor(kBlue);  
  // gr_pAl_4072_N_sys->SetLineColor(kBlue);  
      
  // gr_pAl_020_S_sys->SetLineColor(kBlue);   
  // gr_pAl_2040_S_sys->SetLineColor(kBlue);  
  // gr_pAl_4072_S_sys->SetLineColor(kBlue);  
 
  // gr_HeAu_020_N_sys->SetLineColor(kBlue); 
  // gr_HeAu_2040_N_sys->SetLineColor(kBlue);  
  // gr_HeAu_4088_N_sys ->SetLineColor(kBlue); 

  // gr_HeAu_020_S_sys->SetLineColor(kBlue); 
  // gr_HeAu_2040_S_sys->SetLineColor(kBlue);  
  // gr_HeAu_4088_S_sys->SetLineColor(kBlue)

  ///////////////////////////////////////////////// Set Fill Color Type B box

 // gr_pAu_05_N_sys->SetFillColor(kBlue);  
 //  gr_pAu_510_N_sys->SetFillColor(kBlack);  
 //  gr_pAu_1020_N_sys->SetFillColor(kViolet);  
 //  gr_pAu_020_N_sys->SetFillColor(kRed);  
 //  gr_pAu_2040_N_sys->SetFillColor(kCyan);  
 //  gr_pAu_4060_N_sys->SetFillColor(kGreen);  
 //  gr_pAu_6084_N_sys->SetFillColor(kMagenta);  

 // gr_pAu_05_S_sys->SetFillColor(kMagenta);  
 // gr_pAu_510_S_sys->SetFillColor(kBlue);  
 // gr_pAu_1020_S_sys->SetFillColor(kBlack);  
  gr_pAu_020_S_sys->SetFillColor(kRed);  
  // gr_pAu_2040_S_sys->SetFillColor(kViolet);  
  // gr_pAu_4060_S_sys->SetFillColor(kGreen);  
  // gr_pAu_6084_S_sys->SetFillColor(kCyan);  

  // gr_pAl_020_N_sys->SetFillColor(kBlue); 
  // gr_pAl_2040_N_sys->SetFillColor(kBlue);  
  // gr_pAl_4072_N_sys->SetFillColor(kBlue);  
      
  // gr_pAl_020_S_sys->SetFillColor(kBlue);   
  // gr_pAl_2040_S_sys->SetFillColor(kBlue);  
  // gr_pAl_4072_S_sys->SetFillColor(kBlue);  
 
  // gr_HeAu_020_N_sys->SetFillColor(kBlue); 
  // gr_HeAu_2040_N_sys->SetFillColor(kBlue);  
  // gr_HeAu_4088_N_sys ->SetFillColor(kBlue); 

  // gr_HeAu_020_S_sys->SetFillColor(kBlue); 
  // gr_HeAu_2040_S_sys->SetFillColor(kBlue);  
  // gr_HeAu_4088_S_sys->SetFillColor(kBlue)

  ////////////////////////////////////////////////////// Select arrays to plot (comment out everything else)
  // gr_pAu_05_N->Draw("AP");
  // gr_pAu_510_N->Draw("P"); 
  // gr_pAu_1020_N->Draw("P"); 
  // gr_pAu_020_N->Draw("P"); 
  // gr_pAu_2040_N->Draw("P"); 
  // gr_pAu_4060_N->Draw("P"); 
  // gr_pAu_6084_N->Draw("P"); 

  // gr_pAu_05_S->Draw("AP"); 
  // gr_pAu_510_S->Draw("P"); 
  // gr_pAu_1020_S->Draw("P"); 
  gr_pAu_020_S->Draw("AP"); 
  // gr_pAu_2040_S->Draw("P"); 
  //  gr_pAu_4060_S->Draw("P"); 
  //  gr_pAu_6084_S->Draw("P"); 

  // gr_pAl_020_N->Draw("AP");
  // gr_pAl_2040_N->Draw("AP"); 
  // gr_pAl_4072_N->Draw("AP"); 
      
  // gr_pAl_020_S->Draw("AP");  
  // gr_pAl_2040_S->Draw("AP"); 
  // gr_pAl_4072_S->Draw("AP"); 
 
  // gr_HeAu_020_N->Draw("AP");
  // gr_HeAu_2040_N->Draw("AP"); 
  // gr_HeAu_4088_N->Draw("AP");
  // gr_HeAu_020_S->Draw("AP");
  // gr_HeAu_2040_S->Draw("AP"); 
  // gr_HeAu_4088_S->Draw("AP"); 
      
  ////////////////////////////////////////////////////////////////////////////////// Add systematic uncertainty boxes

  // gr_pAu_05_N_sys->Draw("e2same");  
  // gr_pAu_510_N_sys->Draw("e2same");  
  // gr_pAu_1020_N_sys->Draw("e2same");  
  // gr_pAu_020_N_sys->Draw("e2same");  
  // gr_pAu_2040_N_sys->Draw("e2same");  
  // gr_pAu_4060_N_sys->Draw("e2same");  
  // gr_pAu_6084_N_sys->Draw("e2same");  

  // gr_pAu_05_S_sys->Draw("e2same");  
  // gr_pAu_510_S_sys->Draw("e2same");  
  // gr_pAu_1020_S_sys->Draw("e2same");  
  gr_pAu_020_S_sys->Draw("e2same");  
  //  gr_pAu_2040_S_sys->Draw("e2same");  
  //  gr_pAu_4060_S_sys->Draw("e2same");  
  //  gr_pAu_6084_S_sys->Draw("e2same");  

  // gr_pAl_020_N_sys->Draw("e2same"); 
  // gr_pAl_2040_N_sys->Draw("e2same");  
  // gr_pAl_4072_N_sys->Draw("e2same");  
      
  // gr_pAl_020_S_sys->Draw("e2same");   
  // gr_pAl_2040_S_sys->Draw("e2same");  
  // gr_pAl_4072_S_sys->Draw("e2same");  
 
  // gr_HeAu_020_N_sys->Draw("e2same"); 
  // gr_HeAu_2040_N_sys->Draw("e2same");  
  // gr_HeAu_4088_N_sys ->Draw("e2same"); 

  // gr_HeAu_020_S_sys->Draw("e2same"); 
  // gr_HeAu_2040_S_sys->Draw("e2same");  
  // gr_HeAu_4088_S_sys->Draw("e2same")

     
  // // draw boxes

  // for(int i = 0; i < 19; i++)
  // 	{
  // 	  draw_box_pAu_N_05(pt_array_pAu[i],  RpAu_05_array[1][i], box_pt,  RpAu_05_sys_frac_array[1][i]* RpAu_05_array[1][i], box_color); 
  // 	}

  ////////////////////////////////////////////////////////////////////////////////// Add legend entries (comment out everything else)

  TLegend *leg = new TLegend(0.1, 0.8, 0.35, 0.9);  //(start x, start y, end x, end y)
  leg->SetFillColor(0); 
  leg->SetTextSize(0.03); 

  //  leg->AddEntry(gr_pAu_05_N, "Run15pAu, N 0-5 ", "p"); 
  // leg->AddEntry(gr_pAu_510_N, "Run15pAu, N 5-10 ", "p"); 
  // leg->AddEntry(gr_pAu_1020_N, "Run15pAu, N 10-20 ", "p"); 
  // leg->AddEntry(gr_pAu_020_N, "Run15pAu, N 0-20 ", "p"); 
  // leg->AddEntry(gr_pAu_2040_N, "Run15pAu, N 20-40 ", "p"); 
  // leg->AddEntry(gr_pAu_4060_N, "Run15pAu, N 40-60 ", "p"); 
  // leg->AddEntry(gr_pAu_6084_N, "Run15pAu, N 60-84 ", "p"); 

  // leg->AddEntry(gr_pAu_05_S, "Run15pAu, S 0-5 ", "p"); 
  // leg->AddEntry(gr_pAu_510_S, "Run15pAu, S 5-10 ", "p"); 
  // leg->AddEntry(gr_pAu_1020_S, "Run15pAu, S 10-20 ", "p"); 
  leg->AddEntry(gr_pAu_020_S, "Run15pAu, S 0-20 ", "p"); 
  // leg->AddEntry(gr_pAu_2040_S, "Run15pAu, S 20-40 ", "p"); 
  // leg->AddEntry(gr_pAu_4060_S, "Run15pAu, S 40-60 ", "p"); 
  // leg->AddEntry(gr_pAu_6084_S, "Run15pAu, S 60-84 ", "p"); 
  
  // leg->AddEntry(gr_pAl_020_N, "Run15pAl, N 0-20 ", "p"); 
  // leg->AddEntry(gr_pAl_2040_N, "Run15pAl, N 20-40 ", "p"); 
  // leg->AddEntry(gr_pAl_4072_N, "Run15pAl, N 40-72 ", "p"); 
  // leg->AddEntry(gr_pAl_020_S, "Run15pAl, S 0-20 ", "p"); 
  // leg->AddEntry(gr_pAl_2040_S, "Run15pAl, S 20-40 ", "p"); 
  // leg->AddEntry(gr_pAl_4072_S, "Run15pAl, S 40-72 ", "p"); 

  // leg->AddEntry(gr_HeAu_020_N, "Run15HeAu, N 0-20 ", "p"); 
  // leg->AddEntry(gr_HeAu_2040_N, "Run15HeAu, N 20-40 ", "p"); 
  // leg->AddEntry(gr_HeAu_4088_N, "Run15HeAu, N 40-88", "p"); 
  // leg->AddEntry(gr_HeAu_020_S, "Run15HeAu, S 0-20 ", "p"); 
  // leg->AddEntry(gr_HeAu_2040_S, "Run15HeAu, S 20-40 ", "p"); 
  // leg->AddEntry(gr_HeAu_4088_S, "Run15HeAu, S 40-88 ", "p"); 

  leg->Draw();  
 
  
  

}




 // gr_pAu_05_N
 //  gr_pAu_510_N 
 //  gr_pAu_1020_N 
 //  gr_pAu_020_N 
 //  gr_pAu_2040_N 
 //  gr_pAu_4060_N 
 //  gr_pAu_6084_N 

 //  gr_pAu_05_S 
 //  gr_pAu_510_S 
 //  gr_pAu_1020_S 
 //  gr_pAu_020_S 
 //  gr_pAu_2040_S 
 //  gr_pAu_4060_S 
 //  gr_pAu_6084_S 

 //  gr_pAl_020_N
 //  gr_pAl_2040_N 
 //  gr_pAl_4072_N 
      
 //  gr_pAl_020_S  
 //  gr_pAl_2040_S 
 //  gr_pAl_4072_S 
 
 //  gr_HeAu_020_N
 //  gr_HeAu_020_N 
 //  gr_HeAu_020_N 
 //  gr_HeAu_020_S
 //  gr_HeAu_020_S 
 //  gr_HeAu_020_S 
