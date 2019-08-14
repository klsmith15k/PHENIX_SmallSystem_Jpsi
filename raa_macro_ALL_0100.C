
///////////////////// For p+Au,p+Al,He+Au /////////////////////////
// have not added correct systematic uncertainties for centrality plots

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
#include <TGraphAsymmErrors.h>

using namespace std;
                                       
void draw_box_N( double x, double y, double delta_x, double delta_y, int color, int flag = 0, int style = 1 )
{
  double x0 = x+delta_x;
  double x1 = x-delta_x;
  double y0 = y+delta_y;
  double y1 = y-delta_y;
  TBox *box1 = new TBox( x0, y0, x1, y1 );
  if( style != 0 ) 
    box1->SetFillStyle( style );
  box1->SetLineColor(kBlue);                          
  box1->SetFillColor(kBlue);                          
  box1->SetLineWidth(2);                           
  box1->SetFillStyle(0);
  box1->Draw();
}

void draw_box_S( double x, double y, double delta_x, double delta_y, int color, int flag = 0, int style = 1 )
{
  double x0 = x+delta_x;
  double x1 = x-delta_x;
  double y0 = y+delta_y;
  double y1 = y-delta_y;
  TBox *box2 = new TBox( x0, y0, x1, y1 );
  if( style != 0 ) 
    box2->SetFillStyle( style );
  box2->SetLineColor(kViolet);                          
  box2->SetFillColor(kViolet);                          
  box2->SetLineWidth(2);                           
  box2->SetFillStyle(0);
  box2->Draw();
}
void draw_box_dAu( double x, double y, double delta_x, double delta_y, int color, int flag = 0, int style = 1 )
{
  double x0 = x+delta_x;
  double x1 = x-delta_x;
  double y0 = y+delta_y;
  double y1 = y-delta_y;
  TBox *box3 = new TBox( x0, y0, x1, y1 );
  if( style != 0 ) 
    box3->SetFillStyle( style );
  box3->SetLineColor(kRed);                          
  box3->SetFillColor(kRed);                          
  box3->SetLineWidth(2);                           
  box3->SetFillStyle(0);
  box3->Draw();
}
void draw_box_black( double x, double y, double delta_x, double delta_y, int color, int flag = 0, int style = 1 )
{
  double x0 = x+delta_x;
  double x1 = x-delta_x;
  double y0 = y+delta_y;
  double y1 = y-delta_y;
  TBox *box3 = new TBox( x0, y0, x1, y1 );
  if( style != 0 ) 
    box3->SetFillStyle( style );
  box3->SetLineColor(kBlack);                          
  box3->SetFillColor(kBlack);                          
  box3->SetLineWidth(2);                           
  box3->SetFillStyle(0);
  box3->Draw();
}


#include "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/RpAu_RdAu_pt_arrays.C"
#include "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run14HeAu_efficiencies_arrays.C"
//include "centrality_TypeB_arrays.C"

void raa_macro_ALL_0100()
{


  int i_system;
  cout << "Enter the system ('0' for Run15pAu, '1' for Run15pAl, '2' for Run14HeAu or '3' for Run08dAu)" << endl;
  cin >> i_system;
  int i_cent;
  cout << "Enter the centrality range ('0' for no centrality" << endl;
  cin >> i_cent;


  ///////////////////////////////////////////////// dAu information
 
  int bins = 26;

  double RAA_dAu_frac_sys[2][5][28];
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0;i < bins-3; i++)
	{
	  RAA_dAu_frac_sys[arm][i_cent][i] =  RAA_dAu_ERR_2[arm][i_cent][i]/RAA_dAu[arm][i_cent][i];
	}
    }


  double rA[2] = {1.0,1.0}; // initialize to 1 in case not using
  double r_pp[2] = {1.0,1.0};

  bool write = false;
  bool binshift_A = true; 
  bool binshift_pp = true; 
  bool print_screen = true;

  // [i_cent][i_system]
  double N_coll[5][3] = {4.667,2.10,10.4, // pAu,pAl,HeAu no centrality 
			 8.2,3.35,22.3,
			 6.1,2.3,14.8,
			 4.4,1.7,5.48, 
			 2.6,0,0}; 
  
  double bias_correction[5][3] = {0.858,0.80,0.89, // no centrality pAu,pAl,HeAu 
				  0.90,0.81,0.95, // 020
				  0.98,0.90,1.01, // 2040 
				  1.02,1.01,1.03,  // 4060/40-72/40-88 
				  1,0,0};  // 6084  
  
  int data_points[5][3] = {26,20,19, // no centrality (pAu,pAl,HeAu)
			   19,16,13, // 020
			   19,16,13, 
			   19,16,13, 
			   19,0,0}; // 6084

  // [i_cent,i_system,data_points]
  double pt_array[5][3][26] = {
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.25,6.75,     
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,5.5,6.5,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // no centrality 

    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},   // 0-20  (pAu, pAl and HeAu share same binning)

    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},   // 20-40

    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, //40-60/40-72/40-88

    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; //60-84

  double pt_width[5][3][28] = {  
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.0,0.0, // Run15pAu
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // Run15pAl 
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // Run14HeAu    NO CENTRALITY

    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  // cent 0-20

    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // cent 20-40

    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // cent 40-60/40-72/40-88

    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; // cent 60-84
  
  double invariant_yield = 0;
  double invariant_yield_err = 0;
  double two_pi = 2*3.141592653;
  double muid_trig[2][28] = {0};
  double acc_reco[2][28] = {0};
  double muid_trig_eff = 0;
  double acc_reco_eff = 0;
  double d_pt = 0;
  double d_rap = 1.0; 
  double denom = 0;
  double jpsi_counts[2] = {0};
  double jpsi_err_1[2] = {0};
  double jpsi_err_2[2] = {0};
  double exp = pow(10,-9);

    
  // [i_cent][arm][i_system]
  double MB[5][2][3] = {2.09769926640*pow(10,11), 2.04722*pow(10,11), 2.087*pow(10,10), // no centrality (pAu,pAl,HeAu)
			 1.84122015616*pow(10,11), 2.00743*pow(10,11),3.32923*pow(10,10),  // 

			 2.10268903808*pow(10,11)*0.2431, 2.04722*pow(10,11)*(0.2900), 2.087*pow(10,10)*(0.2/0.88), // cent 0-20
			 1.9356537088*pow(10,11)*0.2426, 2.00743*pow(10,11)*(0.2899), 3.32923*pow(10,10)*(0.2/0.88),  
			 
			 2.10268903808*pow(10,11)*0.2381, 2.04722*pow(10,11)*(0.277), 2.087*pow(10,10)*(0.2/0.88), // cent 20-40
			 1.9356537088*pow(10,11)*0.2382, 2.00743*pow(10,11)*(0.277), 3.32923*pow(10,10)*(0.2/0.88),  
			 
			 2.10268903808*pow(10,11)*0.2361, 2.04722*pow(10,11)*(0.4327), 2.087*pow(10,10)*(0.48/0.88),  // cent 40-60/40-72/40-88
			 1.9356537088*pow(10,11)*0.2363, 2.00743*pow(10,11)*(0.433), 3.32923*pow(10,10)*(0.48/0.88),  
			 
			 2.10268903808*pow(10,11)*0.2826,1,1,  // cent 60-84  no 4th centrality bin for pAl or HeAu
			 1.9356537088*pow(10,11)*0.2829,1,1};  

 


  double r_N = 0;
  double r_S = 0;
  double rA_N = 0;
  double rA_S = 0;

  double p0,p1;
  double f0,f1;
  
  TF1 *a[5][2][3]; //[i_cent][arm][system]
  TF1 *t[5][2][3];
 
  double pt_err[28] = {0};
  double run15[2][28] = {0};
  double run15_err[2][28] = {0};
  double pp_inv_corrected[2][28] = {0};
  
  double RAA_array[5][2][28] = {0};
  double RAA_ERR_1[5][2][28] = {0};

  double pp_inv_yield[2][28] = {0};
  double pp_inv_yield_err[2][28] = {0};
    
  double tmp_N = 0;
  double tmp_S = 0;
  double err_N = 0;
  double pt = 0;
  
  TFile *file_data_acc;
  TFile *file_data_trig;
  TFile *file_theory_curves = new TFile();

  TGraphAsymmErrors *epp1;
  TGraphAsymmErrors *epp0;
  TGraphAsymmErrors *ct1;
  TGraphAsymmErrors *ct0;
      
  file_theory_curves = TFile::Open("/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/RpA_Jpsi_RpA_EPPS16_nCTEQ15.root");
  
  file_theory_curves->GetObject("gepps16_pAu_pT_fwd",epp1);
  file_theory_curves->GetObject("gepps16_pAu_pT_bwd",epp0);
  file_theory_curves->GetObject("gncteq15_pAu_pT_fwd",ct1);
  file_theory_curves->GetObject("gncteq15_pAu_pT_bwd",ct0);
  
    // includes all centrality filenames --- . don't need to do query replace
  std::string acc_filename[5][3] = {
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0084_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT0072_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT0088_acceff_pT_InclusiveJpsi_combined.root"},
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0020_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT0020_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT0020_acceff_pT_InclusiveJpsi_combined.root"},
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT2040_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT2040_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT2040_acceff_pT_InclusiveJpsi_combined.root"},
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT4060_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT4072_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT4060_acceff_pT_InclusiveJpsi_combined.root"},
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT6084_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT4072_acceff_pT_InclusiveJpsi_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT6088_acceff_pT_InclusiveJpsi_combined.root"} };
   
  std::string obj_acc_filename[2][3] = {"frecoeff_pT_pAu_arm0","frecoeff_pT_pAl_arm0","frecoeff_pT_HeAu_arm0",
					"frecoeff_pT_pAu_arm1","frecoeff_pT_pAl_arm1","frecoeff_pT_HeAu_arm1"};
   
  std::string system_name[3] = {"p + Au","p + Al","He + Au"};

  //fill contents in a[arm]
  for(int cent = 0; cent < 5; cent++)
    {
      for(int arm = 0; arm < 2; arm++)
	{ 
	  for(int j_system = 0; j_system < 3; j_system++)
	    {
	      file_data_acc = TFile::Open(acc_filename[cent][j_system].c_str()); 
	      file_data_acc->GetObject(obj_acc_filename[arm][j_system].c_str(),a[cent][arm][j_system]);
	    }
	}
    }

  std::string trig_filename[5][3] = {
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0084_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT0072_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT0088_trigeff_MUID2D_pT_output_combined.root"}, //0-100     // centrality set files
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0020_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT0020_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT0020_trigeff_MUID2D_pT_output_combined.root"}, //0-20
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT2040_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT2040_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT2040_trigeff_MUID2D_pT_output_combined.root"}, //20-40
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT4060_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT4072_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT4060_trigeff_MUID2D_pT_output_combined.root"}, //40-60,40-72
    {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT6084_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAl200_CENT4072_trigeff_MUID2D_pT_output_EMBED_C.root", "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run14HeAu200_CENT6088_trigeff_MUID2D_pT_output_combined.root"} }; //60-84
    
  
  std::string trig_obj_filename[2][3] = {"ftrigeff_pT_dAu_arm0","ftrigeff_pT_pAl_arm0","ftrigeff_pT_HeAu_arm0",
  					 "ftrigeff_pT_dAu_arm1","ftrigeff_pT_pAl_arm1","ftrigeff_pT_HeAu_arm1"};
  
  //fill contents in t[arm]
  for(int cent = 0; cent < 5; cent++)
    {
      for(int arm = 0; arm < 2; arm++)
	{ 
	  for(int j_system = 0; j_system < 3; j_system++)
	    {
	      file_data_trig = TFile::Open(trig_filename[cent][j_system].c_str());
	      file_data_trig->GetObject(trig_obj_filename[arm][j_system].c_str(),t[cent][arm][j_system]);
	      
	    }
	}
    }

  ///////////////////////////////////////////////////  
  
  //number of groups of data(centrality),rows per group(arm),elements per row (system type)
  // read in pA jpsi counts
  char basename[30][800] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run15pAu_S_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run15pAl_S_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run14HeAu_S_bestfit_parameters_tony_1_", // 0-100 *****yuehang******* centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run15pAu_N_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run15pAl_N_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/Run14HeAu_N_bestfit_parameters_tony_1_", // 0-100 *****yuehang******** centrality, arm 1

			     "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 0-20 centrality, arm 0

			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 0-20 centrality, arm 1

			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 20-40 centrality, arm 0


			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 20-40 centrality, arm 1

			      "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4060_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/TH2D/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 40-60 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4060_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/TH2D/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 40-60 centrality, arm 1

			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6072_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 60-84 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_"}; // 60-84 centrality, arm 1
   
  // pA binshift correction files 
  char basename2[30][800] = {
"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 0-100 *****yuehang****** centrality, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check//HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 0-100 *******yuehang******* centrality, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 020, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 020, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 2040, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/2040_Centrality/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 2040, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4060_Centrality/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 4060, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4060_Centrality/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 4060, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 6084, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/6084_Centrality/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4072_Centrality/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/4088_Centrality/HeAu_binshift_corrections/Run14HeAu_N_correction_r_"}; // 6084, arm 1 
	
  // for pp inv. cross section [arm][system]
  char basename3[30][800] = {"/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_cross_sections/Run15pp_cross_S_", // arm 0 ******yuehang********0-100
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_cross_sections/Run15pp_cross_N_", // arm 1******yuehang******** 0-100

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_S_", // arm 0 0-20
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_N_", //arm 0 0-20

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_S_", // arm 0 20-40
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_N_", // arm 1 20-40

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_S_", // arm 0 40-60/40-72/40-88
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_N_", // arm 1 40-60/40-72/40-88

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_S_", //arm 0 60-84
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run15pp_cross_N_"}; // arm 1 60-84

  // for writing the  AB cross sections [cent][arm][system] to then get binshift corrcetions from binshift_TGraph
  char basename4[30][800] = {"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/HeAu_cross_sections/Run14HeAu_cross_S_", // 0-100 ****yuehang*****arm 0
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/cent_int_prelim_check/HeAu_cross_sections/Run14HeAu_cross_N_", // 0-100 ****yuehang***** arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_S_", // 0-20 arm 0
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_N_", // 0-20 arm 1

"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_S_", // 20-40 arm 0
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_N_", // 20-40 arm 1

"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_S_", // 40-60/40-72/40-88 arm 0
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_N_", // 40-60/40-72/40-88 arm 1

"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_S_", // 60-84 arm 0
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_cross_sections/Run14HeAu_cross_N_"};  // 60-84 arm 1

   
  // pp binshift correction files 
  char basename5[30][800] = {
"/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 0-100 centrality, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 0-100 centrality, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 020, arm 0

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 020, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 2040, arm 0

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 2040, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 4060, arm 0

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 4060, arm 1

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 6084, arm 0

			     "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/HeAu_binshift_corrections/Run14HeAu_N_correction_r_"}; // 6084, arm 1 
	

  vector < vector < vector < vector <std::string> > > > pp_cross_filename;  // for reading in (basename 3)
  vector < vector < vector < vector <std::string> > > > AA_cross_filename;  // for writing out (basename 4)
     
  for(int j_cent = 0; j_cent < 5; j_cent++)
    {
      vector < vector < vector <std::string> > > filename_type;
      vector < vector < vector <std::string> > > filename_type2;
      
      for(int j_arm = 0; j_arm < 2; j_arm++)
	{
	  vector < vector< std::string > > filename_sys;
	  vector < vector< std::string > > filename_sys2;
	  
	  for(int j_system = 0; j_system < 3; j_system++)
	    {
	      int j = j_cent*6 + j_arm*3 + j_system; // how to access the elements for the 1D array "basename".  The first 6 are bestfit files for no centrality (i_cent 0).  The first 3 elements are the South arm, the next 3 are the North arm, then 3 are South, 3 are North.  Then the system changes for each element in the row arm 0. 
	      //cout << "j: " << j << endl;
	      vector <std::string> filename_pt; 
	      vector <std::string> filename_pt2; 
	      
	      for(int i = 0; i < data_points[j_cent][j_system]; i++)
		{
		  char name[800];
		  sprintf(name,"%s%i%s",basename3[j],i+1,".dat");	     
		  std::string filename_tmp = name;
		  filename_pt.push_back(filename_tmp);
		  //cout << "i+1 = " << i+1 << endl;
		  
		  char name2[800];
		  sprintf(name2,"%s%i%s",basename4[j],i+1,".dat");	     
		  std::string filename_tmp2 = name2;
		  filename_pt2.push_back(filename_tmp2);
		} 
	      filename_sys.push_back(filename_pt);
	      filename_sys2.push_back(filename_pt2);
	    }
	  filename_type.push_back(filename_sys);
	  filename_type2.push_back(filename_sys2);
	}
      pp_cross_filename.push_back(filename_type);
      AA_cross_filename.push_back(filename_type2);
    } // 
  
  vector < vector < vector < vector <std::string> > > > bestfit_filename;  // [cent][arm][system][pt]
  vector < vector < vector < vector <std::string> > > > binshift_A_filename;
  vector < vector < vector < vector <std::string> > > > binshift_pp_filename;
  
  for(int j_cent = 0; j_cent < 5; j_cent++)
    {
      vector < vector < vector <std::string> > > filename_type;
      vector < vector < vector <std::string> > > filename_type2;
      vector < vector < vector <std::string> > > filename_type3;
      
      for(int j_arm = 0; j_arm < 2; j_arm++)
	{
	  vector < vector< std::string > > filename_sys;
	  vector < vector< std::string > > filename_sys2;
	  vector < vector< std::string > > filename_sys3;
	  
	  for(int j_system = 0; j_system < 3; j_system++)
	    {
	      int j = j_cent*6 + j_arm*3 + j_system;
	      vector <std::string> filename_pt; 
	      vector <std::string> filename_pt2; 
	      vector <std::string> filename_pt3;
	      
	      for(int i = 0; i < data_points[j_cent][j_system]; i++)
		{
		  char name[800];
		  sprintf(name,"%s%i%s",basename[j],i+1,".dat");	     
		  std::string filename_tmp = name;
		  filename_pt.push_back(filename_tmp);
		  
		  char name2[800];
		  sprintf(name2,"%s%i%s",basename2[j],i+1,".dat");	     
		  std::string filename_tmp2 = name2;
		  filename_pt2.push_back(filename_tmp2);

		  char name3[800];
		  sprintf(name3,"%s%i%s",basename5[j],i+1,".dat");	     
		  std::string filename_tmp3 = name3;
		  filename_pt3.push_back(filename_tmp3);
		} 
	      filename_sys.push_back(filename_pt);
	      filename_sys2.push_back(filename_pt2);
	      filename_sys3.push_back(filename_pt3);
	    }
	  filename_type.push_back(filename_sys);
	  filename_type2.push_back(filename_sys2);
	  filename_type3.push_back(filename_sys3);
	}
      bestfit_filename.push_back(filename_type);
      binshift_A_filename.push_back(filename_type2);
      binshift_pp_filename.push_back(filename_type3);
    } // 
  
  for(int l = 0; l < 2; l++) // arm
    {
      for(int k = 0; k < 5; k++) // centrality
        {
  	 for(int i = 0; i < data_points[k][i_system]; i++)
  	   {
  	     //cout << bestfit_filename[k][l][i_system][i].c_str() << endl;
	     cout << pp_cross_filename[k][l][i_system][i].c_str() << endl;
  	   }
        }
    } 
 

  // print out efficiencies

  for(int cent = 0; cent < 5; cent++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int system = 0; system < 3; system++)
	    {
	      for(int i = 0; i < data_points[i_cent][i_system];i++)
		{
		  pt = pt_array[cent][system][i];
		  acc_reco_eff = a[cent][arm][system]->Eval(pt);
		  muid_trig_eff = t[cent][arm][system]->Eval(pt);
		  if(cent == i_cent  && system == i_system)
		    {
		      //cout << muid_trig_eff << endl;
		      // cout << acc_reco_eff << endl;
		      muid_trig[arm][i] = muid_trig_eff;
		      acc_reco[arm][i] = acc_reco_eff;
		    }
		}
	    }
	}
    }

  // BEGIN RAA Calc

  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < data_points[i_cent][i_system];i++)
	{
	  pt = pt_array[i_cent][i_system][i];
	  d_pt = pt_width[i_cent][i_system][i];
	  pt_err[i] = 0.0;
	  
	  acc_reco_eff = a[i_cent][arm][i_system]->Eval(pt);
	  muid_trig_eff = t[i_cent][arm][i_system]->Eval(pt);
	
	  if(i_system == 2 && i_cent == 3)  // read in the weighted efficiencies for HeAu 4088 from .C file
	    {
	      acc_reco_eff = HeAu_acc_eff_4088[arm][i];
	      muid_trig_eff = HeAu_trig_eff_4088[arm][i];
	    }

	  ifstream bestfit_parameters(bestfit_filename[i_cent][arm][i_system][i].c_str());
	  if(bestfit_parameters)
	    {
	      do 
		{
		  double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, NJpsi, err_NJpsi;	  
		  bestfit_parameters >> f0 >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> NJpsi >> err_NJpsi; 
		  jpsi_counts[arm] = NJpsi;
		  jpsi_err_1[arm] = err_NJpsi;
		  
		} while(bestfit_parameters.good());
	    }
	  else
	    cout << "Did not find fit files for i_system " << i_system << " and i_centrality :" << i_cent << endl;;
	  
	  if(arm == 0)
	    cout << round(jpsi_counts[0]) << endl;

	  denom = two_pi*pt*d_pt*d_rap*acc_reco_eff*muid_trig_eff*MB[i_cent][arm][i_system];
	  invariant_yield = (jpsi_counts[arm]*bias_correction[i_cent][i_system])/denom; 
	  invariant_yield_err = (jpsi_err_1[arm]*bias_correction[i_cent][i_system])/denom;
	  
	  if(binshift_A == true)
	    {
	      ifstream binshift(binshift_A_filename[i_cent][arm][i_system][i].c_str());
	      if(binshift)
		{
		  do 
		    {
		      double binshift_correction;	  
		      binshift >> binshift_correction; 
		      rA[arm] = binshift_correction;
		      
		    } while(binshift.good());
		}
	      else
		cout << "Did not find pA binshift files for system " << i_system << " and centrality :" << i_cent << endl;;
	    }

	  if(binshift_pp == true)
	    {
	      ifstream binshift(binshift_pp_filename[i_cent][arm][i_system][i].c_str());
	      if(binshift)
		{
		  do 
		    {
		      double binshift_correction;	  
		      binshift >> binshift_correction; 
		      r_pp[arm] = binshift_correction;
		      
		    } while(binshift.good());
		}
	      else
		cout << "Did not find pp binshift files for system " << i_system << " and centrality :" << i_cent << endl;;
	    }
	  
	  run15[arm][i] = invariant_yield/rA[arm];
	  run15_err[arm][i] = invariant_yield_err/rA[arm];;
	  
	  std::ifstream Run15pp_cross(pp_cross_filename[i_cent][arm][i_system][i].c_str(),std::fstream::in); 
	  Run15pp_cross >>  p0  >>  p1; 
	  pp_inv_yield[arm][i] = p0/(42*0.001);
	  pp_inv_yield_err[arm][i] = p1/(42*0.001); 
	  Run15pp_cross.close();
	  
	  pp_inv_yield[arm][i] /= r_pp[arm];
	  pp_inv_yield_err[arm][i] /= r_pp[arm];

	  RAA_array[i_cent][arm][i] = 1/N_coll[i_cent][i_system]*(run15[arm][i]/pp_inv_yield[arm][i]); 
	  double tmp2 =(pow(1/N_coll[i_cent][i_system]*run15_err[arm][i]/run15[arm][i],2)+pow(1/N_coll[i_cent][i_system]*pp_inv_yield_err[arm][i]/pp_inv_yield[arm][i],2))*(pow(run15[arm][i]/pp_inv_yield[arm][i],2));
	  RAA_ERR_1[i_cent][arm][i] = sqrt(tmp2);
	  
	  if(write == true)
	    {
	      f0 = run15[arm][i]*42*0.001;
	      f1 = run15_err[arm][i]*42*0.001; 
	      
	      std::fstream AA_cross(AA_cross_filename[i_cent][arm][i_system][i].c_str(),std::ofstream::out); 
	      AA_cross <<  f0  <<  " " << f1; 
	      AA_cross.close();
	    }

	}//i for loop
    } // end arm loop
   
  /////////////////////////////////////////////////// print out arrays

  int bins_system = data_points[i_cent][i_system];
  
  if(print_screen == true)
    {
      cout << "_______________________________ " << endl;
      cout << "RAA_N_ARRAY[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++)
	{
	  cout << RAA_array[i_cent][1][i] << ", ";
	}
      cout <<  RAA_array[i_cent][1][bins_system-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "RAA_N_ERR_1[" <<bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++)
      	{
	  cout << RAA_ERR_1[i_cent][1][i] << ", ";
      	}
      cout <<  RAA_ERR_1[i_cent][1][bins_system-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "RAA_S_ARRAY[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++)
    	{
    	  cout << RAA_array[i_cent][0][i] << ", ";
    	}
      cout <<  RAA_array[i_cent][0][bins_system-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "RAA_S_ERR_1[" << bins_system  << "] = {";
      for(int i = 0;i < bins_system -1; i++)
      	{
      	  cout << RAA_ERR_1[i_cent][0][i] << ", ";
      	}
      cout <<  RAA_ERR_1[i_cent][0][bins_system -1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "pt_array[" << bins_system  << "] = {";
      for(int i = 0;i < bins_system -1; i++)
      	{
      	  cout << pt_array[i_cent][i_system][i] << ", ";
      	}
      cout <<  pt_array[i_cent][i_system][bins_system -1] << "};" << endl;
      /////////////



      /////////////// inv yields printed out:
  cout << "INV_YIELD_N[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // north AA inv yield
    	{
    	  cout << run15[1][i] << ", ";
    	}
      cout <<  run15[1][bins_system-1] << "};" << endl;
  cout << "INV_YIELD_N_ERR[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // north AA inv yield err
    	{
    	  cout << run15_err[1][i] << ", ";
    	}
      cout <<  run15_err[1][bins_system-1] << "};" << endl;
  cout << "INV_YIELD_S[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield
	{
	  cout << run15[0][i] << ", ";
	}
      cout <<  run15[0][bins_system-1] << "};" << endl;
  cout << "INV_YIELD_S_ERR[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield err
    	{
    	  cout << run15_err[0][i] << ", ";
    	}
      cout <<  run15_err[0][bins_system-1] << "};" << endl;
      
      cout << "MUID_TRIG_N[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield err
    	{
    	  cout << muid_trig[1][i] << ", ";
    	}
      cout <<  muid_trig[1][bins_system-1] << "};" << endl;

      cout << "MUID_TRIG_S[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield err
    	{
    	  cout << muid_trig[0][i] << ", ";
    	}
      cout <<  muid_trig[0][bins_system-1] << "};" << endl;

      cout << "ACC_RECO_N[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield err
    	{
    	  cout << acc_reco[1][i] << ", ";
    	}
      cout <<  acc_reco[1][bins_system-1] << "};" << endl;

      cout << "ACC_RECO_S[" << bins_system << "] = {";
      for(int i = 0;i < bins_system-1; i++) // south AA inv yield err
    	{
    	  cout << acc_reco[0][i] << ", ";
    	}
      cout <<  acc_reco[0][bins_system-1] << "};" << endl;

      
    } // print out
  
  

  /////////////////////////////////////////////////// make plots


  std::string gr_legend[2][3]  = {"RAA_{pAu}, South","RAA_{pAl}, South","RAA_{^{3}HeAu}, South",
				  "RAA_{pAu}, North","RAA_{pAl}, North","RAA_{^{3}HeAu}, North"};
  
  std::string gr_legend_dAu[2] = {"RAA_{dAu}, South",
				  "RAA_{dAu}, North"};
  
  std::string gr_title[3] = {"Nuclear Mod. Factor vs. pT, Run15pAu",
			     "Nuclear Mod. Factor vs. pT, Run15pAl",
			     "Nuclear Mod. Factor vs. pT, Run14^{3}HeAu"};
  
  std::string gr_title_dAu[2][3] = {
    {"R_{AB} vs. pT, Run15pAu and Run8dAu, South","R_{AB} vs. pT, Run15pAl and Run8dAu, South","R_{AB} vs. pT, Run14HeAu and Run8dAu, South"},
    {"R_{AB} vs. pT, Run15pAu and Run8dAu, North","R_{AB} vs. pT, Run15pAl and Run8dAu, North","R_{AB} vs. pT, Run14HeAu and Run8dAu, North"}};
  
  std::string y_title[3] = {"RAB_{pAu}", "RAB_{pAl}", "RAB_{^{3}HeAu}"};
  
  // Centrality text box
  TLatex l3;
  l3.SetTextSize(0.06);
  l3.SetTextAlign(13);
  l3.SetTextColor(4);
  
  char text4[100];
  if(i_cent == 1)
    sprintf(text4,"0-20 Centrality");
  if(i_cent == 2)
    sprintf(text4,"20-40 Centrality");
  if(i_cent == 3 && i_system == 0)
    sprintf(text4,"40-60 Centrality");
if(i_cent == 3 && i_system == 1)
    sprintf(text4,"40-72 Centrality");
if(i_cent == 3 && i_system == 2)
    sprintf(text4,"40-88 Centrality");
  if(i_cent == 4)
    sprintf(text4,"60-84 Centrality");
  if(i_cent == 0)
    sprintf(text4,"0-100 Centrality");
  l3.SetTextAlign(12);

  TCanvas *c1 = new TCanvas("c1","title",200,10,700,500);
  c1->SetGrid();
  
  TGraphErrors *gr[2];
  
  gr[1] = new TGraphErrors(data_points[i_cent][i_system],pt_array[i_cent][i_system],RAA_array[i_cent][1],pt_err,RAA_ERR_1[i_cent][1]);  // RAA North
  gr[0] = new TGraphErrors(data_points[i_cent][i_system],pt_array[i_cent][i_system],RAA_array[i_cent][0],pt_err,RAA_ERR_1[i_cent][0]); // RAA South
  
  // Plot for North and South pAu,pAl,HeAu
  gr[1]->SetTitle(gr_title[i_system].c_str());
  c1->cd(1);  
  gPad->SetLeftMargin(0.11); // 0.3
  gPad->SetBottomMargin(0.15);
  gr[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr[1]->GetXaxis()->SetLabelSize(0.06);
  gr[1]->GetXaxis()->SetTitleSize(0.07);
  gr[1]->GetXaxis()->SetTitleOffset(0.9);
  gr[1]->GetYaxis()->SetTitle(y_title[i_system].c_str());
  gr[1]->GetYaxis()->SetLabelSize(0.06);
  gr[1]->GetYaxis()->SetTitleSize(0.07);
  gr[1]->GetYaxis()->SetTitleOffset(0.52);
  gr[1]->SetMarkerColor(kBlue); 
  gr[1]->SetMarkerSize(1.0);
  gr[1]->SetMarkerStyle(20); 
  gr[1]->GetXaxis()->SetNdivisions(404);
  gr[1]->GetYaxis()->SetNdivisions(404);
  gr[1]->Draw("AP"); 
  gr[1]->GetYaxis()->SetRangeUser(0.0,2.5);
      
  gr[0]->SetMarkerColor(kViolet);  
  gr[0]->SetMarkerSize(1.0);
  gr[0]->SetMarkerStyle(20); 
  gr[0]->Draw("P");

  ///////////
  double box_pt = 0.1;
  int box_color = 0; 
  
  double sys_frac[3][2][28] = {
    {0.0709831, 0.0710421, 0.0711406, 0.0713135, 0.0715864, 0.0719718, 0.0724668, 0.0730594, 0.0737375, 0.0744808, 0.0752404, 0.075946, 0.0765415, 0.0770054, 0.0773514, 0.0776176, 0.0778499, 0.0780881, 0.0783553, 0.078651, 0.0789504, 0.0792134, 0.0793872, 0.0794182, 0.0792815, 0.0789795, // pAu S 0-100
     0.0709831, 0.0710421, 0.0711406, 0.0713135, 0.0715864, 0.0719718, 0.0724668, 0.0730594, 0.0737375, 0.0744808, 0.0752404, 0.075946, 0.0765415, 0.0770054, 0.0773514, 0.0776176, 0.0778499, 0.0780881, 0.0783553, 0.078651, 0.0789504, 0.0792134, 0.0793872, 0.0794182, 0.0792815, 0.0789795}, // pAu N 0-100
    {0.0472466, 0.0472316, 0.0472529, 0.0473817, 0.0476657, 0.048109, 0.0486986, 0.0494277, 0.0503036, 0.0513228, 0.0524435, 0.053594, 0.054713, 0.0557731, 0.0567573, 0.0576239, 0.0585511, 0.0591197, 0.0598549, 0.0618512,
     {0.0402385, 0.0402491, 0.0402669, 0.0402987, 0.040353, 0.0404405, 0.0405726, 0.040758, 0.0409971, 0.041278, 0.0415797, 0.0418787, 0.0421559, 0.0423997, 0.0426071, 0.0427816, 0.0429867, 0.0432041, 0.0436274, 0.0440612},
     {0.0742518, 0.0739469, 0.073652, 0.0733786, 0.0731536, 0.0730028, 0.0729377, 0.0729583, 0.0730468, 0.0731673, 0.0732946, 0.0734158, 0.0735135, 0.07357, 0.0735867, 0.0735852, 0.0735809, 0.073622, 0.0734297,
      {0.0387342, 0.0387325, 0.0387309, 0.0387293, 0.038727, 0.0387237, 0.0387197, 0.0387166, 0.0387163, 0.0387198, 0.0387258, 0.0387317, 0.0387355, 0.038737, 0.0387364, 0.0387347, 0.0387326, 0.0387298, 0.0387233}};
     
       for(int i = 0; i < data_points[i_cent][i_system];i++)
    {
      draw_box_N(pt_array[i_cent][i_system][i], RAA_array[i_cent][1][i], box_pt, sys_frac[i_system][1][i]*RAA_array[i_cent][1][i], box_color);          
      draw_box_S( pt_array[i_cent][i_system][i], RAA_array[i_cent][0][i], box_pt, sys_frac[i_system][0][i]*RAA_array[i_cent][0][i], box_color); 
    }  
     
  l3.DrawLatexNDC(0.65, 0.92, text4); //4.4,150  

  TLegend *leg = new TLegend(0.11, 0.7, 0.4, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.035);
  leg->AddEntry(gr[1], gr_legend[1][i_system].c_str() , "p"); 
  leg->AddEntry(gr[0], gr_legend[0][i_system].c_str(), "p"); 
  leg->Draw();
 
  epp1->SetLineStyle(6);
  epp1->SetLineColor(kBlue);
  epp1->SetLineWidth(2);

  epp0->SetLineStyle(6);
  epp0->SetLineColor(kBlack);
  epp0->SetLineWidth(2);

  if(i_system == 0)
    {
      // epp1->Draw("SAME");
      //epp0->Draw("SAME");
    }
  
  ////////////////////////////////////Plot for i_system and dAu
  TGraphErrors *gr_dAu[2];

  TCanvas *c[2];
  std::string c_title[2] = {"c[1]","c[0]"}; 

  TLegend *leg_dAu[2];
  
  for(int arm = 0; arm < 2; arm++)
    {
      c[arm] = new TCanvas(c_title[arm].c_str(),"title",200,10,700,500);
      gr_dAu[arm] = new TGraphErrors(23,pt_array_dAu,RAA_dAu[arm][i_cent],pt_err,RAA_dAu_ERR_1[arm][i_cent]); 
      gr_dAu[arm]->SetTitle(gr_title_dAu[arm][i_system].c_str());
      
      c[arm]->SetGrid();
      c[arm]->cd(1);  
	  
      gPad->SetLeftMargin(0.11); // 0.3
      gPad->SetBottomMargin(0.15);
      gr_dAu[arm]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      gr_dAu[arm]->GetXaxis()->SetLabelSize(0.06);
      gr_dAu[arm]->GetXaxis()->SetTitleSize(0.07);
      gr_dAu[arm]->GetXaxis()->SetTitleOffset(0.9);
      gr_dAu[arm]->GetYaxis()->SetTitle("RAB");
      gr_dAu[arm]->GetYaxis()->SetLabelSize(0.06);
      gr_dAu[arm]->GetYaxis()->SetTitleSize(0.07);
      gr_dAu[arm]->GetYaxis()->SetTitleOffset(0.52);
      gr_dAu[arm]->SetMarkerColor(kRed); 
      gr_dAu[arm]->SetMarkerSize(1.0);
      gr_dAu[arm]->SetMarkerStyle(20); 
      gr_dAu[arm]->GetXaxis()->SetNdivisions(404);
      gr_dAu[arm]->GetYaxis()->SetNdivisions(404);
      gr_dAu[arm]->Draw("AP"); 
      gr_dAu[arm]->GetYaxis()->SetRangeUser(0.0,2.5);
      
      gr[arm]->Draw("P");
      l3.DrawLatexNDC(0.65, 0.92, text4); //4.4,150
      
      for(int i = 0; i < 23; i++)
	{
	  if(arm == 0)
	    draw_box_black(pt_array_dAu[i], RAA_dAu[0][i_cent][i], box_pt, RAA_dAu_frac_sys[0][i_cent][i]*RAA_dAu[0][i_cent][i], box_color);
	  if(arm == 1)
	    draw_box_dAu(pt_array_dAu[i], RAA_dAu[1][i_cent][i], box_pt, RAA_dAu_frac_sys[1][i_cent][i]*RAA_dAu[1][i_cent][i], box_color);
	}

      for(int i = 0; i < data_points[i_cent][i_system]; i++)
	{
	  if(arm == 1)
	    draw_box_N(pt_array[i_cent][i_system][i], RAA_array[i_cent][1][i], box_pt, sys_frac[i_system][1][i]*RAA_array[i_cent][1][i], box_color); 
	  else
	    draw_box_S( pt_array[i_cent][i_system][i], RAA_array[i_cent][0][i], box_pt, sys_frac[i_system][0][i]*RAA_array[i_cent][0][i], box_color); 
	}
     
      leg_dAu[arm] = new TLegend(0.11, 0.7, 0.4, 0.9);
      leg_dAu[arm]->SetFillColor(0);
      leg_dAu[arm]->SetTextSize(0.035);
      leg_dAu[arm]->AddEntry(gr[arm], gr_legend[arm][i_system].c_str(), "p"); 
      leg_dAu[arm]->AddEntry(gr_dAu[arm], gr_legend_dAu[arm].c_str(), "p"); 
      leg_dAu[arm]->Draw();	     
    }
  

  bool dAu = true;
  if(dAu)
    {
      TCanvas *c3 = new TCanvas("c3","title",200,10,700,500);
      //c3->cd(1);
      c3->SetGrid();

      gr_dAu[0]->SetTitle("R_{dAu} North and South");
      gr_dAu[0]->SetMarkerColor(kBlack);
      gr_dAu[0]->Draw("AP");
      gr_dAu[1]->Draw("P");

      for(int i = 0; i < 23; i++)
	{
	  draw_box_black(pt_array_dAu[i], RAA_dAu[0][i_cent][i], box_pt, RAA_dAu_frac_sys[0][i_cent][i]*RAA_dAu[0][i_cent][i], box_color);
	  draw_box_dAu(pt_array_dAu[i], RAA_dAu[1][i_cent][i], box_pt, RAA_dAu_frac_sys[1][i_cent][i]*RAA_dAu[1][i_cent][i], box_color);
	}
    }
  l3.DrawLatexNDC(0.65, 0.92, text4); //4.4,150



}// end void macro
