
///////////////////// For p+Au,p+Al,He+Au /////////////////////////
//// NEED TO FILL IN BIAS CORRECTION FACTORS AND N_COLL FOR CENTRALITY STUDY

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

void raa_macro_ALL()
{
  ///////////////////////////////////////////////// dAu information
 
  int bins = 26;

  int i_system;
  cout << "Enter the system ('0' for Run15pAu, '1' for Run15pAl, '2' for Run14HeAu or '3' for Run08dAu)" << endl;
  cin >> i_system;
  int i_cent;
  cout << "Enter the centrality range ('0' for no centrality, '1' for 0-20, '2' for 20-40, '3' for 40-60, '4' for 60-84)" << endl;
  cin >> i_cent;

  double RAA_dAu[2][5][28] = {
    {0.759,0.772,0.853,0.899,0.955,0.934,0.988,1.000,1.043,1.182,1.159,1.161,1.150,1.059,1.234,1.043,1.285,1.133,1.556,1.265,1.186,1.227,1.228,1.643,0.0,0.0,0.0,0.0,
     0.702,0.760,0.840,0.894,0.954,1.008,1.000,1.058,1.058,1.252,1.240,1.249,1.230,1.088,1.406,1.131,1.244,1.209,1.824,1.508,1.182,1.171,1.848,2.079,0.0,0.0,0.0,0.0,
     0.706,0.739,0.815,0.925,0.984,0.874,1.021,1.008,1.067,1.169,1.103,1.137,1.160,1.099,1.234,0.957,1.376,1.077,1.401,1.028,1.481,1.188,0.793,0.996,0.0,0.0,0.0,0.0,
     0.953,0.865,0.995,0.920,0.976,0.929,0.979,1.018,1.015,1.159,1.227,1.036,1.046,1.114,1.152,1.124,1.313,1.427,1.080,1.365,0.985,1.354,1.297,4.130,0.0,0.0,0.0,0.0,
     0.829,0.823,0.922,0.965,0.947,0.965,0.990,0.888,1.097,1.073,1.047,1.175,1.132,0.908,0.926,0.804,1.277,0.642,1.827,0.700,1.193,3.141,1.122,0.443,0.0,0.0,0.0,0.0}, // ARM 0   0-100,0-20,20-40,40-60,60-84
    {0.693,0.690,0.664,0.659,0.652,0.671,0.739,0.703,0.732,0.772,0.764,0.821,0.844,0.716,0.765,0.786,1.032,1.030,1.118,1.047,0.836,0.984,1.122,1.275,0.0,0.0,0.0,0.0,
     0.566,0.557,0.557,0.563,0.559,0.594,0.634,0.587,0.698,0.691,0.691,0.713,0.812,0.631,0.708,0.772,0.922,0.834,0.953,0.922,0.735,0.997,0.971,0.989,0.0,0.0,0.0,0.0,
     0.649,0.735,0.680,0.697,0.661,0.653,0.712,0.701,0.716,0.818,0.831,0.831,0.844,0.728,0.752,0.819,1.191,1.252,1.276,1.206,0.859,0.960,1.111,1.519,0.0,0.0,0.0,0.0,
     0.897,0.731,0.816,0.778,0.802,0.862,0.916,0.865,0.859,0.854,0.781,0.933,0.907,0.807,0.878,0.868,1.109,0.911,1.082,1.246,1.114,1.417,1.268,0.747,0.0,0.0,0.0,0.0,
     1.038,1.100,0.961,0.873,0.859,0.800,1.028,0.963,0.864,0.939,0.966,1.146,0.999,0.886,0.897,0.741,0.985,1.228,1.477,0.906,0.913,1.622,1.381,1.527,0.0,0.0,0.0,0.0}}; // ARM 1 0-100,0-20,20-40,40-60,60-84

  double RAA_dAu_ERR_1[2][5][28] = {
    {0.1, 0.094, 0.092, 0.06, 0.048, 0.048, 0.048, 0.049, 0.054, 0.061, 0.059, 0.067, 0.075, 0.076, 0.101, 0.098, 0.15, 0.152, 0.252, 0.256, 0.303, 0.518, 0.537, 1.036, 0.0, 0.0, 0.0, 0.0,
     0.121,0.104,0.105,0.076,0.059,0.066,0.061,0.066,0.062,0.072,0.076,0.083,0.095,0.093,0.13,0.124,0.165,0.186,0.314,0.329,0.33,0.522,0.83,1.364,0.0,0.0,0.0,0.0,
     0.123,0.108,0.113,0.069,0.051,0.059,0.054,0.054,0.06,0.072,0.073,0.084,0.101,0.102,0.131,0.121,0.195,0.19,0.28,0.266,0.496,0.641,0.469,0.836,0.0,0.0,0.0,0.0,
     0.103,0.104,0.094,0.06,0.055,0.056,0.067,0.06,0.074,0.091,0.086,0.088,0.101,0.114,0.14,0.153,0.212,0.253,0.263,0.356,0.345,0.796,0.671,2.905,0.0,0.0,0.0,0.0,
     0.14,0.091,0.074,0.068,0.057,0.059,0.058,0.059,0.077,0.084,0.09,0.107,0.122,0.115,0.138,0.136,0.234,0.187,0.401,0.251,0.42,1.737,0.714,0.0,0.0,0.0,0.0,0.0},      // 0-100,0-20,20-40,40-60,60-84
    {0.064,0.059,0.048,0.037,0.032,0.032,0.029,0.029,0.031,0.035,0.037,0.043,0.049,0.048,0.056,0.07,0.107,0.129,0.166,0.186,0.195,0.298,0.387,0.83,0.0,0.0,0.0,0.0,
     0.054,0.046,0.034,0.031,0.029,0.029,0.031,0.029,0.034,0.037,0.041,0.047,0.056,0.052,0.064,0.081,0.112,0.124,0.163,0.183,0.192,0.396,0.378,0.739,0.0,0.0,0.0,0.0,
     0.076,0.066,0.065,0.05,0.042,0.035,0.037,0.037,0.039, 0.047,0.052,0.057,0.065,0.066,0.075,0.095,0.148,0.192,0.223,0.247,0.237,0.356,0.448,1.038,0.0,0.0,0.0,0.0,
     0.098,0.087,0.093,0.046,0.042,0.046,0.045,0.047,0.051,0.056,0.057,0.072,0.079,0.081,0.096,0.115,0.161,0.167,0.23,0.281,0.343,0.632,0.55,0.0,0.0,0.0,0.0,0.0,
     0.123,0.089,0.064,0.065,0.05,0.051,0.054,0.058,0.062,0.068,0.076,0.094,0.097,0.097,0.112,0.116,0.169,0.227,0.305,0.257,0.309,0.816,0.628,1.12,0.0,0.0,0.0,0.0}};// none,0-20,20-40,40-60,60-84

  double RAA_dAu_ERR_2[2][5][28] = {
    {0.053,0.054,0.059,0.062,0.065,0.064,0.068,0.069,0.072,0.081,0.08,0.08,0.079,0.073,0.085,0.072,0.089,0.078,0.108,0.089,0.083,0.085,0.085,0.0,0.0,0.0,0.0,0.0,
     0.049,0.053,0.058,0.062,0.065,0.069,0.069,0.073,0.073,0.086,0.085,0.086,0.085,0.075,0.097,0.078,0.086,0.084,0.127,0.106,0.082,0.082,0.128,0.146,0.0,0.0,0.0,0.0,
     0.049,0.042,0.057,0.064,0.068,0.060,0.07,0.069,0.073,0.08,0.076,0.078,0.08,0.076,0.085,0.066,0.095,0.075,0.097,0.072,0.103,0.083,0.055,0.07,0.0,0.0,0.0,0.0,
     0.066,0.061,0.069,0.063,0.067,0.064,0.067,0.07,0.07,0.08,0.085,0.071,0.072,0.077,0.08,0.078,0.091,0.099,0.075,0.096,0.069,0.094,0.09,0.291,0.0,0.0,0.0,0.0,
     0.057,0.058,0.064,0.066,0.065,0.066,0.068,0.061,0.075,0.074,0.072,0.081,0.078,0.063,0.064,0.056,0.088,0.044,0.127,0.049,0.083,0.219,0.078,0.031,0.0,0.0,0.0,0.0}, /// ARM 0
    {0.052,0.052,0.05,0.049,0.048,0.05,0.055,0.052,0.054,0.057,0.057,0.061,0.063,0.053,0.057,0.059,0.077,0.077,0.084,0.079,0.063,0.074,0.084,0.0,0.0,0.0,0.0,0.0,
     0.042,0.042,0.042,0.042,0.042,0.044,0.047,0.044,0.052,0.051,0.052,0.053,0.061,0.047,0.053,0.058,0.069,0.062,0.071,0.07,0.055,0.075,0.073,0.075,0.0,0.0,0.0,0.0,
     0.049,0.056,0.051,0.052,0.049,0.048,0.053,0.052,0.053,0.061,0.062,0.062,0.063,0.054,0.056,0.061,0.89,0.094,0.096,0.091,0.065,0.072,0.083,0.115,0.0,0.0,0.0,0.0,
     0.067,0.055,0.061,0.058,0.06,0.064,0.068,0.064,0.064,0.064,0.058,0.069,0.068,0.06,0.066,0.065,0.083,0.068,0.081,0.094,0.084,0.107,0.095,0.057,0.0,0.0,0.0,0.0,
     0.078,0.083,0.072,0.065,0.064,0.059,0.076,0.072,0.064,0.07,0.072,0.085,0.075,0.066,0.067,0.055,0.074,0.092,0.111,0.068,0.069,0.122,0.103,0.116,0.0,0.0,0.0,0.0}}; //ARM 1

  double RAA_dAu_frac_sys[2][5][28];

  // calculate the fractional systematic uncertainties for dAu
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0;i < bins-3; i++)
	{
	  RAA_dAu_frac_sys[arm][i_cent][i] =  RAA_dAu_ERR_2[arm][i_cent][i]/RAA_dAu[arm][i_cent][i];
	  // cout << "Arm " << arm << RAA_dAu_sys: " << RAA_dAu_frac_sys[arm][i_cent][i] << endl;
	}
    }

  /// UNCOMMENT WHEN PLOTTING pAu RAA
      /*
  cout << "RAA_dAu_frac_sys[" << bins-3 << "] = {";
  for(int arm = 0; arm < 2; arm++)
    {  
      for(int i=0;i<bins-4;i++)
	{
	  cout << RAA_dAu_frac_sys[arm][i_cent][i] << ", ";
	}
      cout <<  RAA_dAu_frac_sys[arm][i_cent][bins-4] << "};" << endl;
    }

  cout << "RAA_dAu_frac_sys[" << bins-3 << "] = {";
  for(int arm = 0; arm < 2; arm++)
    { 
      for(int i=0;i<bins-4;i++)
	{
	  cout << RAA_dAu_frac_sys[arm][i_cent][i] << ", ";
	}
    }
  cout <<  RAA_dAu_frac_sys[arm][i_cent][bins-4] << "};" << endl;
      */

  ////////

  double pt_array_dAu[28] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.25,5.75,6.50,0.0,0.0,0.0,0.0,0.0};
  
  ///////////////////////////////////////////  end dAu information  
  
  int bin_width;
  double rA[2] = {1.0,1.0}; // initialize to 1.  Will be overwritten when binshift is set to true
  bool write = false;
  bool binshift = true; // run macro first with this set to false to generate new cross section results.  Then run binshift_Tgraph_ALL.C in pp directory to get pAu binshift corrections written.  Then turn this to true to read in the new files.
  bool print_screen = true;

  // [i_cent][i_system]
  double N_coll[5][3] = {4.667,2.10,10.4, // pAu,pAl,HeAu no centrality (FILL IN LATER)
			 8.2,1.1,1.1,
			 6.1,1.1,1.1,
			 4.4,1.1,1.1, 
			 2.6,1.1,1.1}; 
  
  double bias_correction[5][3] = {0.858,0.80,0.89, // no centrality pAu,pAl,HeAu (FILL IN LATER)
				  0.90,0.98,1.02, // 020
				  1.1,1.1,1.1, // 2040 (FILL IN LATER)
				  1.1,1.1,1.1,  // 4060  (FILL IN LATER)
				  1.1,1.1,1.1};  // 6084  (FILL IN LATER)
  
  int data_points[5][3] = {26,20,19, // no centrality (pAu,pAl,HeAu)
			   19,19,19, // 020
			   19,19,19, 
			   19,19,19, 
			   19,19,19}; // 6084

  // [i_cent,i_system,data_points]
  double pt_array[5][3][26] = {
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.25,6.75,     
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,5.5,6.5,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // no centrality 
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},   // 0-20  (pAu, pAl and HeAu share same binning)
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},   // 20-40
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, //40-60
    {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; //60-84

  double pt_width[5][3][28] = {  
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.0,0.0, // Run15pAu
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // Run15pAl 
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // Run14HeAu    NO CENTRALITY
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  // cent 0-20
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // cent 20-40
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // cent 40-60
    {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; // cent 60-84
  
  double invariant_yield;
  double two_pi = 2*3.141592653;
  double eff_BBC = 1.0;  
  if(i_cent == 7)
    eff_BBC = 0.55;
  double eff_jpsi = 1.0;
  if(i_cent == 7)
    eff_jpsi = 1/0.79;

  double muid_trig_eff;
  double d_pt = 0.25; // pt bin width
  double d_rap = 1.0; // rapidity bin width
  double pt_bin_center = 0.125;
  double denom;
  double acc_reco_eff;
  double jpsi_counts[2];
  double jpsi_err_1[2];
  double jpsi_err_2[2];
  double exp = pow(10,-9);
  double ppg_cross = 0.0;
  double ppg_err_1 = 0.0;
  double ppg_err_2 = 0.0;

  // [i_cent][arm][i_system]
  double MB[5][2][3] = { 2.10268903808*pow(10,11),2.04722*pow(10,11),2.087*pow(10,10), // no centrality (pAu,pAl,HeAu)
			 1.9356537088*pow(10,11),2.00743*pow(10,11),3.32923*pow(10,10),  // 

			 2.10268903808*pow(10,11)*0.2,1.1,1.1, // cent 0-20
			 1.9356537088*pow(10,11)*0.2,1.1,1.1,  
			 
			 2.10268903808*pow(10,11)*0.2,1.1,1.1, // cent 20-40
			 1.9356537088*pow(10,11)*0.2,1.1,1.1,  
			 
			 2.10268903808*pow(10,11)*0.2,1.1,1.1,  // cent 40-60
			 1.9356537088*pow(10,11)*0.2,1.1,1.1,  
			 
			 2.10268903808*pow(10,11)*0.24,1.1,1.1,  // cent 60-84
			 1.9356537088*pow(10,11)*0.2,1.1,1.1};  

  double invariant_yield_err;
  double r_N;
  double r_S;
  double p0,p1;
  double f0,f1;
  double sum[2] = {0.0,0.0};
  double large_errors[2] = {0.0,0.0};
  double percent_error[2];

  TF1 *a[2][3]; 
  TF1 *t[2][3];
 
  double pt_err[28];
  double run15[2][28];
  double run15_err[2][28];
  double pp_inv_corrected[2][28];
  
  double ratio_array[28]; // Run15pp N/S
  double ratio_err_1[28];
  double ave_array[28];
  double ave_array_ERR_1[28];
  double ratio[28];
  double ratio_ERR_1[28];

  double RAA_array[2][28];
  double RAA_ERR_1[2][28]; 
  double pp_inv_yield[2][28];
  double pp_inv_yield_err[2][28];
    
  double tmp_N;
  double tmp_S;
  double err_N;
  double pt;
  
  TFile *file_data_acc;
  TFile *file_data_trig;

  std::string acc_filename[3] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAu200_acceff_pT_InclusiveJpsi_20180914.root",
				 // "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0020_acceff_pT_InclusiveJpsi_EMBED_C.root",
				 "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAl200_acceff_pT_InclusiveJpsi_20180914.root",
				 "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run14HeAu200_acceff_pT_InclusiveJpsi_20180914.root"};	
   
  std::string obj_acc_filename[2][3] = {"frecoeff_pT_pAu_arm0","frecoeff_pT_pAl_arm0","frecoeff_pT_HeAu_arm0",
					"frecoeff_pT_pAu_arm1","frecoeff_pT_pAl_arm1","frecoeff_pT_HeAu_arm1"};
   
  std::string system_name[3] = {"p + Au","p + Al","He + Au"};

  //fill contents in a[arm]
  for(int arm = 0; arm < 2; arm++)
    { 
      for(int j_system = 0; j_system < 3; j_system++)
	{
	  file_data_acc = TFile::Open(acc_filename[j_system].c_str()); 
	  file_data_acc->GetObject(obj_acc_filename[arm][j_system].c_str(),a[arm][j_system]);
	}
    }

  std::string trig_filename[3] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAu200_trigeff_MUID2D_pT_output_20180913.root",
				  // "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_corrbg/Run15pAu200_CENT0020_trigeff_MUID2D_pT_output_EMBED_C.root",
				  "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAl200_trigeff_MUID2D_pT_output_20180913.root",
				  "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run14HeAu200_trigeff_MUID2D_pT_output_20180913.root"};
  
  std::string trig_obj_filename[2][3] = {"ftrigeff_pT_dAu_arm0","ftrigeff_pT_pAl_arm0","ftrigeff_pT_HeAu_arm0",
  					 "ftrigeff_pT_dAu_arm1","ftrigeff_pT_pAl_arm1","ftrigeff_pT_HeAu_arm1"};
  
  //fill contents in t[arm]
  for(int arm = 0; arm < 2; arm++)
    { 
      for(int j_system = 0; j_system < 3; j_system++)
	{
	  file_data_trig = TFile::Open(trig_filename[j_system].c_str());
	  file_data_trig->GetObject(trig_obj_filename[arm][j_system].c_str(),t[arm][j_system]);
	}
    }
  
  ///////////////////////////////////////////////////  
  
  //number of groups of data(centrality),rows per group(arm),elements per row (system type)
  char basename[30][800] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAu_S_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAl_S_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run14HeAu_S_bestfit_parameters_tony_1_", // no centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAu_N_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run15pAl_N_bestfit_parameters_tony_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/Run14HeAu_N_bestfit_parameters_tony_1_", // no centrality, arm 1
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/020_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/020_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 0-20 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrrbg/020_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/020_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/020_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 0-20 centrality, arm 1
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 20-40 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/2040_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 20-40 centrality, arm 1
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 40-60 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run15pp_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/4060_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_", // 40-60 centrality, arm 1
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run15pAu_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run15pAl_S_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run14HeAu_S_bestfit_parameters_tony_cent_1_", // 60-84 centrality, arm 0
			    "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run15pAu_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run15pAl_N_bestfit_parameters_tony_cent_1_","/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/6084_Centrality/Run14HeAu_N_bestfit_parameters_tony_cent_1_"}; // 60-84 centrality, arm 1
   
  // for binshift correction files (MAY NEED TO REVISE)
  char basename2[30][800] = {"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // no centrality, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // no centrality, arm 1
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 020, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/020_Centrality/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 020, arm 1
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 2040, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 2040, arm 1
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 4060, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_N_correction_r_", // 4060, arm 1
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_S_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_S_correction_r_", // 6084, arm 0
			     "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_binshift_corrections/Run15pAu_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_binshift_corrections/Run15pAl_N_correction_r_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_binshift_corrections/Run14HeAu_N_correction_r_"}; // 6084, arm 1 
	
  // for pp inv. cross section 
  char basename3[6][800] = {"/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_cross_sections/Run15pp_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_cross_sections/Run15pp_cross_S_",
			    "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_cross_sections/Run15pp_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_cross_sections/Run15pp_cross_N_"};

  // [arm][i_system] for AA cross sections
  char basename4[6][500] = {"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_cross_sections/Run15pAu_cross_S_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_cross_sections/Run15pAl_cross_S_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_cross_sections/Run14HeAu_cross_S_",
			    "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAu_cross_sections/Run15pAu_cross_N_","/phenix/hhj/klsmith/Jpsis/Run15pAl/fit_NOFVTX/krista_fits/tony_bestfit_parameters/pAl_cross_sections/Run15pAl_cross_N_","/phenix/hhj/klsmith/Jpsis/Run14HeAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/HeAu_cross_sections/Run14HeAu_cross_N_"};

  vector < vector < vector <std::string> > > pp_filename;
  vector < vector < vector <std::string> > > binshifted_cross_filename;
  
  for(int arm = 0; arm < 2; arm++)
    {
      vector < vector <std::string> > filename_type;
      vector < vector <std::string> > filename_type2;
      
      for(int j_system = 0; j_system < 3; j_system++)
	{
	  vector <std::string> filename_pt;
	  vector <std::string> filename_pt2;
	  
	  int j = arm*3 + j_system;
	  cout << "j: " << j << endl;
	  
	  for(int i = 0; i < data_points[i_cent][j_system]; i++)
	    {
	      char name[500];
	      sprintf(name,"%s%i%s",basename3[j],i+1,".dat");
	      std::string filename_tmp = name;
	      filename_pt.push_back(filename_tmp);
	      
	      char name2[500];
	      sprintf(name2,"%s%i%s",basename4[j],i+1,".dat");
	      std::string filename_tmp2 = name2;
	      filename_pt2.push_back(filename_tmp2);
	    }
	  filename_type.push_back(filename_pt);
	  filename_type2.push_back(filename_pt2);
	}
      pp_filename.push_back(filename_type);
      binshifted_cross_filename.push_back(filename_type2);
    }
  
  vector < vector < vector < vector <std::string> > > > bestfit_filename;
  vector < vector < vector < vector <std::string> > > > binshift_r_filename;
  
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
		  sprintf(name,"%s%i%s",basename[j],i+1,".dat");	     
		  std::string filename_tmp = name;
		  filename_pt.push_back(filename_tmp);
		  //cout << "i+1 = " << i+1 << endl;
		  
		  char name2[800];
		  sprintf(name2,"%s%i%s",basename2[j],i+1,".dat");	     
		  std::string filename_tmp2 = name2;
		  filename_pt2.push_back(filename_tmp2);
		} 
	      filename_sys.push_back(filename_pt);
	      filename_sys2.push_back(filename_pt2);
	    }
	  filename_type.push_back(filename_sys);
	  filename_type2.push_back(filename_sys2);
	}
      bestfit_filename.push_back(filename_type);
      binshift_r_filename.push_back(filename_type2);
    } // 
  
  // for(int l = 0; l < 2; l++)
  //   {
  //     for(int k = 0; k < 5; k++)
  //       {
  // 	 for(int i = 0; i < data_points[k][i_system]; i++)
  // 	   {
  // 	     cout << bestfit_filename[k][l][i_system][i].c_str() << endl;
  // 	   }
  //       }
  //   } // prints out the desired system, not all systems (don't have pAl centrality or HeAu centrality files yet)
 

  // BEGIN RAA Calc
  for(int arm = 0; arm < 2; arm++)
    {

      if(i_cent == 0)
	MB[i_cent][arm][i_cent] /= eff_BBC;
		
      for(int i = 0; i < data_points[i_cent][i_system];i++)
	{
	  pt = pt_array[i_cent][i_system][i];
	  pt_err[i] = 0.0;
	  
	  acc_reco_eff = a[arm][i_system]->Eval(pt);
	  muid_trig_eff = t[arm][i_system]->Eval(pt);

	  cout << "muid " << muid_trig_eff << ", acc " << acc_reco_eff << endl;

	  bin_width = pt_width[i_cent][i_system][i]/(0.25);
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
	    return;
	  sum[arm]+= jpsi_counts[arm];
	  
	  percent_error[arm] = jpsi_err_1[arm]/jpsi_counts[arm]*100; 
	  cout << "Arm " << arm << " " << system_name[i_system].c_str() << " jpsi counts " << round(jpsi_counts[arm]) << " , sqrt error " << sqrt(jpsi_counts[arm]) << " , error " << jpsi_err_1[arm] << " and % error " << round(percent_error[arm]) << endl;
	  
	  if(percent_error[arm] > 10.00)
	    large_errors[arm]+= 1;
	  
	  denom = two_pi*pt*d_pt*d_rap*acc_reco_eff*muid_trig_eff*MB[i_cent][arm][i_system];
	  invariant_yield = (jpsi_counts[arm]*eff_jpsi*bias_correction[i_cent][i_system])/denom; 
	  invariant_yield_err = (jpsi_err_1[arm]*eff_jpsi*bias_correction[i_cent][i_system])/denom;
	  
	  invariant_yield /= bin_width;
	  invariant_yield_err /= bin_width;
	  
	  if(binshift == true)
	    {
	      ifstream binshift(binshift_r_filename[i_cent][arm][i_system][i].c_str());
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
		return;
	    }
	  
	  run15[arm][i] = invariant_yield/rA[arm];
      	  run15_err[arm][i] = invariant_yield_err;
	 
    	  std::ifstream Run15pp_cross(pp_filename[arm][i_system][i].c_str(),std::fstream::in); 
	  Run15pp_cross >>  p0  >>  p1; 
	  pp_inv_yield[arm][i] = p0/(42*0.001);
	  pp_inv_yield_err[arm][i] = p1/(42*0.001); 
	  Run15pp_cross.close();
      
	  RAA_array[arm][i] = 1/N_coll[i_cent][i_system]*(run15[arm][i]/pp_inv_yield[arm][i]); 
	  double tmp2 =(pow(1/N_coll[i_cent][i_system]*run15_err[arm][i]/run15[arm][i],2)+pow(1/N_coll[i_cent][i_system]*pp_inv_yield_err[arm][i]/pp_inv_yield[arm][i],2))*(pow(run15[arm][i]/pp_inv_yield[arm][i],2));
	  RAA_ERR_1[arm][i] = sqrt(tmp2);
	
	  if(write == true)
	    {
	      f0 = run15[arm][i]*42*0.001;;
	      f1 = run15_err[arm][i]*42*0.001; 
	      
	      std::fstream binshifted_cross(binshifted_cross_filename[arm][i_system][i].c_str(),std::ofstream::out); 
	      binshifted_cross <<  f0  <<  " " << f1; 
	      binshifted_cross.close();
	    }
	  
	}//i for loop
    }// end arm loop  
  
  cout << "Number of fits with error > 10% in N: " << large_errors[1] << " and S: " << large_errors[0] << endl;
  cout << "Total Jpsi counts North: " << int(sum[1]) << " and South: " << int(sum[0]) << endl;
  
  
///////////////////////////////////////////////////
 
  if(print_screen == true)
    {
      cout << "_______________________________ " << endl;
      cout << "RAA_N_ERR_1[" <<bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	   cout << RAA_ERR_1[1][i] << ", ";
      	}
      cout <<  RAA_ERR_1[1][bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "RAA_S_ERR_1[" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	  cout << RAA_ERR_1[0][i] << ", ";
      	}
      cout <<  RAA_ERR_1[0][bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "pt_array[" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	  cout << pt_array[i_cent][i_system][i] << ", ";
      	}
      cout <<  pt_array[i_cent][i_system][bins-1] << "};" << endl;
      /////////////
      bins = data_points[i_cent][i_system];
      cout << "_______________________________ " << endl;
      cout << "RAA_N_ARRAY[" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
	{
	  cout << RAA_array[1][i] << ", ";
	}
      cout <<  RAA_array[1][bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
      cout << "RAA_S_ARRAY[" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
    	{
    	  cout << RAA_array[0][i] << ", ";
    	}
      cout <<  RAA_array[0][bins-1] << "};" << endl;
      bins = 26; //reset
    } // print out
  
  
  std::string gr_legend[2][3]  = {"RAA_{pAu}, South","RAA_{pAl}, South","RAA_{^{3}HeAu}, South",
				  "RAA_{pAu}, North","RAA_{pAl}, North","RAA_{^{3}HeAu}, North"};
  
  std::string gr_legend_dAu[2] = {"RAA_{dAu}, South",
				  "RAA_{dAu}, North"};
  
  std::string gr_title[3] = {"Nuclear Mod. Factor vs. pT, Run15pAu",
			     "Nuclear Mod. Factor vs. pT, Run15pAl",
			     "Nuclear Mod. Factor vs. pT, Run14^{3}HeAu"};
  
  std::string gr_title_dAu[2][3] = {
    {"Nuclear Mod. Factor vs. pT, Run15pAu and Run8dAu, South","Nuclear Mod. Factor vs. pT, Run15pAl and Run8dAu, South","Nuclear Mod. Factor vs. pT, Run14HeAu and Run8dAu, South"},
    {"Nuclear Mod. Factor vs. pT, Run15pAu and Run8dAu, North","Nuclear Mod. Factor vs. pT, Run15pAl and Run8dAu, North","Nuclear Mod. Factor vs. pT, Run14HeAu and Run8dAu, North"}};
  
  std::string y_title[3] = {"RAA_{pAu}", "RAA_{pAl}", "RAA_{^{3}HeAu}"};
  
  TCanvas *c1 = new TCanvas("c1","title",200,10,700,500);
  c1->SetGrid();
  
  TGraphErrors *gr[2];
  
  gr[1] = new TGraphErrors(data_points[i_cent][i_system],pt_array[i_cent][i_system],RAA_array[1],pt_err,RAA_ERR_1[1]);  // RAA North
  gr[0] = new TGraphErrors(data_points[i_cent][i_system],pt_array[i_cent][i_system],RAA_array[0],pt_err,RAA_ERR_1[0]); // RAA South
  
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
   
  double box_pt = 0.1;
  int box_color = 0; 
  
  double sys_frac[3][2][28] = {
    {0.090,0.090,0.091,0.091,0.091,0.091,0.092,0.092,0.093,0.094,0.094,0.095,0.096,0.096,0.097,0.098,0.098,0.098,0.098,0.098,0.099,0.099,0.1,0.1,0.1,0.1,0.0,0.0, // pAu South no centrality
     0.075,0.075,0.075,0.074,0.074,0.074,0.074,0.074,0.073,0.073,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.074,0.0,0.0}, // pAu North no centrality
    {0.078,0.078,0.078,0.078,0.078,0.078,0.079,0.079,0.081,0.082,0.084,0.084,0.085,0.086,0.087,0.087,0.087,0.088,0.089,0.090,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // pAl South no centrality
     0.076,0.076,0.076,0.077,0.077,0.077,0.077,0.077,0.077,0.077,0.078,0.078,0.078,0.078,0.079,0.079,0.079,0.079,0.079,0.079,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, //pAl North no centrality
    {0.088,0.088,0.088,0.088,0.087,0.087,0.087,0.087,0.088,0.088,0.088,0.089,0.089,0.089,0.089,0.090,0.090,0.090,0.090,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // HeAu South no centrality
     0.065,0.065,0.065,0.065,0.065,0.065,0.065,0.065,0.065,0.065,0.066,0.066,0.066,0.066,0.066,0.066,0.066,0.066,0.066,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};  // HeAu North no centrality

  for(int i = 0; i < data_points[i_cent][i_system];i++)
    {
      draw_box_N(pt_array[i_cent][i_system][i], RAA_array[1][i], box_pt, sys_frac[i_system][1][i]*RAA_array[1][i], box_color);          
      draw_box_S( pt_array[i_cent][i_system][i], RAA_array[0][i], box_pt, sys_frac[i_system][0][i]*RAA_array[0][i], box_color); 
    }       
  
  TLegend *leg = new TLegend(0.11, 0.7, 0.4, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.035);
  leg->AddEntry(gr[1], gr_legend[1][i_system].c_str() , "p"); 
  leg->AddEntry(gr[0], gr_legend[0][i_system].c_str(), "p"); 
  leg->Draw();
 
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
      gr_dAu[arm]->GetYaxis()->SetTitle("RAA");
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
      
      for(int i = 0; i < 23; i++)
	draw_box_dAu(pt_array_dAu[i], RAA_dAu[arm][i_cent][i], box_pt, RAA_dAu_frac_sys[arm][i_cent][i]*RAA_dAu[arm][i_cent][i], box_color);
      for(int i = 0; i < data_points[i_cent][i_system]; i++)
	{
	  if(arm == 1)
	    draw_box_N(pt_array[i_cent][i_system][i], RAA_array[1][i], box_pt, sys_frac[i_system][1][i]*RAA_array[1][i], box_color); 
	  else
	    draw_box_S( pt_array[i_cent][i_system][i], RAA_array[0][i], box_pt, sys_frac[i_system][0][i]*RAA_array[0][i], box_color); 
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

      gr_dAu[0]->Draw("AP");
      gr_dAu[1]->Draw("P");

      for(int arm = 0; arm < 2; arm++)
	{
	  for(int i = 0; i < 23; i++)
	    draw_box_dAu(pt_array_dAu[i], RAA_dAu[arm][i_cent][i], box_pt, RAA_dAu_frac_sys[arm][i_cent][i]*RAA_dAu[arm][i_cent][i], box_color);
	}
    }
  
}// end void macro

