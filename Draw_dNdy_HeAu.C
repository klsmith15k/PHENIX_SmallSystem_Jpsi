#include "Style.h"

void Draw_dNdy_HeAu(){

	const bool bWRITE = false;

	SetStyle();

	const int ncent = 3;
	const int cent_per[ncent+1] = {0, 20, 40, 88};
	const int cent_perv2[ncent+2] = {0, 20, 40, 60, 88};

	const int nColor[ncent] = {2, 4, 6};

	const int narm = 2;
	const int nybin = 4;

	double xx[narm][nybin] = {{-1.325, -1.575, -1.825, -2.075}, {1.325, 1.575, 1.825, 2.075}};
	double xx_err[narm][nybin] = {{0.125, 0.125, 0.125, 0.125}, {0.125, 0.125, 0.125, 0.125}};

	double xx_int[narm] = {-1.7, 1.7};
	double xx_int_err[narm] = {0.1, 0.1};

	const double Ncoll_cent[ncent] = {22.3, 14.8, 5.48};
	const double Ncoll_cent_sys[ncent] = {0.076233184, 0.074324324, 0.080272996};
	const double BiasF_cent[ncent] = {0.950, 1.010, 1.025};
	const double BiasF_cent_sys[ncent] = {0.010526316, 0.0099009901, 0.035018614};
	const double Ncoll_syserr = 0.25;

	const double Njpsi[ncent][narm][nybin] = {
		{{236, 732, 658, 196}, {122, 556, 560, 281}},
		{{133, 450, 423, 146}, {92, 370, 493, 188}},
		{{148, 491, 393, 122}, {116, 389, 454, 179}}
	};

	const double Njpsi_err[ncent][narm][nybin] = {
		{{22, 38, 43, 21}, {14, 29, 34, 26}},
		{{16, 30, 24, 15}, {11, 25, 31, 16}},
		{{13, 30, 31, 14}, {12, 24, 28, 15}}
	};

	const double Njpsi_cent[ncent][narm] = {
		{1814, 1509}, {1158, 1154}, {1128, 1123}
	};

	const double Njpsi_cent_err[ncent][narm] = {
		{79, 51}, {45, 47}, {46, 46}
	};

	const double sys_fit[narm] = {0.030, 0.028};
	const double sys_corrbkg[ncent][narm] = {{0.0095, 0.0297}, {0.0411, 0.0279}, {0.0228, 0.0385}};
	const double sys_var[narm] = {0.050, 0.015};
	const double sys_acceff[narm] = {0.032, 0.037};
	double sys_trigeff_y_int[ncent+1][narm] = {0.}; //will be updated later
	double sys_trigeff_y[ncent+1][narm][nybin] = {0.}; //will be updated later

	//Run15pp200
	const double eff_pp_bbc = 0.55;
	const double eff_pp_trk = 0.79;

	const double Njpsi_pp15[narm][nybin] = {{3151, 12079, 10645, 2978}, {2082, 10583, 12170, 7040}};
	const double Njpsi_pp15_err[narm][nybin] = {{66, 132, 124, 63}, {54, 120, 138, 109}};

	const double Njpsi_pp15_int[narm] = {28511, 31487};
	const double Njpsi_pp15_int_err[narm] = {205, 215};

	const double sys_fit_pp[narm] = {0.011, 0.010};
	const double sys_corrbkg_pp[narm] = {0.020, 0.020};
	const double sys_var_pp[narm] = {0.047, 0.040};
	const double sys_acceff_pp[narm] = {0.054, 0.061};
	double sys_trigeff_y_int_pp[narm] = {0.026, 0.023}; //will be updated later
	double sys_trigeff_y_pp[narm][4] = {{0.010, 0.017, 0.040, 0.073}, {0.012, 0.010, 0.027, 0.054}}; //will be updated later
	const double sys_pol_pp[narm][4] = {{0.138, 0.128, 0.116, 0.101}, {0.129, 0.1333, 0.116, 0.104}};
	const double sys_pol_pp_int[narm] = {0.12, 0.12};

	TFile *infile_pp = new TFile("event_files/Run15pp200_eventhist_for_dimuon_20180402.root","read");

	double NMB_pp[narm];
	TH1F *hBBCZ_pp[narm];
	for (int iarm=0; iarm<narm; iarm++){
		hBBCZ_pp[iarm] = (TH1F*)infile_pp->Get(Form("hBBCZ_arm%d",iarm));
		NMB_pp[iarm] = hBBCZ_pp[iarm]->Integral();
	}
	infile_pp->Close();

	TFile *infile1_pp = new TFile("eff_files/Run15pp200_acceff_InclusiveJpsi_C.root","read");

	TH1F *heff_y_pp[narm];
	TH1F *heff_y_int_pp[narm];
	double eff_y_pp[narm][nybin], eff_y_int_pp[narm];

	for (int iarm=0; iarm<narm; iarm++){
		heff_y_pp[iarm] = (TH1F*)infile1_pp->Get(Form("heff_yy_w_pp_set1_arm%d",iarm));
		heff_y_int_pp[iarm] = (TH1F*)infile1_pp->Get(Form("heff_yy_int_w_pp_set1_arm%d",iarm));

		for (int iy=0; iy<nybin; iy++){
			eff_y_pp[iarm][iy] = heff_y_pp[iarm]->GetBinContent(iy+1);
			cout << "pp ARM: " << iarm << ", Y: " << iy << ", eff0: " << eff_y_pp[iarm][iy] << endl;
		}

		eff_y_int_pp[iarm] = heff_y_int_pp[iarm]->GetBinContent(1);
	}
	infile1_pp->Close();

	TFile *infile2_pp = new TFile("eff_files/Run15pp200_C_trigeff_MUID2D_output.root","read");

	TH1F *htrigeff_pp[narm];
	TH1F *htrigeff_pp_int[narm];
	TH1F *htrigeff_wt_pp[narm];
	TH1F *htrigeff_wt_pp_int[narm];
	double trigeff_y_pp[narm][nybin];
	double trigeff_y_int_pp[narm];

	for (int iarm=0; iarm<narm; iarm++){
		htrigeff_pp[iarm] = (TH1F*)infile2_pp->Get(Form("htrigeff_y_w_pp_set1_arm%d",iarm));
		htrigeff_wt_pp[iarm] = (TH1F*)infile2_pp->Get(Form("htrigeff_y_wt_w_pp_set1_arm%d",iarm));
		htrigeff_pp_int[iarm] = (TH1F*)infile2_pp->Get(Form("htrigeff_y_int_w_pp_set1_arm%d",iarm));
		htrigeff_wt_pp_int[iarm] = (TH1F*)infile2_pp->Get(Form("htrigeff_y_int_wt_w_pp_set1_arm%d",iarm));

		for (int iy=0; iy<nybin; iy++){
			float eff0 = htrigeff_pp[iarm]->GetBinContent(iy+1);
			float eff1 = htrigeff_wt_pp[iarm]->GetBinContent(iy+1);

			trigeff_y_pp[iarm][iy] = (eff0 + eff1)/2;
			sys_trigeff_y_pp[iarm][iy] = fabs(eff0 - eff1)/(eff0 + eff1);

			cout << "pp ARM: " << iarm << ", Y: " << iy << ", TRIG EFF: " << trigeff_y_pp[iarm][iy] << ", Err: " << sys_trigeff_y_pp[iarm][iy] << endl;
		}

		float eff0 = htrigeff_pp_int[iarm]->GetBinContent(1);
		float eff1 = htrigeff_wt_pp_int[iarm]->GetBinContent(1);

		trigeff_y_int_pp[iarm] = (eff0 + eff1)/2;
		sys_trigeff_y_int_pp[iarm] = fabs(eff0 - eff1)/(eff0 + eff1);
	}
	infile2_pp->Close();

	double Xec_pp15[narm][nybin], Xec_pp15_err[narm][nybin], Xec_pp15_syserr[narm][nybin];
	double Y_pp15[narm][nybin], Y_pp15_err[narm][nybin], Y_pp15_syserr[narm][nybin];
	double Xec_pp15_int[narm], Xec_pp15_int_err[narm], Xec_pp15_int_syserr[narm];
	double Y_pp15_int[narm], Y_pp15_int_err[narm], Y_pp15_int_syserr[narm];

	for (int iarm=0; iarm<narm; iarm++){
		for (int iy=0; iy<nybin; iy++){
			Y_pp15[iarm][iy] = (Njpsi_pp15[iarm][iy]/eff_pp_trk) / (NMB_pp[iarm]/eff_pp_bbc) / 0.25 / eff_y_pp[iarm][iy] / trigeff_y_pp[iarm][iy];
			Y_pp15_err[iarm][iy] = (Njpsi_pp15_err[iarm][iy]/eff_pp_trk) / (NMB_pp[iarm]/eff_pp_bbc) / 0.25 / eff_y_pp[iarm][iy] / trigeff_y_pp[iarm][iy];
			Y_pp15_syserr[iarm][iy] = sys_corrbkg_pp[iarm]*sys_corrbkg_pp[iarm];
			Y_pp15_syserr[iarm][iy] += sys_var_pp[iarm]*sys_var_pp[iarm];
			Y_pp15_syserr[iarm][iy] += sys_acceff_pp[iarm]*sys_acceff_pp[iarm];
			Y_pp15_syserr[iarm][iy] += sys_trigeff_y_pp[iarm][iy]*sys_trigeff_y_pp[iarm][iy];
			Y_pp15_syserr[iarm][iy] += sys_pol_pp[iarm][iy]*sys_pol_pp[iarm][iy];
			Y_pp15_syserr[iarm][iy] = Y_pp15[iarm][iy] * sqrt(Y_pp15_syserr[iarm][iy]); 

			Xec_pp15[iarm][iy] = Y_pp15[iarm][iy] * 42.0 * 1e6;
			Xec_pp15_err[iarm][iy] = Y_pp15_err[iarm][iy] * 42.0 * 1e6;
			Xec_pp15_syserr[iarm][iy] = Y_pp15_syserr[iarm][iy] * 42.0 * 1e6;
		}
		Y_pp15_int[iarm] = (Njpsi_pp15_int[iarm]/eff_pp_trk) / (NMB_pp[iarm]/eff_pp_bbc) / 1.0 / eff_y_int_pp[iarm] / trigeff_y_int_pp[iarm];
		Y_pp15_int_err[iarm] = (Njpsi_pp15_int_err[iarm]/eff_pp_trk) / (NMB_pp[iarm]/eff_pp_bbc) / 1.0 / eff_y_int_pp[iarm] / trigeff_y_int_pp[iarm];
		Y_pp15_int_syserr[iarm] = sys_var_pp[iarm]*sys_var_pp[iarm];
		Y_pp15_int_syserr[iarm] += sys_acceff_pp[iarm]*sys_acceff_pp[iarm];
		Y_pp15_int_syserr[iarm] += sys_trigeff_y_int_pp[iarm]*sys_trigeff_y_int_pp[iarm];
		Y_pp15_int_syserr[iarm] += sys_pol_pp_int[iarm]*sys_pol_pp_int[iarm];
		Y_pp15_int_syserr[iarm] = Y_pp15_int[iarm] * sqrt(Y_pp15_int_syserr[iarm]); 

		Xec_pp15_int[iarm] = Y_pp15_int[iarm] * 42.0 * 1e6;
		Xec_pp15_int_err[iarm] = Y_pp15_int_err[iarm] * 42.0 * 1e6;
		Xec_pp15_int_syserr[iarm] = Y_pp15_int_syserr[iarm] * 42.0 * 1e6;
	}

	//return;

	TGraphErrors *gy_pp15[narm];
	TGraphErrors *gy_pp15_sys[narm];
	TGraphErrors *gy_xec_pp15[narm];
	TGraphErrors *gy_xec_pp15_sys[narm];

	TGraphErrors *gy_pp15_int[narm];
	TGraphErrors *gy_pp15_int_sys[narm];
	TGraphErrors *gy_xec_pp15_int[narm];
	TGraphErrors *gy_xec_pp15_int_sys[narm];

	for (int iarm=0; iarm<narm; iarm++){
		gy_pp15[iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Y_pp15[iarm][0], 0, &Y_pp15_err[iarm][0]);
		gy_pp15[iarm]->SetMarkerStyle(20);
		gy_pp15[iarm]->SetMarkerColor(1);
		gy_pp15[iarm]->SetLineWidth(1);
		gy_pp15[iarm]->SetLineColor(1);

		gy_pp15_sys[iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Y_pp15[iarm][0], &xx_err[iarm][0], &Y_pp15_syserr[iarm][0]);
		gy_pp15_sys[iarm]->SetMarkerStyle(20);
		gy_pp15_sys[iarm]->SetMarkerColor(1);
		gy_pp15_sys[iarm]->SetLineWidth(1);
		gy_pp15_sys[iarm]->SetLineColor(1);
		gy_pp15_sys[iarm]->SetFillStyle(0);

		gy_xec_pp15[iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Xec_pp15[iarm][0], 0, &Xec_pp15_err[iarm][0]);
		gy_xec_pp15[iarm]->SetMarkerStyle(20);
		gy_xec_pp15[iarm]->SetMarkerColor(1);
		gy_xec_pp15[iarm]->SetLineWidth(1);
		gy_xec_pp15[iarm]->SetLineColor(1);

		gy_xec_pp15_sys[iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Xec_pp15[iarm][0], &xx_err[iarm][0], &Xec_pp15_syserr[iarm][0]);
		gy_xec_pp15_sys[iarm]->SetMarkerStyle(29);
		gy_xec_pp15_sys[iarm]->SetMarkerColor(1);
		gy_xec_pp15_sys[iarm]->SetLineWidth(1);
		gy_xec_pp15_sys[iarm]->SetLineColor(1);
		gy_xec_pp15_sys[iarm]->SetFillStyle(0);

		gy_pp15_int[iarm] = new TGraphErrors(1, &xx_int[iarm], &Y_pp15_int[iarm], 0, &Y_pp15_int_err[iarm]);
		gy_pp15_int[iarm]->SetMarkerStyle(20);
		gy_pp15_int[iarm]->SetMarkerColor(1);
		gy_pp15_int[iarm]->SetLineWidth(1);
		gy_pp15_int[iarm]->SetLineColor(1);

		gy_pp15_int_sys[iarm] = new TGraphErrors(1, &xx_int[iarm], &Y_pp15_int[iarm], &xx_int_err[iarm], &Y_pp15_int_syserr[iarm]);
		gy_pp15_int_sys[iarm]->SetMarkerStyle(20);
		gy_pp15_int_sys[iarm]->SetMarkerColor(1);
		gy_pp15_int_sys[iarm]->SetLineWidth(1);
		gy_pp15_int_sys[iarm]->SetLineColor(1);
		gy_pp15_int_sys[iarm]->SetFillStyle(0);

		gy_xec_pp15_int[iarm] = new TGraphErrors(1, &xx_int[iarm], &Xec_pp15_int[iarm], 0, &Xec_pp15_int_err[iarm]);
		gy_xec_pp15_int[iarm]->SetMarkerStyle(20);
		gy_xec_pp15_int[iarm]->SetMarkerColor(2);
		gy_xec_pp15_int[iarm]->SetLineWidth(1);
		gy_xec_pp15_int[iarm]->SetLineColor(2);

		gy_xec_pp15_int_sys[iarm] = new TGraphErrors(1, &xx_int[iarm], &Xec_pp15_int[iarm], &xx_int_err[iarm], &Xec_pp15_int_syserr[iarm]);
		gy_xec_pp15_int_sys[iarm]->SetMarkerStyle(29);
		gy_xec_pp15_int_sys[iarm]->SetMarkerColor(2);
		gy_xec_pp15_int_sys[iarm]->SetLineWidth(1);
		gy_xec_pp15_int_sys[iarm]->SetLineColor(2);
		gy_xec_pp15_int_sys[iarm]->SetFillStyle(0);
	}

	TCanvas *cpp = new TCanvas("cpp","cpp",1.2*400,400);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3,0,3,60);
	SetHistoStyle("y","B_{ll}d#sigma/dy (nb)");

	for (int iarm=0; iarm<narm; iarm++){
		gy_xec_pp15_sys[iarm]->Draw("2");
		gy_xec_pp15[iarm]->Draw("p");

		gy_xec_pp15_int_sys[iarm]->SetMarkerColor(2);
		gy_xec_pp15_int_sys[iarm]->SetLineColor(2);
		gy_xec_pp15_int_sys[iarm]->Draw("2");
		gy_xec_pp15_int[iarm]->SetMarkerColor(2);
		gy_xec_pp15_int[iarm]->SetLineColor(2);
		gy_xec_pp15_int[iarm]->Draw("p");
	}

	//return;

	TFile *infile = new TFile("event_files/Run14HeAu200_eventhist_for_dimuon_20180314.root","read");

	TH1D *hBBCZ[ncent+1][narm];
	double NMB[ncent+1][narm];
	for (int icent=0; icent<ncent+1; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			hBBCZ[icent][iarm] = (TH1D*)infile->Get(Form("hBBCZ_CENT%02d%02d_arm%d",cent_perv2[icent],cent_perv2[icent+1],iarm));
			NMB[icent][iarm] = hBBCZ[icent][iarm]->Integral();
		}
	}

	//combine 40-60% and 60-88%
	for (int iarm=0; iarm<narm; iarm++){
		NMB[2][iarm] += NMB[3][iarm];
	}

	TH1D *hBBCZ_ALL[narm];
	TH1D *hBBCZ_A[narm];
	double NMB_ALL[narm], NMB_A[narm];
	for (int iarm=0; iarm<narm; iarm++){
		hBBCZ_ALL[iarm] = (TH1D*)infile->Get(Form("hBBCZ_arm%d",iarm));
		hBBCZ_A[iarm] = (TH1D*)infile->Get(Form("hBBCZ_Z_arm%d",iarm));
		NMB_ALL[iarm] = hBBCZ_ALL[iarm]->Integral();
		NMB_A[iarm] = hBBCZ_A[iarm]->Integral();
	}

	//return;

	TH1F *heff_y[narm];
	TH1F *heff_y_int[narm];
	TH1F *heff_y_A[narm];
	TH1F *heff_y_int_A[narm];

	double eff_y[ncent+1][narm][nybin], eff_y_A[ncent+1][narm][nybin];
	double eff_y_int[ncent+1][narm], eff_y_int_A[ncent+1][narm];

	for (int icent=0; icent<ncent+1; icent++){

		TFile *infile1 = new TFile(Form("eff_files/Run14HeAu200_CENT%02d%02d_acceff_InclusiveJpsi_B.root",cent_perv2[icent],cent_perv2[icent+1]),"read");
		TFile *infile1_A = new TFile(Form("eff_files/Run14HeAu200_CENT%02d%02d_acceff_InclusiveJpsi_A.root",cent_perv2[icent],cent_perv2[icent+1]),"read");

		for (int iarm=0; iarm<narm; iarm++){
			heff_y[iarm] = (TH1F*)infile1->Get(Form("heff_yy_w_dAu_set1_arm%d",iarm));
			heff_y_A[iarm] = (TH1F*)infile1_A->Get(Form("heff_yy_w_dAu_set1_arm%d",iarm));

			heff_y_int[iarm] = (TH1F*)infile1->Get(Form("heff_yy_int_w_dAu_set1_arm%d",iarm));
			heff_y_int_A[iarm] = (TH1F*)infile1_A->Get(Form("heff_yy_int_w_dAu_set1_arm%d",iarm));

			for (int iy=0; iy<nybin; iy++){
				eff_y[icent][iarm][iy] = heff_y[iarm]->GetBinContent(iy+1);
				eff_y_A[icent][iarm][iy] = heff_y_A[iarm]->GetBinContent(iy+1);

				cout << "HeAu CENT: " << icent<< ", ARM: " << iarm << ", Y: " << iy << ", eff0: " << eff_y[icent][iarm][iy] << ", eff1: " << eff_y_A[icent][iarm][iy] << endl;

				if ( iarm==0 ){
					double F = NMB_A[iarm]/NMB_ALL[iarm];
					eff_y[icent][iarm][iy] = eff_y[icent][iarm][iy]*(1-F) + eff_y_A[icent][iarm][iy]*F;
					cout << "HeAu CENT: " << icent << ", ARM: " << iarm << ", Y: " << iy << ", eff mod: " << eff_y[icent][iarm][iy] << endl;
				}
			}

			eff_y_int[icent][iarm] = heff_y_int[iarm]->GetBinContent(1);
			eff_y_int_A[icent][iarm] = heff_y_int_A[iarm]->GetBinContent(1);

			if ( iarm==0 ){
				double F = NMB_A[iarm]/NMB_ALL[iarm];
				eff_y_int[icent][iarm] = eff_y_int[icent][iarm]*(1-F) + eff_y_int_A[icent][iarm]*F;
			}

		}//iarm

		infile1->Close();

	}//icent

	for (int iarm=0; iarm<narm; iarm++){
		for (int iy=0; iy<nybin; iy++){
			eff_y[2][iarm][iy] = (20*eff_y[2][iarm][iy] + 28*eff_y[3][iarm][iy])/48.0;
		}//iy
		eff_y_int[2][iarm] = (20*eff_y_int[2][iarm] + 28*eff_y_int[3][iarm])/48.0;
	}//iarm

	//return;

	TH1F *htrigeff[narm];
	TH1F *htrigeff_int[narm];
	TH1F *htrigeff_A[narm];
	TH1F *htrigeff_int_A[narm];
	TH1F *htrigeff_wt[narm];
	TH1F *htrigeff_wt_int[narm];
	TH1F *htrigeff_wt_A[narm];
	TH1F *htrigeff_wt_int_A[narm];

	double trigeff_y[ncent+1][narm][nybin];
	double trigeff_y_int[ncent+1][narm];
	double trigeff_y_A[ncent+1][narm][nybin];
	double trigeff_y_int_A[ncent+1][narm];

	for (int icent=0; icent<ncent+1; icent++){

		TFile *infile2 = new TFile(Form("eff_files/Run14HeAu200_CENT%02d%02d_B_trigeff_MUID2D_output.root",cent_perv2[icent],cent_perv2[icent+1]),"read");
		TFile *infile2_A = new TFile(Form("eff_files/Run14HeAu200_CENT%02d%02d_A_trigeff_MUID2D_output.root",cent_perv2[icent],cent_perv2[icent+1]),"read");

		for (int iarm=0; iarm<narm; iarm++){

			htrigeff[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_set1_arm%d",iarm));
			htrigeff_A[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_half_set1_arm%d",iarm));
			htrigeff_int[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_int_set1_arm%d",iarm));
			htrigeff_int_A[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_int_half_set1_arm%d",iarm));

			htrigeff_wt[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_wt_set1_arm%d",iarm));
			htrigeff_wt_A[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_half_wt_set1_arm%d",iarm));
			htrigeff_wt_int[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_int_wt_set1_arm%d",iarm));
			htrigeff_wt_int_A[iarm] = (TH1F*)infile2->Get(Form("htrigeff_y_int_half_wt_set1_arm%d",iarm));

			for (int iy=0; iy<nybin; iy++){
				float eff0 = htrigeff[iarm]->GetBinContent(iy+1);
				float eff1 = htrigeff_wt[iarm]->GetBinContent(iy+1);
				trigeff_y[icent][iarm][iy] = (eff0 + eff1)/2;
				sys_trigeff_y[icent][iarm][iy] = fabs(eff0 - eff1)/(eff0 + eff1);

				eff0 = htrigeff_A[iarm]->GetBinContent(iy+1);
				eff1 = htrigeff_wt_A[iarm]->GetBinContent(iy+1);
				trigeff_y_A[icent][iarm][iy] = (eff0 + eff1)/2;

				cout << "HeAu CENT: " << icent<< ", ARM: " << iarm << ", Y: " << iy << ", TRIG EFF: " << trigeff_y[icent][iarm][iy] << ", TRIG EFF A: " << trigeff_y_A[icent][iarm][iy] << ", Err: " << sys_trigeff_y[icent][iarm][iy] << endl;

				if ( iarm==0 ){
					double F = NMB_A[iarm]/NMB_ALL[iarm];
					trigeff_y[icent][iarm][iy] = trigeff_y[icent][iarm][iy]*(1-F) + trigeff_y_A[icent][iarm][iy]*F;
					cout << "TRIG EFF MOD: " << trigeff_y[icent][iarm][iy] << endl; 
				}
			}

			float eff0 = htrigeff_int[iarm]->GetBinContent(1);
			float eff1 = htrigeff_wt_int[iarm]->GetBinContent(1);

			trigeff_y_int[icent][iarm] = (eff0 + eff1)/2;
			sys_trigeff_y_int[icent][iarm] = fabs(eff0 - eff1)/(eff0 + eff1);

		}//iarm

		infile2->Close();
	}//icent

	for (int iarm=0; iarm<narm; iarm++){
		for (int iy=0; iy<nybin; iy++){
			trigeff_y[2][iarm][iy] = (20*trigeff_y[2][iarm][iy] + 28*trigeff_y[3][iarm][iy])/48.0;
		}//iy
		trigeff_y_int[2][iarm] = (20*trigeff_y_int[2][iarm] + 28*trigeff_y_int[3][iarm])/48.0;
	}//iarm

	//return;

	double Y[ncent][narm][nybin], Y_err[ncent][narm][nybin], Y_syserr[ncent][narm][nybin];
	double rpA[ncent][narm][nybin], rpA_err[ncent][narm][nybin], rpA_syserr[ncent][narm][nybin];

	double Y_int[ncent][narm], Y_int_err[ncent][narm], Y_int_syserr[ncent][narm];
	double rpA_int[ncent][narm], rpA_int_err[ncent][narm], rpA_int_syserr[ncent][narm];

	TGraphErrors *gy[ncent][narm];
	TGraphErrors *gy_sys[ncent][narm];

	TGraphErrors *gR[ncent][narm];
	TGraphErrors *gR_sys[ncent][narm];

	TGraphErrors *gR_int[narm];
	TGraphErrors *gR_int_sys[narm];

	for (int iarm=0; iarm<narm; iarm++){
		gR_int[iarm] = new TGraphErrors;
		gR_int_sys[iarm] = new TGraphErrors;

		gR_int[iarm]->SetMarkerStyle(20+4*iarm);
		gR_int[iarm]->SetMarkerColor(iarm+1);
		gR_int[iarm]->SetLineColor(iarm+1);
		gR_int[iarm]->SetLineWidth(1);

		gR_int_sys[iarm]->SetFillStyle(0);
		gR_int_sys[iarm]->SetMarkerStyle(20+4*iarm);
		gR_int_sys[iarm]->SetMarkerColor(iarm+1);
		gR_int_sys[iarm]->SetLineColor(iarm+1);
		gR_int_sys[iarm]->SetLineWidth(1);
	}

	int np_cent[2] = {0};

	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			for (int iy=0; iy<nybin; iy++){

				cout << icent << " " << iarm << " " << iy << " " << Ncoll_cent[icent] << " " << eff_y[icent][iarm][iy] << " " << trigeff_y[icent][iarm][iy] << endl;

				Y[icent][iarm][iy] = Njpsi[icent][iarm][iy] / NMB[icent][iarm] / 0.25 / Ncoll_cent[icent] * BiasF_cent[icent] / eff_y[icent][iarm][iy] / trigeff_y[icent][iarm][iy];
				Y_err[icent][iarm][iy] = Njpsi_err[icent][iarm][iy] / NMB[icent][iarm] / 0.25 / Ncoll_cent[icent] * BiasF_cent[icent] / eff_y[icent][iarm][iy] / trigeff_y[icent][iarm][iy];
				Y_syserr[icent][iarm][iy] = sys_corrbkg[icent][iarm]*sys_corrbkg[icent][iarm];
				Y_syserr[icent][iarm][iy] += sys_var[iarm]*sys_var[iarm];
				Y_syserr[icent][iarm][iy] += sys_acceff[iarm]*sys_acceff[iarm];
				Y_syserr[icent][iarm][iy] += sys_trigeff_y[icent][iarm][iy]*sys_trigeff_y[icent][iarm][iy];
				Y_syserr[icent][iarm][iy] += sys_pol_pp[iarm][iy]*sys_pol_pp[iarm][iy];
				Y_syserr[icent][iarm][iy] = Y[icent][iarm][iy] * sqrt(Y_syserr[icent][iarm][iy]); 

				rpA[icent][iarm][iy] = Y[icent][iarm][iy] / Y_pp15[iarm][iy];
				float err0 = Y_pp15_err[iarm][iy]/Y_pp15[iarm][iy];
				float err1 = Y_err[icent][iarm][iy]/Y[icent][iarm][iy];

				//TypeA err
				rpA_err[icent][iarm][iy] = rpA[icent][iarm][iy]*sqrt(err0*err0 + err1*err1 + sys_corrbkg_pp[iarm]*sys_corrbkg_pp[iarm] + sys_corrbkg[icent][iarm]*sys_corrbkg[icent][iarm]); 

				//TypeB err
				err0 = sqrt(sys_var_pp[iarm]*sys_var_pp[iarm] + sys_acceff_pp[iarm]*sys_acceff_pp[iarm] + sys_trigeff_y_pp[iarm][iy]*sys_trigeff_y_pp[iarm][iy]);
				err1 = sqrt(sys_var[iarm]*sys_var[iarm] + sys_acceff[iarm]*sys_acceff[iarm] + sys_trigeff_y[icent][iarm][iy]*sys_trigeff_y[icent][iarm][iy]);
				rpA_syserr[icent][iarm][iy] = rpA[icent][iarm][iy]*sqrt(err0*err0 + err1*err1);

			}//iy

			Y_int[icent][iarm] = Njpsi_cent[icent][iarm] / NMB[icent][iarm] / 1.0 / Ncoll_cent[icent] * BiasF_cent[icent] / eff_y_int[icent][iarm] / trigeff_y_int[icent][iarm];
			Y_int_err[icent][iarm] = Njpsi_cent_err[icent][iarm] / NMB[icent][iarm] / 1.0 / Ncoll_cent[icent] * BiasF_cent[icent] / eff_y_int[icent][iarm] / trigeff_y_int[icent][iarm];
			Y_int_syserr[icent][iarm] = sys_corrbkg[icent][iarm]*sys_corrbkg[icent][iarm];
			Y_int_syserr[icent][iarm] += sys_var[iarm]*sys_var[iarm];
			Y_int_syserr[icent][iarm] += sys_acceff[iarm]*sys_acceff[iarm];
			Y_int_syserr[icent][iarm] += sys_trigeff_y_int[icent][iarm]*sys_trigeff_y_int[icent][iarm];
			Y_int_syserr[icent][iarm] = Y_int[icent][iarm] * sqrt(Y_int_syserr[icent][iarm]); 

			rpA_int[icent][iarm] = Y_int[icent][iarm] / Y_pp15_int[iarm];
			float err0 = Y_pp15_int_err[iarm]/Y_pp15_int[iarm];
			float err1 = Y_int_err[icent][iarm]/Y_int[icent][iarm];

			//TypeA err
			rpA_int_err[icent][iarm] = rpA_int[icent][iarm]*sqrt(err0*err0 + err1*err1 + sys_corrbkg_pp[iarm]*sys_corrbkg_pp[iarm] + sys_corrbkg[icent][iarm]*sys_corrbkg[icent][iarm]); 

			//TypeB err
			err0 = sqrt(sys_var_pp[iarm]*sys_var_pp[iarm] + sys_acceff_pp[iarm]*sys_acceff_pp[iarm] + sys_trigeff_y_int_pp[iarm]*sys_trigeff_y_int_pp[iarm]);
			err1 = sqrt(sys_var[iarm]*sys_var[iarm] + sys_acceff[iarm]*sys_acceff[iarm] + sys_trigeff_y_int[icent][iarm]*sys_trigeff_y_int[icent][iarm]);
			rpA_int_syserr[icent][iarm] = rpA_int[icent][iarm]*sqrt(err0*err0 + err1*err1);

			gR_int[iarm]->SetPoint(np_cent[iarm], Ncoll_cent[icent],  rpA_int[icent][iarm]);
			gR_int[iarm]->SetPointError(np_cent[iarm], 0, rpA_int_err[icent][iarm]);
			gR_int_sys[iarm]->SetPoint(np_cent[iarm], Ncoll_cent[icent],  rpA_int[icent][iarm]);
			gR_int_sys[iarm]->SetPointError(np_cent[iarm], Ncoll_syserr, rpA_int_syserr[icent][iarm]);
			np_cent[iarm]++;

			gy[icent][iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Y[icent][iarm][0], 0, &Y_err[icent][iarm][0]);
			gy[icent][iarm]->SetMarkerStyle(24);
			gy[icent][iarm]->SetMarkerColor(nColor[icent]);
			gy[icent][iarm]->SetLineWidth(1);
			gy[icent][iarm]->SetLineColor(nColor[icent]);

			gy_sys[icent][iarm] = new TGraphErrors(nybin, &xx[iarm][0], &Y[icent][iarm][0], &xx_err[iarm][0], &Y_syserr[icent][iarm][0]);
			gy_sys[icent][iarm]->SetMarkerStyle(24);
			gy_sys[icent][iarm]->SetMarkerColor(nColor[icent]);
			gy_sys[icent][iarm]->SetLineWidth(1);
			gy_sys[icent][iarm]->SetLineColor(nColor[icent]);
			gy_sys[icent][iarm]->SetFillStyle(0);

			gR[icent][iarm] = new TGraphErrors(nybin, &xx[iarm][0], &rpA[icent][iarm][0], 0, &rpA_err[icent][iarm][0]);
			gR[icent][iarm]->SetMarkerStyle(24);
			gR[icent][iarm]->SetMarkerColor(nColor[icent]);
			gR[icent][iarm]->SetLineWidth(1);
			gR[icent][iarm]->SetLineColor(nColor[icent]);

			gR_sys[icent][iarm] = new TGraphErrors(nybin, &xx[iarm][0], &rpA[icent][iarm][0], &xx_err[iarm][0], &rpA_syserr[icent][iarm][0]);
			gR_sys[icent][iarm]->SetMarkerStyle(24);
			gR_sys[icent][iarm]->SetMarkerColor(nColor[icent]);
			gR_sys[icent][iarm]->SetLineWidth(1);
			gR_sys[icent][iarm]->SetLineColor(nColor[icent]);
			gR_sys[icent][iarm]->SetFillStyle(0);
		}//iarm
	}//icent

	//return;

	TCanvas *c10 = new TCanvas("c10","c10",1.3*3*300,1*300);
	c10->Divide(3,1);

	TCanvas *c12 = new TCanvas("c12","c12",1.3*3*300,1*300);
	c12->Divide(3,1);

	for (int icent=0; icent<ncent; icent++){

		c10->cd(icent+1);
		SetPadStyle();
		gPad->SetTopMargin(0.08);
		gPad->SetRightMargin(0.015);
		gPad->SetLeftMargin(0.14);
		gPad->SetBottomMargin(0.14);
		htmp = (TH1F*)gPad->DrawFrame(-3,0,3,2.0e-6);
		SetHistoStyle("y","B_{ll}dN/dy","",22,18);
		htmp->GetYaxis()->SetTitleOffset(0.8);
		htmp->GetXaxis()->SetTitleOffset(0.8);
		for (int iarm=0; iarm<narm; iarm++){
			gy_pp15_sys[iarm]->Draw("2");
			gy_pp15[iarm]->Draw("p");

			gy_sys[icent][iarm]->Draw("2");
			gy[icent][iarm]->Draw("p");
		}

		{
			TLegend *leg = new TLegend(0.20,0.68,0.9,0.88);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","^{3}He+Au #sqrt{s_{NN}}=200 GeV","h");
			le = leg->AddEntry("","Inclusive J/#psi","h");
			le = leg->AddEntry("",Form("%d%c-%d%c",cent_per[icent],'%',cent_per[icent+1],'%'),"h");
			leg->Draw();
		}

		{
			TLegend *leg = new TLegend(0.6,0.68,0.93,0.88);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","","");
			le = leg->AddEntry(gy[icent][0],"^{3}He+Au/#LTN_{coll}#GT","P");
			le = leg->AddEntry(gy_pp15[0],"p+p","P");
			leg->Draw();
		}

		c12->cd(icent+1);
		SetPadStyle();
		gPad->SetRightMargin(0.015);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.14);
		htmp = (TH1F*)gPad->DrawFrame(-3,0,3,2.0);
		SetHistoStyle("y","R_{AB}","",22,18);
		htmp->GetYaxis()->SetTitleOffset(0.9);
		htmp->GetXaxis()->SetTitleOffset(0.8);

		{
			TLegend *leg = new TLegend(0.20,0.7,0.9,0.90);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","^{3}He+Au #sqrt{s_{NN}}=200 GeV","h");
			le = leg->AddEntry("","Inclusive J/#psi","h");
			le = leg->AddEntry("",Form("%d%c-%d%c",cent_per[icent],'%',cent_per[icent+1],'%'),"h");
			leg->Draw();

			TLine *line = new TLine(-3, 1, 3, 1);
			line->SetLineStyle(2);
			line->Draw();

			float global_sys = sqrt(Ncoll_cent_sys[icent]*Ncoll_cent_sys[icent] + BiasF_cent_sys[icent]*BiasF_cent_sys[icent] + 0.101*0.101);
			TBox *box = new TBox(2.8,1-global_sys,3.0,1+global_sys);
			box->SetFillStyle(1000);
			box->SetFillColorAlpha(nColor[icent],0.5);
			box->Draw();
		}


		for (int iarm=0; iarm<narm; iarm++){
			gR_sys[icent][iarm]->Draw("2");
			gR[icent][iarm]->Draw("p");
		}

	}

	TCanvas *c13 = new TCanvas("c13","c13",1.3*400,400);
	SetPadStyle();
	gPad->SetRightMargin(0.03);
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.16);
	htmp = (TH1F*)gPad->DrawFrame(0,0,25,2.0);
	SetHistoStyle("#LTN_{coll}#GT","R_{AB}","",24,20);
	htmp->GetYaxis()->SetTitleOffset(1.0);
	htmp->GetXaxis()->SetTitleOffset(1.0);

	{
		TLegend *leg = new TLegend(0.40,0.65,0.9,0.90);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		le = leg->AddEntry("","^{3}He+Au #sqrt{s_{NN}}=200 GeV","");
		le = leg->AddEntry("","J/#psi#rightarrow#mu#mu","");
		le = leg->AddEntry(gR_int[0],"-2.2<y<-1.2","P");
		le = leg->AddEntry(gR_int[1],"1.2<y<2.2","P");
		leg->Draw();

		TLine *line = new TLine(0, 1, 25, 1);
		line->SetLineStyle(2);
		line->Draw();

		TBox *box = new TBox(25*0.97,1-0.101,25,1+0.101);
		box->SetFillStyle(1000);
		box->SetFillColorAlpha(1,0.5);
		box->Draw();
	}

	for (int iarm=0; iarm<narm; iarm++){
		gR_int_sys[iarm]->Draw("2");
		gR_int[iarm]->Draw("p");

	}

	if ( bWRITE ){
		TFile *outfile = new TFile("RAB_y_Ncoll_HeAu200.root","recreate");

		for (int iarm=0; iarm<narm; iarm++){
			gR_int_sys[iarm]->Write(Form("RAB_Ncoll_syserr_arm%d",iarm));
			gR_int[iarm]->Write(Form("RAB_Ncoll_staerr_arm%d",iarm));
		}//iarm

		for (int icent=0; icent<ncent; icent++){
			for (int iarm=0; iarm<narm; iarm++){
				gR_sys[icent][iarm]->Write(Form("RAB_y_syserr_cent%d_arm%d",icent,iarm));
				gR[icent][iarm]->Write(Form("RAB_y_staerr_cent%d_arm%d",icent,iarm));
			}
		}//icent

	}//bWRITE

}
