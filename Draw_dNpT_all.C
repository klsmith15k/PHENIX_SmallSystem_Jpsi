#include "Style.h"

void Draw_dNpT_all()
{
	const bool bSAVE = false;

	SetStyle();

  char name1[500];
  char name2[500];
  char name3[500];

	TFile *file_pAu;	  
	TFile *file_pAl;	  
	TFile *file_HeAu;	  

	const int narm = 2;
	const int nptbin_pAu = 19;
	const int nptbin_pAl = 18;
	const int nptbin_HeAu = 15;

	const int ncentbin_pAu = 8;
	const int ncentbin_pAl = 4;
	const int ncentbin_HeAu = 4;

	const int cent_per_pAu_min[ncentbin_pAu] = {0, 5, 10, 0, 20, 40, 60, 0};
	const int cent_per_pAu_max[ncentbin_pAu] = {5, 10, 20, 20, 40, 60, 84, 100};
	const float scaler_pAu[ncentbin_pAu] = {1, 1e-1, 1e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
	const int nMarker_pAu[ncentbin_pAu] = {20, 24, 21, 32, 32, 22, 26, 25};
	const int nColor_pAu[ncentbin_pAu] = {2, 4, 6, 6, 1, kGreen+2, kOrange+6, 1};

	const int cent_per_pAl_min[ncentbin_pAl] = {0, 20, 40, 0};
	const int cent_per_pAl_max[ncentbin_pAl] = {20, 40, 72, 100};
	const float scaler_pAl[ncentbin_pAl] = {1, 1e-1, 1e-2, 1e-3};
	const int nMarker_pAl[ncentbin_pAl] = {20, 24, 21, 25};
	const int nColor_pAl[ncentbin_pAl] = {2, 4, 6, 1};

	const int cent_per_HeAu_min[ncentbin_HeAu] = {0, 20, 40, 0};
	const int cent_per_HeAu_max[ncentbin_HeAu] = {20, 40, 88, 100};
	const float scaler_HeAu[ncentbin_HeAu] = {1, 1e-1, 1e-2, 1e-3};
	const int nMarker_HeAu[ncentbin_HeAu] = {20, 24, 21, 25};
	const int nColor_HeAu[ncentbin_HeAu] = {2, 4, 6, 1};

	const char c_arm[narm] = {'S', 'N'};

	const double Ncoll_pAu_sys[ncentbin_pAu] = {0.6/9.7, 0.5/8.4, 0.5/7.4, 0.5/8.2, 0.4/6.1, 0.3/4.4, 0.2/2.6};
	const double BiasF_pAu_sys[ncentbin_pAu] = {0.01/0.86, 0.01/0.90, 0.01/0.94, 0.01/0.90, 0.01/0.98, 0.01/1.02, 0.06/1.00};

	const double Ncoll_pAl_sys[ncentbin_pAl] = {0.3/3.5, 0.1/2.3, 0.1/1.7};
	const double BiasF_pAl_sys[ncentbin_pAl] = {0.01/0.81, 0.02/0.90, 0.04/1.01};

	const double Ncoll_HeAu_sys[ncentbin_HeAu] = {0.076233184, 0.074324324, 0.080272996};
	const double BiasF_HeAu_sys[ncentbin_HeAu] = {0.010526316, 0.0099009901, 0.035018614};

	// RAB arrays for TGraphs
	double YpAu_array[ncentbin_pAu][narm][nptbin_pAu] = {0};
	double YpAl_array[ncentbin_pAl][narm][nptbin_pAl] = {0};
	double YHeAu_array[ncentbin_HeAu][narm][nptbin_HeAu] = {0};

  // Type A error arrays for TGraphs
  double x_errors[nptbin_pAu] = {0};
  // Width of type B Errors
  double width[nptbin_pAu] = {0};

	double YpAu_ERR_array[ncentbin_pAu][narm][nptbin_pAu] = {0};
	double YpAl_ERR_array[ncentbin_pAl][narm][nptbin_pAl] = {0};
	double YHeAu_ERR_array[ncentbin_HeAu][narm][nptbin_HeAu] = {0};

	// Type B absolute error arrays for TGraphs
	double YpAu_sys_array[ncentbin_pAu][narm][nptbin_pAu] = {0};
	double YpAl_sys_array[ncentbin_pAl][narm][nptbin_pAl] = {0};
	double YHeAu_sys_array[ncentbin_HeAu][narm][nptbin_HeAu] = {0};

	TH1F *YpAu[ncentbin_pAu][narm];
	TH1F *YpAl[ncentbin_pAl][narm];
	TH1F *YHeAu[ncentbin_HeAu][narm];

	TH1F *YpAu_SYS_FRAC_ERR[ncentbin_pAu][narm];
	TH1F *YpAl_SYS_FRAC_ERR[ncentbin_pAl][narm];
	TH1F *YHeAu_SYS_FRAC_ERR[ncentbin_HeAu][narm];

	file_pAu = new TFile("RpAu_YpAu_vs_pT_cent.root","read");
	file_pAl = new TFile("RpAl_YpAl_vs_pT_cent.root","read");
	file_HeAu = new TFile("RHeAu_YHeAu_vs_pT_cent.root","read");

	TFile *file_pAl_mb = new TFile("RpAl_YpAl_vs_pT_0100.root","read");
	TFile *file_pAu_mb = new TFile("RpAu_YpAu_vs_pT_0100.root","read");
	TFile *file_HeAu_mb = new TFile("RHeAu_YHeAu_vs_pT_0100.root","read");

	for (int iarm=0; iarm<narm; iarm++){
		for (int icent=0; icent<ncentbin_pAu-1; icent++){
			YpAu[icent][iarm] = (TH1F*)file_pAu->Get(Form("YpAu_vs_pT_%c_%d%d",c_arm[iarm],cent_per_pAu_min[icent],cent_per_pAu_max[icent]));
			YpAu_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_pAu->Get(Form("YpAu_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_pAu_min[icent],cent_per_pAu_max[icent]));
		}
		YpAu[ncentbin_pAu-1][iarm] = (TH1F*)file_pAu_mb->Get(Form("YpAu_vs_pT_%c_0100",c_arm[iarm]));
		YpAu_SYS_FRAC_ERR[ncentbin_pAu-1][iarm] = (TH1F*)file_pAu_mb->Get(Form("YpAu_vs_pT_%c_0100_sys_frac",c_arm[iarm]));

		for (int icent=0; icent<ncentbin_pAl-1; icent++){
			YpAl[icent][iarm] = (TH1F*)file_pAl->Get(Form("YpAl_vs_pT_%c_%d%d",c_arm[iarm],cent_per_pAl_min[icent],cent_per_pAl_max[icent]));
			YpAl_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_pAl->Get(Form("YpAl_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_pAl_min[icent],cent_per_pAl_max[icent]));
		}
		YpAl[ncentbin_pAl-1][iarm] = (TH1F*)file_pAl_mb->Get(Form("YpAl_vs_pT_%c_0100",c_arm[iarm]));
		YpAl_SYS_FRAC_ERR[ncentbin_pAl-1][iarm] = (TH1F*)file_pAl_mb->Get(Form("YpAl_vs_pT_%c_0100_sys_frac",c_arm[iarm]));

		for (int icent=0; icent<ncentbin_HeAu-1; icent++){
			YHeAu[icent][iarm] = (TH1F*)file_HeAu->Get(Form("YHeAu_vs_pT_%c_%d%d",c_arm[iarm],cent_per_HeAu_min[icent],cent_per_HeAu_max[icent]));
			YHeAu_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_HeAu->Get(Form("YHeAu_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_HeAu_min[icent],cent_per_HeAu_max[icent]));
		}
		YHeAu[ncentbin_HeAu-1][iarm] = (TH1F*)file_HeAu_mb->Get(Form("YHeAu_vs_pT_%c_0100",c_arm[iarm]));
		YHeAu_SYS_FRAC_ERR[ncentbin_HeAu-1][iarm] = (TH1F*)file_HeAu_mb->Get(Form("YHeAu_vs_pT_%c_0100_sys_frac",c_arm[iarm]));
	}

	//return;

	TGraphErrors *gr_pAu[ncentbin_pAu][narm];
	TGraphErrors *gr_pAl[ncentbin_pAl][narm];
	TGraphErrors *gr_HeAu[ncentbin_HeAu][narm];

	TGraphErrors *gr_pAu_sys[ncentbin_pAu][narm];
	TGraphErrors *gr_pAl_sys[ncentbin_pAl][narm];
	TGraphErrors *gr_HeAu_sys[ncentbin_HeAu][narm];

	//plot cent arrays
	double pt_array_pAu[nptbin_pAu] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0};
	double pt_array_pAl[nptbin_pAl] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.5,6.0};
	double pt_array_HeAu[nptbin_HeAu] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,4.5,6.0};

	const static int nptbin_mb_pAu = 26;
	const static int nptbin_mb_pAl = 20;
	const static int nptbin_mb_HeAu = 19;

	//plot MB arrays
	double pt_array_mb_pAu[nptbin_mb_pAu] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.25,6.75};
	double pt_array_mb_pAl[nptbin_mb_pAl] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,5.5,6.5};
	double pt_array_mb_HeAu[nptbin_mb_HeAu] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0};

	for (int icent=0; icent<ncentbin_pAu; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			gr_pAu[icent][iarm] = new TGraphErrors;
			gr_pAu[icent][iarm]->SetMarkerStyle(nMarker_pAu[icent]);
			gr_pAu[icent][iarm]->SetMarkerColor(nColor_pAu[icent]);
			gr_pAu[icent][iarm]->SetMarkerSize(0.8);
			gr_pAu[icent][iarm]->SetLineColor(nColor_pAu[icent]);

			gr_pAu_sys[icent][iarm] = new TGraphErrors;
			gr_pAu_sys[icent][iarm]->SetMarkerStyle(nMarker_pAu[icent]);
			gr_pAu_sys[icent][iarm]->SetMarkerColor(nColor_pAu[icent]);
			gr_pAu_sys[icent][iarm]->SetMarkerSize(0.8);
			gr_pAu_sys[icent][iarm]->SetLineColor(nColor_pAu[icent]);
			gr_pAu_sys[icent][iarm]->SetFillStyle(0);

			for (int ipt=0; ipt<nptbin_mb_pAu; ipt++){
			 double xx;
			  double yy;
			  double yy_err;
			  double yy_syserr;
				if(ipt < nptbin_pAu)
				  {
				    xx = pt_array_pAu[ipt];
				    yy = YpAu[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YpAu[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YpAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_pAu[icent][iarm]->SetPoint(ipt, pt_array_pAu[ipt], yy*scaler_pAu[icent]);
				    gr_pAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAu[icent]);
				    gr_pAu_sys[icent][iarm]->SetPoint(ipt, pt_array_pAu[ipt], yy*scaler_pAu[icent]);
				    gr_pAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAu[icent]);
				  }
				if(icent == ncentbin_pAu-1)
				  {
				    xx = pt_array_mb_pAu[ipt];
				    yy = YpAu[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YpAu[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YpAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_pAu[icent][iarm]->SetPoint(ipt, pt_array_mb_pAu[ipt], yy*scaler_pAu[icent]);
				    gr_pAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAu[icent]);
				    gr_pAu_sys[icent][iarm]->SetPoint(ipt, pt_array_mb_pAu[ipt], yy*scaler_pAu[icent]);
				    gr_pAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAu[icent]);
				  }
			}//ipt
		}//iarm
	}//icent

	//return;

	for (int icent=0; icent<ncentbin_pAl; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			gr_pAl[icent][iarm] = new TGraphErrors;
			gr_pAl[icent][iarm]->SetMarkerStyle(nMarker_pAl[icent]);
			gr_pAl[icent][iarm]->SetMarkerColor(nColor_pAl[icent]);
			gr_pAl[icent][iarm]->SetMarkerSize(0.8);
			gr_pAl[icent][iarm]->SetLineColor(nColor_pAl[icent]);

			gr_pAl_sys[icent][iarm] = new TGraphErrors;
			gr_pAl_sys[icent][iarm]->SetMarkerStyle(nMarker_pAl[icent]);
			gr_pAl_sys[icent][iarm]->SetMarkerColor(nColor_pAl[icent]);
			gr_pAl_sys[icent][iarm]->SetMarkerSize(0.8);
			gr_pAl_sys[icent][iarm]->SetLineColor(nColor_pAl[icent]);
			gr_pAl_sys[icent][iarm]->SetFillStyle(0);

			for (int ipt=0; ipt<nptbin_mb_pAl; ipt++){
			  double xx;
			  double yy;
			  double yy_err;
			  double yy_syserr;
				if(ipt <nptbin_pAl)
				  {
				    xx = pt_array_pAl[ipt];
				    yy = YpAl[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YpAl[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YpAl_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_pAl[icent][iarm]->SetPoint(ipt, pt_array_pAl[ipt], yy*scaler_pAl[icent]);
				    gr_pAl[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAl[icent]);
				    gr_pAl_sys[icent][iarm]->SetPoint(ipt, pt_array_pAl[ipt], yy*scaler_pAl[icent]);
				    gr_pAl_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAl[icent]);
				  }
				if(icent == ncentbin_pAl-1)
				  {
				    xx = pt_array_mb_pAl[ipt];
				    yy = YpAl[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YpAl[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YpAl_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_pAl[icent][iarm]->SetPoint(ipt, pt_array_mb_pAl[ipt], yy*scaler_pAl[icent]);
				    gr_pAl[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAl[icent]);
				    gr_pAl_sys[icent][iarm]->SetPoint(ipt, pt_array_mb_pAl[ipt], yy*scaler_pAl[icent]);
				    gr_pAl_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAl[icent]);
				  }
			}//ipt
		}//iarm
	}//icent

	for (int icent=0; icent<ncentbin_HeAu; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			gr_HeAu[icent][iarm] = new TGraphErrors;
			gr_HeAu[icent][iarm]->SetMarkerStyle(nMarker_HeAu[icent]);
			gr_HeAu[icent][iarm]->SetMarkerColor(nColor_HeAu[icent]);
			gr_HeAu[icent][iarm]->SetMarkerSize(0.8);
			gr_HeAu[icent][iarm]->SetLineColor(nColor_HeAu[icent]);

			gr_HeAu_sys[icent][iarm] = new TGraphErrors;
			gr_HeAu_sys[icent][iarm]->SetMarkerStyle(nMarker_HeAu[icent]);
			gr_HeAu_sys[icent][iarm]->SetMarkerColor(nColor_HeAu[icent]);
			gr_HeAu_sys[icent][iarm]->SetMarkerSize(0.8);
			gr_HeAu_sys[icent][iarm]->SetLineColor(nColor_HeAu[icent]);
			gr_HeAu_sys[icent][iarm]->SetFillStyle(0);

			for (int ipt=0; ipt<nptbin_mb_HeAu; ipt++){
			  double xx;
			  double yy;
			  double yy_err;
			  double yy_syserr;
				if(ipt < nptbin_HeAu)
				  {
				    xx = pt_array_HeAu[ipt];
				    yy = YHeAu[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YHeAu[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YHeAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_HeAu[icent][iarm]->SetPoint(ipt, pt_array_HeAu[ipt], yy*scaler_HeAu[icent]);
				    gr_HeAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_HeAu[icent]);
				    gr_HeAu_sys[icent][iarm]->SetPoint(ipt, pt_array_HeAu[ipt], yy*scaler_HeAu[icent]);
				    gr_HeAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_HeAu[icent]);
				  }
				if(icent == ncentbin_HeAu-1)
				  {
				    xx = pt_array_mb_HeAu[ipt];
				    yy = YHeAu[icent][iarm]->GetBinContent(ipt+1);
				    yy_err = YHeAu[icent][iarm]->GetBinError(ipt+1);
				    yy_syserr = yy*YHeAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				    gr_HeAu[icent][iarm]->SetPoint(ipt, pt_array_mb_HeAu[ipt], yy*scaler_HeAu[icent]);
				    gr_HeAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_HeAu[icent]);
				    gr_HeAu_sys[icent][iarm]->SetPoint(ipt, pt_array_mb_HeAu[ipt], yy*scaler_HeAu[icent]);
				    gr_HeAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_HeAu[icent]);
				  }
				
			}//ipt
		}//iarm
	}//icent

	TCanvas *c00[narm];;

	for (int iarm=0; iarm<narm; iarm++){
		c00[iarm] = new TCanvas(Form("c00_arm%d",iarm),Form("c00_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetTopMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.15);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-21,7.8,5e-5);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.2);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.18,0.17,0.93,0.34);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetNColumns(3);
			for (int icent=0; icent<ncentbin_pAu; icent++){
				if ( icent==3 ){
					continue;
				}
				gr_pAu_sys[icent][iarm]->Draw("2");
				le = leg->AddEntry(gr_pAu_sys[icent][iarm],Form("%d%c-%d%c#times10^{%d}",cent_per_pAu_min[icent],'%',cent_per_pAu_max[icent],'%',int(TMath::Log10(scaler_pAu[icent])-0.00001)),"P");
				le->SetTextSize(17);
				leg->Draw();
			}//
			leg->Draw();
		}

		{
			for (int icent=0; icent<ncentbin_pAu; icent++){
				if ( icent==3 ){
					continue;
				}
				gr_pAu[icent][iarm]->Draw("P");
			}//
		}

		{
			TLegend *leg = new TLegend(0.6,0.78,0.9,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","p+Au #sqrt{s_{NN}}=200 GeV","h");
			le->SetTextSize(20);
			le = leg->AddEntry("","Inclusive J/#psi","h");
			le->SetTextSize(20);
			if ( iarm==0 ){
				le = leg->AddEntry("","-2.2<y<-1.2","h");
				le->SetTextSize(20);
			}else{
				le = leg->AddEntry("","1.2<y<2.2","h");
				le->SetTextSize(20);
			}
			leg->Draw();
		}

		{
			TLatex *tex = new TLatex(0.5,3e-6,"PHENIX");
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}

		{
			TLatex *tex = new TLatex(7.,4e-8,Form("(%c)",97+iarm));
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}
	}

	//return;

	TCanvas *c01[narm];

	for (int iarm=0; iarm<narm; iarm++){
		c01[iarm] = new TCanvas(Form("c01_arm%d",iarm),Form("c01_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetTopMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.15);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-16,7.8,5e-6);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.2);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.18,0.20,0.73,0.32);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetNColumns(2);
			for (int icent=0; icent<ncentbin_pAl; icent++){
				gr_pAl_sys[icent][iarm]->Draw("2");
				le = leg->AddEntry(gr_pAl_sys[icent][iarm],Form("%d%c-%d%c#times10^{%d}",cent_per_pAl_min[icent],'%',cent_per_pAl_max[icent],'%',int(TMath::Log10(scaler_pAl[icent])-0.00001)),"P");
				leg->Draw();
			}//
			leg->Draw();
		}

		{
			for (int icent=0; icent<ncentbin_pAl; icent++){
				gr_pAl[icent][iarm]->Draw("P");
			}//
		}

		{
			TLegend *leg = new TLegend(0.6,0.78,0.9,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","p+Al #sqrt{s_{NN}}=200 GeV","h");
			le->SetTextSize(20);
			le = leg->AddEntry("","Inclusive J/#psi","h");
			le->SetTextSize(20);
			if ( iarm==0 ){
				le = leg->AddEntry("","-2.2<y<-1.2","h");
				le->SetTextSize(20);
			}else{
				le = leg->AddEntry("","1.2<y<2.2","h");
				le->SetTextSize(20);
			}
			leg->Draw();
		}

		{
			TLatex *tex = new TLatex(0.3,1e-6,"PHENIX");
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}

		{
			TLatex *tex = new TLatex(7.,9e-8,Form("(%c)",97+iarm));
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}
	}

	//return;


	TCanvas *c02[narm];;

	for (int iarm=0; iarm<narm; iarm++){
		c02[iarm] = new TCanvas(Form("c02_arm%d",iarm),Form("c02_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetTopMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.15);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-15,7.8,5e-5);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.2);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.18,0.20,0.73,0.32);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetNColumns(2);
			for (int icent=0; icent<ncentbin_HeAu; icent++){
				gr_HeAu_sys[icent][iarm]->Draw("2");
				leg->AddEntry(gr_HeAu_sys[icent][iarm],Form("%d%c-%d%c#times10^{%d}",cent_per_HeAu_min[icent],'%',cent_per_HeAu_max[icent],'%',int(TMath::Log10(scaler_HeAu[icent])-0.00001)),"P");
				leg->Draw();
			}//
			leg->Draw();
		}

		{
			for (int icent=0; icent<ncentbin_HeAu; icent++){
				gr_HeAu[icent][iarm]->Draw("P");
			}//
		}

		{
			TLegend *leg = new TLegend(0.55,0.78,0.9,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","^{3}He+Au #sqrt{s_{NN}}=200 GeV","h");
			le->SetTextSize(20);
			le = leg->AddEntry("","Inclusive J/#psi","h");
			le->SetTextSize(20);
			if ( iarm==0 ){
				le = leg->AddEntry("","-2.2<y<-1.2","h");
				le->SetTextSize(20);
			}else{
				le = leg->AddEntry("","1.2<y<2.2","h");
				le->SetTextSize(20);
			}
			leg->Draw();
		}

		{
			TLatex *tex = new TLatex(0.3,1e-5,"PHENIX");
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}

		{
			TLatex *tex = new TLatex(7.,9e-7,Form("(%c)",97+iarm));
			tex->SetTextFont(43);
			tex->SetTextSize(20);
			tex->Draw();
		}
	}

	if ( bSAVE ){

		for (int iarm=0; iarm<narm; iarm++){

			c00[iarm]->cd();
			c00[iarm]->SaveAs(Form("pdf/fig_dNdpT_pAu_arm%d.pdf",iarm));

			c01[iarm]->cd();
			c01[iarm]->SaveAs(Form("pdf/fig_dNdpT_pAl_arm%d.pdf",iarm));

			c02[iarm]->cd();
			c02[iarm]->SaveAs(Form("pdf/fig_dNdpT_HeAu_arm%d.pdf",iarm));

		}//iarm

	}//bSAVE

}


