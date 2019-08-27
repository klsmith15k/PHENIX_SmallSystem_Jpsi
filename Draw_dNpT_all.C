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
	const int nptbin_pAl = 16;
	const int nptbin_HeAu = 13;

	const int ncentbin_pAu = 7;
	const int ncentbin_pAl = 3;
	const int ncentbin_HeAu = 3;

	const int cent_per_pAu_min[ncentbin_pAu] = {0, 5, 10, 0, 20, 40, 60};
	const int cent_per_pAu_max[ncentbin_pAu] = {5, 10, 20, 20, 40, 60, 84};
	const float scaler_pAu[ncentbin_pAu] = {1, 1e-1, 1e-2, 1e-2, 1e-3, 1e-4, 1e-5};
	const int nMarker_pAu[ncentbin_pAu] = {20, 24, 21, 25, 25, 22, 26};
	const int nColor_pAu[ncentbin_pAu] = {1, 2, 4, 6, 6, kGreen+2, kOrange+6};

	const int cent_per_pAl_min[ncentbin_pAl] = {0, 20, 40};
	const int cent_per_pAl_max[ncentbin_pAl] = {20, 40, 72};
	const float scaler_pAl[ncentbin_pAl] = {1, 1e-1, 1e-2};
	const int nMarker_pAl[ncentbin_pAl] = {20, 24, 21};
	const int nColor_pAl[ncentbin_pAl] = {1, 2, 4};

	const int cent_per_HeAu_min[ncentbin_HeAu] = {0, 20, 40};
	const int cent_per_HeAu_max[ncentbin_HeAu] = {20, 40, 88};
	const float scaler_HeAu[ncentbin_HeAu] = {1, 1e-1, 1e-2};
	const int nMarker_HeAu[ncentbin_HeAu] = {20, 24, 21};
	const int nColor_HeAu[ncentbin_HeAu] = {1, 2, 4};

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

	for (int iarm=0; iarm<narm; iarm++){
		for (int icent=0; icent<ncentbin_pAu; icent++){
			YpAu[icent][iarm] = (TH1F*)file_pAu->Get(Form("YpAu_vs_pT_%c_%d%d",c_arm[iarm],cent_per_pAu_min[icent],cent_per_pAu_max[icent]));
			YpAu_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_pAu->Get(Form("YpAu_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_pAu_min[icent],cent_per_pAu_max[icent]));
		}

		for (int icent=0; icent<ncentbin_pAl; icent++){
			YpAl[icent][iarm] = (TH1F*)file_pAl->Get(Form("YpAl_vs_pT_%c_%d%d",c_arm[iarm],cent_per_pAl_min[icent],cent_per_pAl_max[icent]));
			YpAl_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_pAl->Get(Form("YpAl_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_pAl_min[icent],cent_per_pAl_max[icent]));
		}

		for (int icent=0; icent<ncentbin_HeAu; icent++){
			YHeAu[icent][iarm] = (TH1F*)file_HeAu->Get(Form("YHeAu_vs_pT_%c_%d%d",c_arm[iarm],cent_per_HeAu_min[icent],cent_per_HeAu_max[icent]));
			YHeAu_SYS_FRAC_ERR[icent][iarm] = (TH1F*)file_HeAu->Get(Form("YHeAu_vs_pT_%c_%d%d_sys_frac",c_arm[iarm],cent_per_HeAu_min[icent],cent_per_HeAu_max[icent]));
		}
	}

	//return;

	TGraphErrors *gr_pAu[ncentbin_pAu][narm];
	TGraphErrors *gr_pAl[ncentbin_pAl][narm];
	TGraphErrors *gr_HeAu[ncentbin_HeAu][narm];

	TGraphErrors *gr_pAu_sys[ncentbin_pAu][narm];
	TGraphErrors *gr_pAl_sys[ncentbin_pAl][narm];
	TGraphErrors *gr_HeAu_sys[ncentbin_HeAu][narm];

	//plot arrays
	double pt_array_pAu[nptbin_pAu] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.25,4.75,6.0};
	double pt_array_pAl[nptbin_pAl] =  {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875};
	double pt_array_HeAu[nptbin_HeAu] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75};

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

			for (int ipt=0; ipt<nptbin_pAu; ipt++){
				double xx = pt_array_pAu[ipt];
				double yy = YpAu[icent][iarm]->GetBinContent(ipt+1);
				double yy_err = YpAu[icent][iarm]->GetBinError(ipt+1);
				double yy_syserr = yy*YpAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				gr_pAu[icent][iarm]->SetPoint(ipt, pt_array_pAu[ipt], yy*scaler_pAu[icent]);
				gr_pAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAu[icent]);

				gr_pAu_sys[icent][iarm]->SetPoint(ipt, pt_array_pAu[ipt], yy*scaler_pAu[icent]);
				gr_pAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAu[icent]);
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

			for (int ipt=0; ipt<nptbin_pAl; ipt++){
				double xx = pt_array_pAl[ipt];
				double yy = YpAl[icent][iarm]->GetBinContent(ipt+1);
				double yy_err = YpAl[icent][iarm]->GetBinError(ipt+1);
				double yy_syserr = yy*YpAl_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				gr_pAl[icent][iarm]->SetPoint(ipt, pt_array_pAl[ipt], yy*scaler_pAl[icent]);
				gr_pAl[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_pAl[icent]);

				gr_pAl_sys[icent][iarm]->SetPoint(ipt, pt_array_pAl[ipt], yy*scaler_pAl[icent]);
				gr_pAl_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_pAl[icent]);
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

			for (int ipt=0; ipt<nptbin_HeAu; ipt++){
				double xx = pt_array_HeAu[ipt];
				double yy = YHeAu[icent][iarm]->GetBinContent(ipt+1);
				double yy_err = YHeAu[icent][iarm]->GetBinError(ipt+1);
				double yy_syserr = yy*YHeAu_SYS_FRAC_ERR[icent][iarm]->GetBinError(ipt+1);
				gr_HeAu[icent][iarm]->SetPoint(ipt, pt_array_HeAu[ipt], yy*scaler_HeAu[icent]);
				gr_HeAu[icent][iarm]->SetPointError(ipt, 0, yy_err*scaler_HeAu[icent]);

				gr_HeAu_sys[icent][iarm]->SetPoint(ipt, pt_array_HeAu[ipt], yy*scaler_HeAu[icent]);
				gr_HeAu_sys[icent][iarm]->SetPointError(ipt, 0.25/2, yy_syserr*scaler_HeAu[icent]);
			}//ipt
		}//iarm
	}//icent

	TCanvas *c00[narm];;

	for (int iarm=0; iarm<narm; iarm++){
		c00[iarm] = new TCanvas(Form("c00_arm%d",iarm),Form("c00_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.13);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-18,7.8,5e-6);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",20,18);
		htmp->GetYaxis()->SetTitleOffset(1.4);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.18,0.15,0.73,0.32);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetNColumns(2);
			for (int icent=0; icent<ncentbin_pAu; icent++){
				if ( icent==3 ){
					continue;
				}
				gr_pAu_sys[icent][iarm]->Draw("2");
				leg->AddEntry(gr_pAu_sys[icent][iarm],Form("%d%c-%d%c#times10^{%d}",cent_per_pAu_min[icent],'%',cent_per_pAu_max[icent],'%',int(TMath::Log10(scaler_pAu[icent])-0.00001)),"P");
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
			TLegend *leg = new TLegend(0.65,0.77,0.9,0.92);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","p+Au #sqrt{s_{NN}}=200 GeV","h");
			le = leg->AddEntry("","Inclusive J/#psi","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y<2.2","h");
			}
			leg->Draw();
		}
	}

	TCanvas *c01[narm];;

	for (int iarm=0; iarm<narm; iarm++){
		c01[iarm] = new TCanvas(Form("c01_arm%d",iarm),Form("c01_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.13);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-12,4.8,5e-6);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",20,18);
		htmp->GetYaxis()->SetTitleOffset(1.4);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.2,0.15,0.65,0.32);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			for (int icent=0; icent<ncentbin_pAl; icent++){
				gr_pAl_sys[icent][iarm]->Draw("2");
				leg->AddEntry(gr_pAl_sys[icent][iarm],Form("%d%c-%d%c#times10^{%d}",cent_per_pAl_min[icent],'%',cent_per_pAl_max[icent],'%',int(TMath::Log10(scaler_pAl[icent])-0.00001)),"P");
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
			TLegend *leg = new TLegend(0.65,0.77,0.9,0.92);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","p+Al #sqrt{s_{NN}}=200 GeV","h");
			le = leg->AddEntry("","Inclusive J/#psi","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y<2.2","h");
			}
			leg->Draw();
		}
	}


	TCanvas *c02[narm];;

	for (int iarm=0; iarm<narm; iarm++){
		c02[iarm] = new TCanvas(Form("c02_arm%d",iarm),Form("c02_arm%d",iarm),1.3*1*400,1*400);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetRightMargin(0.02);
		gPad->SetLeftMargin(0.17);
		gPad->SetBottomMargin(0.13);
		htmp = (TH1F*)gPad->DrawFrame(0,5e-12,4.8,5e-6);
		SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",20,18);
		htmp->GetYaxis()->SetTitleOffset(1.4);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.2,0.15,0.65,0.32);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
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
			TLegend *leg = new TLegend(0.63,0.77,0.9,0.92);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","^{3}He+Au #sqrt{s_{NN}}=200 GeV","h");
			le = leg->AddEntry("","Inclusive J/#psi","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y<2.2","h");
			}
			leg->Draw();
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



