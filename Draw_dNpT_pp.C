#include "Style.h"

void Draw_dNpT_pp()
{

	const bool bSAVE = false;

	SetStyle();

	TFile *file_pp = new TFile("Ypp_vs_pT.root","read");

	const int narm = 2;
	const char c_arm[narm] = {'S', 'N'};
	const int nMarker_pp[narm] = {20, 24};
	const int nColor_pp[narm] = {1, 2};

	const float pol_sys_pp[narm] = {0.12, 0.12};

	TH1F *Ypp[narm];
	TH1F *Ypp_SYS_FRAC_ERR[narm];

	for (int iarm=0; iarm<narm; iarm++){
		Ypp[iarm] = (TH1F*)file_pp->Get(Form("Ypp_vs_pT_%c",c_arm[iarm]));
		Ypp_SYS_FRAC_ERR[iarm] = (TH1F*)file_pp->Get(Form("Ypp_vs_pT_%c_frac_sys",c_arm[iarm]));
	}

	TGraphErrors *gr_pp[narm];
	TGraphErrors *gr_pp_sys[narm];

	for (int iarm=0; iarm<narm; iarm++){
		gr_pp[iarm] = new TGraphErrors;
		gr_pp[iarm]->SetMarkerStyle(nMarker_pp[iarm]);
		gr_pp[iarm]->SetMarkerColor(nColor_pp[iarm]);
		gr_pp[iarm]->SetMarkerSize(0.8);
		gr_pp[iarm]->SetLineColor(nColor_pp[iarm]);

		gr_pp_sys[iarm] = new TGraphErrors;
		gr_pp_sys[iarm]->SetMarkerStyle(nMarker_pp[iarm]);
		gr_pp_sys[iarm]->SetMarkerColor(nColor_pp[iarm]);
		gr_pp_sys[iarm]->SetMarkerSize(0.8);
		gr_pp_sys[iarm]->SetLineColor(nColor_pp[iarm]);
		gr_pp_sys[iarm]->SetFillStyle(0);

		for (int ipt=0; ipt<Ypp[iarm]->GetNbinsX(); ipt++){
			double xx = Ypp[iarm]->GetBinCenter(ipt+1);
			double yy = Ypp[iarm]->GetBinContent(ipt+1);
			double yy_err = Ypp[iarm]->GetBinError(ipt+1);
			double yy_syserr = yy*Ypp_SYS_FRAC_ERR[iarm]->GetBinError(ipt+1);
			gr_pp[iarm]->SetPoint(ipt, xx, yy);
			gr_pp[iarm]->SetPointError(ipt, 0, yy_err);

			gr_pp_sys[iarm]->SetPoint(ipt, xx, yy);
			gr_pp_sys[iarm]->SetPointError(ipt, 0.25/2, yy_syserr);
		}//ipt
	}//iarm

	TGraphErrors *gr_ratio = new TGraphErrors;
	gr_ratio->SetMarkerStyle(20);
	gr_ratio->SetMarkerColor(1);
	gr_ratio->SetLineColor(1);

	TGraphErrors *gr_ratio_sys = new TGraphErrors;
	gr_ratio_sys->SetMarkerStyle(20);
	gr_ratio_sys->SetMarkerColor(1);
	gr_ratio_sys->SetLineColor(1);
	gr_ratio_sys->SetFillStyle(0);

	for (int ipt=0; ipt<Ypp[0]->GetNbinsX(); ipt++){
		double xx = Ypp[0]->GetBinCenter(ipt+1);
		double yy0 = Ypp[0]->GetBinContent(ipt+1);
		double yy1 = Ypp[1]->GetBinContent(ipt+1);
		double yy0_err = Ypp[0]->GetBinError(ipt+1);
		double yy1_err = Ypp[1]->GetBinError(ipt+1);

		double yy0_syserr = Ypp_SYS_FRAC_ERR[0]->GetBinError(ipt+1);
		double yy1_syserr = Ypp_SYS_FRAC_ERR[1]->GetBinError(ipt+1);

		double ratio = yy0 / yy1;
		double err = ratio * sqrt(pow(yy0_err/yy0,2) + pow(yy1_err/yy1,2));
		double sys_err = ratio * sqrt(pow(yy0_syserr,2) + pow(yy1_syserr,2) - pow(pol_sys_pp[0],2) - pow(pol_sys_pp[1],2));

		gr_ratio->SetPoint(ipt, xx, ratio); 
		gr_ratio->SetPointError(ipt, 0, err); 

		gr_ratio_sys->SetPoint(ipt, xx, ratio); 
		gr_ratio_sys->SetPointError(ipt, 0.25/2, sys_err); 
	}

	TCanvas *c00 = new TCanvas("c00","c00",1.3*400,600);
	TPad *p00 = new TPad("p00","p00",0,0.3333,1.0,1.0);
	p00->Draw();
	p00->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.02);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.0001);
	htmp = (TH1F*)gPad->DrawFrame(0,5e-12,7.8,5e-7);
	SetHistoStyle("p_{T} (GeV/c)","B_{ll} d^{2}N/(2#pip_{T}dp_{T}dy) (GeV/c)^{-2}","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.8);
	htmp->GetXaxis()->SetTitleOffset(1.0);

	gr_pp_sys[0]->Draw("2");
	gr_pp_sys[1]->Draw("2");
	gr_pp[0]->Draw("p");
	gr_pp[1]->Draw("p");

	{
		TLegend *leg = new TLegend(0.55,0.65,0.9,0.90);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		le = leg->AddEntry("","p+p #sqrt{s_{NN}}=200 GeV","");
		le = leg->AddEntry("","Inclusive J/#psi","");
		leg->AddEntry(gr_pp[0],"-2.2<y<-1.2","p");
		leg->AddEntry(gr_pp[1],"1.2<y<2.2","p");
		leg->Draw();
	}

	c00->cd();
	TPad *p01 = new TPad("p01","p01",0,0.0,1.0,0.33333);
	p01->Draw();
	p01->cd();
	gPad->SetRightMargin(0.02);
	gPad->SetTopMargin(0.000000);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.24);
	htmp = (TH1F*)gPad->DrawFrame(0,0.5,7.8,1.5);
	SetHistoStyle("p_{T} (GeV/c)","Ratio","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.8);
	htmp->GetYaxis()->CenterTitle();
	htmp->GetXaxis()->SetTitleOffset(3.0);

	gr_ratio_sys->Draw("2");
	gr_ratio->Draw("p");

	if ( bSAVE ){
		c00->cd();
		c00->SaveAs("pdf/fig_dNdpT_pp.pdf");
	}


}
