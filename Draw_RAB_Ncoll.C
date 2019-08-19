#include "Style.h"

void Draw_RAB_Ncoll(){

	SetStyle();

	TFile *infile_pAu = new TFile("RAB_y_Ncoll_pAu200.root","read");
	TFile *infile_pAl = new TFile("RAB_y_Ncoll_pAl200.root","read");
	TFile *infile_HeAu = new TFile("RAB_y_Ncoll_HeAu200.root","read");

	const int narm = 2;

	TH1D *hRAB_Ncoll_pAu[narm];
	TH1D *hRAB_Ncoll_pAl[narm];
	TH1D *hRAB_Ncoll_HeAu[narm];

	TH1D *hRAB_Ncoll_pAu_sys[narm];
	TH1D *hRAB_Ncoll_pAl_sys[narm];
	TH1D *hRAB_Ncoll_HeAu_sys[narm];

	for (int iarm=0; iarm<narm; iarm++){

		hRAB_Ncoll_pAu_sys[iarm] = (TH1D*)infile_pAu->Get(Form("RAB_Ncoll_syserr_arm%d",iarm));
		hRAB_Ncoll_pAu[iarm] = (TH1D*)infile_pAu->Get(Form("RAB_Ncoll_staerr_arm%d",iarm));
		hRAB_Ncoll_pAu_sys[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_pAu_sys[iarm]->SetMarkerColor(1);
		hRAB_Ncoll_pAu_sys[iarm]->SetLineColor(1);
		hRAB_Ncoll_pAu[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_pAu[iarm]->SetMarkerColor(1);
		hRAB_Ncoll_pAu[iarm]->SetLineColor(1);

		hRAB_Ncoll_pAl_sys[iarm] = (TH1D*)infile_pAl->Get(Form("RAB_Ncoll_syserr_arm%d",iarm));
		hRAB_Ncoll_pAl[iarm] = (TH1D*)infile_pAl->Get(Form("RAB_Ncoll_staerr_arm%d",iarm));
		hRAB_Ncoll_pAl_sys[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_pAl_sys[iarm]->SetMarkerColor(2);
		hRAB_Ncoll_pAl_sys[iarm]->SetLineColor(2);
		hRAB_Ncoll_pAl[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_pAl[iarm]->SetMarkerColor(2);
		hRAB_Ncoll_pAl[iarm]->SetLineColor(2);

		hRAB_Ncoll_HeAu_sys[iarm] = (TH1D*)infile_HeAu->Get(Form("RAB_Ncoll_syserr_arm%d",iarm));
		hRAB_Ncoll_HeAu[iarm] = (TH1D*)infile_HeAu->Get(Form("RAB_Ncoll_staerr_arm%d",iarm));
		hRAB_Ncoll_HeAu_sys[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_HeAu_sys[iarm]->SetMarkerColor(4);
		hRAB_Ncoll_HeAu_sys[iarm]->SetLineColor(4);
		hRAB_Ncoll_HeAu[iarm]->SetMarkerStyle(20);
		hRAB_Ncoll_HeAu[iarm]->SetMarkerColor(4);
		hRAB_Ncoll_HeAu[iarm]->SetLineColor(4);

	}//iarm

	TCanvas *c13[narm];
	for (int iarm=0; iarm<narm; iarm++){
		c13[iarm] = new TCanvas(Form("c13_arm%d",iarm),Form("c13_arm%d",iarm),1.3*1*400,400);
		SetPadStyle();
		gPad->SetRightMargin(0.03);
		gPad->SetLeftMargin(0.14);
		gPad->SetBottomMargin(0.16);
		htmp = (TH1F*)gPad->DrawFrame(0,0,25,2.0);
		SetHistoStyle("#LTN_{coll}#GT","R_{AB}","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.0);
		htmp->GetXaxis()->SetTitleOffset(1.0);

		{
			TLegend *leg = new TLegend(0.10,0.73,0.4,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("","Inclusive J/#psi","");
			le = leg->AddEntry("","#sqrt{s_{NN}}=200 GeV","");
			if ( iarm==0 ){
				le = leg->AddEntry("","-2.2<y<-1.2, Au/Al-going","");
			}else{
				le = leg->AddEntry("","1.2<y<2.2, p/^{3}He-going","");
			}
			leg->Draw();

			TLine *line = new TLine(0, 1, 25, 1);
			line->SetLineStyle(2);
			line->Draw();

			TBox *box = new TBox(25*0.97,1-0.101,25,1+0.101);
			box->SetFillStyle(1000);
			box->SetFillColorAlpha(1,0.5);
			box->Draw();
		}
		{
			TLegend *leg = new TLegend(0.65,0.73,0.9,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->AddEntry(hRAB_Ncoll_pAl[iarm],"p+Al","P");
			leg->AddEntry(hRAB_Ncoll_pAu[iarm],"p+Au","P");
			leg->AddEntry(hRAB_Ncoll_HeAu[iarm],"^{3}He+Au","P");
			leg->Draw();
		}

		hRAB_Ncoll_pAu_sys[iarm]->Draw("2");
		hRAB_Ncoll_pAu[iarm]->Draw("P");

		hRAB_Ncoll_pAl_sys[iarm]->Draw("2");
		hRAB_Ncoll_pAl[iarm]->Draw("P");

		hRAB_Ncoll_HeAu_sys[iarm]->Draw("2");
		hRAB_Ncoll_HeAu[iarm]->Draw("P");
	}


}
