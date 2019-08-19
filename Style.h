#include <TH1.h>
#include <TPad.h>
#include <TStyle.h>
#include <string>

TH1 *htmp;
TLegendEntry *le;

void SetStyle(){

	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLegendFont(43);
	gStyle->SetLegendTextSize(18);

}

void SetPadStyle(bool log = false){
	gPad->SetLeftMargin(0.18);
	if ( log ){
		gPad->SetRightMargin(0.12);
	}else{
		gPad->SetRightMargin(0.04);
	}
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.14);

	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);

	return;
}

void SetHistoStyle(std::string xtitle="", std::string ytitle="", std::string ztitle="",float size_title=22, float size_label=18){
	htmp->SetTitleFont(62);
	htmp->GetYaxis()->SetLabelFont(63);
	htmp->GetYaxis()->SetTitleFont(63);
	htmp->GetYaxis()->SetTitleOffset(1.25);
	htmp->GetYaxis()->SetLabelSize(size_label);
	htmp->GetYaxis()->SetTitleSize(size_title);
	htmp->GetYaxis()->SetTitle(ytitle.c_str());
	htmp->GetYaxis()->SetNdivisions(9,5,0);
	htmp->GetXaxis()->SetLabelFont(63);
	htmp->GetXaxis()->SetTitleFont(63);
	htmp->GetXaxis()->SetTitleOffset(1.0);
	htmp->GetXaxis()->SetLabelSize(size_label);
	htmp->GetXaxis()->SetTitleSize(size_title);
	htmp->GetXaxis()->SetTitle(xtitle.c_str());
	htmp->GetXaxis()->SetNdivisions(9,2,0);
	htmp->GetZaxis()->SetLabelFont(63);
	htmp->GetZaxis()->SetTitleFont(63);
	htmp->GetZaxis()->SetTitleOffset(1.0);
	htmp->GetZaxis()->SetLabelSize(size_label);
	htmp->GetZaxis()->SetTitleSize(size_title);
	htmp->GetZaxis()->SetTitle(ztitle.c_str());

	return;
}

