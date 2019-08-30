// to change colors of data points, change the TGRaph Error definitions for RpAu_bkwd group and RpAu_fwrd group at lines 300 and 767 
#include "Style.h"


void Draw_ptsq()
{

	SetStyle();

	int i_cent = 0; // 0 = 0-20 (or for pT squared plot), 1 = 20-40,2 = 40-60 (or 40-88), 3 = 60-84 

	bool draw_pAu = 0;
	bool draw_dAu = 0;
	bool draw_HeAu = 0;
	bool draw_pAl = 0;
	bool draw_NS = 0;
	bool pt_squared = 1;

	//	double ncoll_pAu[7] = {9.7, 8.4, 7.4, 8.2, 6.1, 4.4, 2.6};
	double ncoll_pAu[6] = {9.7, 8.4, 7.4, 6.1, 4.4, 2.6};
	double ncoll_HeAu[3] = {22.3,14.8,5.4833};
	double ncoll_pAl[3] = {3.35, 2.3, 1.7};

	// pt < 4
	// double pAu_N_ptsq[6] = {3.57363, 3.53415, 3.55327, 3.48693, 3.3936, 3.25398};
	// double pAu_S_ptsq[6] = {3.60433, 3.30734, 3.3847, 3.42589, 3.33675, 3.38347};

	// double pAu_N_stat[6] = {0.212155, 0.212712, 0.157712, 0.124755, 0.134495, 0.147847};
	// double pAu_S_stat[6] = {0.241367, 0.226709, 0.184654, 0.148624, 0.197885, 0.234599};

	// double pAu_N_sys_up[6] = {0.234261, 0.235853, 0.23978, 0.226316, 0.233597, 0.224206};
	// double pAu_S_sys_up[6] = {0.239434, 0.224555, 0.244746, 0.242928, 0.22707, 0.246628};

	// double pAu_N_sys_down[6] = {0.309948, 0.311179, 0.316446, 0.297948, 0.3056, 0.290851};
	// double pAu_S_sys_down[6] = {0.324666, 0.299797, 0.32719, 0.323413, 0.302212, 0.327299};
 
	// double HeAu_N_ptsq[3] = {3.48699, 3.31614, 3.33336};
	// double HeAu_S_ptsq[3] = { 3.52125, 3.45499, 3.40047};

	// double HeAu_N_stat[3] = {0.217861, 0.238833, 0.244705};
	// double HeAu_S_stat[3] = {0.193704, 0.24534, 0.241158};

	// double HeAu_N_sys_up[3] = {0.264157, 0.254954, 0.244728};
	// double HeAu_S_sys_up[3] = {0.254302, 0.262157, 0.263699};

	// double HeAu_N_sys_down[3] = {0.360241, 0.344067, 0.331603};
	// double HeAu_S_sys_down[3] = {0.353651, 0.361043, 0.362704};

	// double pAl_N_ptsq[3] = {3.22431, 3.32582, 3.07016};
	// double pAl_S_ptsq[3] = {3.12468, 2.96197, 3.12653};

	// double pAl_N_stat[3] = {0.129216, 0.139618, 0.132623};
	// double pAl_S_stat[3] = {0.161672, 0.169745, 0.171909};

	// double pAl_N_sys_up[3] = {0.23819, 0.238708, 0.221};
	// double pAl_S_sys_up[3] = {0.224095, 0.214787, 0.238387};

	// double pAl_N_sys_down[3] = {0.314817, 0.317681, 0.29006};
	// double pAl_S_sys_down[3] = {0.295753, 0.280596, 0.313873};

	// pt < 5
	double pAu_N_ptsq[6] = {4.01324, 4.03177, 4.00273, 3.92388, 3.85576, 3.65299};
	double pAu_S_ptsq[6] = {4.08, 3.66752, 3.76513, 3.94679, 3.69379, 3.82217};

	double pAu_N_stat[6] = {0.231953, 0.235493, 0.172977, 0.136692, 0.148697, 0.162113};
	double pAu_S_stat[6] = {0.265216, 0.246243, 0.200759, 0.165873, 0.214349, 0.25816};

	double pAu_N_sys_up[6] = {0.277673, 0.261262, 0.272526, 0.275656, 0.254359, 0.267257};
	double pAu_S_sys_up[6] = {0.274766, 0.248793, 0.266501, 0.247403, 0.277673, 0.261262};

	double pAu_N_sys_down[6] = {0.329786, 0.343737, 0.339176, 0.319799, 0.334544, 0.310229};
	double pAu_S_sys_down[6] = {0.354267, 0.315867, 0.338935, 0.363873, 0.313306, 0.351577};
 
	double HeAu_N_ptsq[3] = {3.98446, 3.73778, 3.71093};
	double HeAu_S_ptsq[3] = {3.89073, 3.82153, 3.72857};

	double HeAu_N_stat[3] = {0.241687, 0.262601, 0.266468};
	double HeAu_S_stat[3] = {0.209392, 0.26554, 0.259419};

	double HeAu_N_sys_up[3] = {0.28509, 0.267279, 0.250905};
	double HeAu_S_sys_up[3] = {0.261739, 0.267693, 0.260072};

	double HeAu_N_sys_down[3] = {0.361938, 0.336137, 0.315885};
	double HeAu_S_sys_down[3] = {0.337256, 0.342448, 0.331127};

	double pAl_N_ptsq[3] =  {3.39899, 3.50583, 3.23112};
	double pAl_S_ptsq[3] = {3.22894, 3.11982, 3.22337};

	 double pAl_N_stat[3] = {0.134829, 0.145627, 0.138283};
	 double pAl_S_stat[3] = {0.166052, 0.177178, 0.176227};

	 double pAl_N_sys_up[3] = {0.216258, 0.217217, 0.200939};
	 double pAl_S_sys_up[3] = {0.193271, 0.196044, 0.202934};

	 double pAl_N_sys_down[3] = {0.267927, 0.270537, 0.247566};
	 double pAl_S_sys_down[3] = {0.238774, 0.240669, 0.250167};

	double HeAu_xlow[3] = {0};
	double HeAu_xhigh[3] = {0};
	double pAl_xlow[3] = {0};
	double pAl_xhigh[3] = {0};
	double pAu_xlow[6] = {0};
	double pAu_xhigh[6] = {0};
	double HeAu_ylow[3] = {0};
	double HeAu_yhigh[3] = {0};
	double pAu_ylow[6] = {0};
	double pAu_yhigh[6] = {0};
	double pAl_ylow[3] = {0};
	double pAl_yhigh[3] = {0};

	double width = 0.5;

	for(int i = 0; i < 3; i++)
	{
		HeAu_xlow[i] = width;
		HeAu_xhigh[i] = width;
		HeAu_ylow[i] = HeAu_S_sys_down[i];
		HeAu_yhigh[i] = HeAu_S_sys_up[i] ;
		cout << HeAu_xlow[i] << "," << HeAu_xhigh[i] << "," << HeAu_ylow[i] << "," << HeAu_yhigh[i] << endl;
	}
	for(int i = 0; i < 3; i++)
	{
		pAl_xlow[i] = width;
		pAl_xhigh[i] = width;
		pAl_ylow[i] = pAl_S_sys_down[i];
		pAl_yhigh[i] = pAl_S_sys_up[i] ;
		cout << pAl_xlow[i] << "," << pAl_xhigh[i] << "," << pAl_ylow[i] << "," << pAl_yhigh[i] << endl;
	}

	for(int i = 0; i < 6; i++)
	{
		pAu_xlow[i] = width;
		pAu_xhigh[i] = width;
		pAu_ylow[i] = pAu_S_sys_down[i];
		pAu_yhigh[i] = pAu_S_sys_up[i];
	}


	double pt_dAu[24] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.25,5.75,6.50,7.5}; 

	double pt_pAu[19] = {0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 6}; // for R_HeAu

	double pt_pAl[18] = {0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.5, 6.5}; // for R_pAl

	double pt_HeAu[15] = {0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.5, 6.0};


	//pt^2
	if ( pt_squared ){
		TGraphAsymmErrors *Ncoll_HeAu_S_sys = new TGraphAsymmErrors(3,ncoll_HeAu,HeAu_S_ptsq,HeAu_xlow,HeAu_xhigh,HeAu_S_sys_down,HeAu_S_sys_up);    TGraphErrors *Ncoll_HeAu_S = new TGraphErrors(3);

		Ncoll_HeAu_S_sys->SetMarkerColor(kBlue);                         Ncoll_HeAu_S->SetMarkerColor(kBlue);
		Ncoll_HeAu_S_sys->SetLineColor(kBlue);                           Ncoll_HeAu_S->SetLineColor(kBlue);
		Ncoll_HeAu_S_sys->SetFillColor(kBlue);                           Ncoll_HeAu_S->SetFillColor(kBlue); 
		Ncoll_HeAu_S_sys->SetLineWidth(1);                           Ncoll_HeAu_S->SetLineWidth(1);
		Ncoll_HeAu_S_sys->SetMarkerStyle(20);                        Ncoll_HeAu_S->SetMarkerStyle(20);
		Ncoll_HeAu_S_sys->SetMarkerSize(1);                        Ncoll_HeAu_S->SetMarkerSize(1.);
		Ncoll_HeAu_S_sys->SetFillStyle(0);

		TGraphAsymmErrors *Ncoll_pAl_S_sys = new TGraphAsymmErrors(3,ncoll_pAl,pAl_S_ptsq,pAl_xlow,pAl_xhigh,pAl_S_sys_down,pAl_S_sys_up);    TGraphErrors *Ncoll_pAl_S = new TGraphErrors(3);

		Ncoll_pAl_S_sys->SetMarkerColor(kBlack);                         Ncoll_pAl_S->SetMarkerColor(kBlack);
		Ncoll_pAl_S_sys->SetLineColor(kBlack);                           Ncoll_pAl_S->SetLineColor(kBlack);
		Ncoll_pAl_S_sys->SetFillColor(kBlack);                           Ncoll_pAl_S->SetFillColor(kBlack); 
		Ncoll_pAl_S_sys->SetLineWidth(1);                           Ncoll_pAl_S->SetLineWidth(1);
		Ncoll_pAl_S_sys->SetMarkerStyle(20);                        Ncoll_pAl_S->SetMarkerStyle(20);
		Ncoll_pAl_S_sys->SetMarkerSize(1.);                        Ncoll_pAl_S->SetMarkerSize(1.);
		Ncoll_pAl_S_sys->SetFillStyle(0);

		TGraphAsymmErrors *Ncoll_pAu_S_sys = new TGraphAsymmErrors(6,ncoll_pAu,pAu_S_ptsq,pAu_xlow,pAu_xhigh,pAu_S_sys_down,pAu_S_sys_up);     TGraphErrors *Ncoll_pAu_S = new TGraphErrors(6);

		Ncoll_pAu_S_sys->SetMarkerColor(kRed);                         Ncoll_pAu_S->SetMarkerColor(kRed);
		Ncoll_pAu_S_sys->SetLineColor(kRed);                           Ncoll_pAu_S->SetLineColor(kRed);
		Ncoll_pAu_S_sys->SetFillColor(kRed);                           Ncoll_pAu_S->SetFillColor(kRed); 
		Ncoll_pAu_S_sys->SetLineWidth(1);                           Ncoll_pAu_S->SetLineWidth(1);
		Ncoll_pAu_S_sys->SetMarkerStyle(20);                        Ncoll_pAu_S->SetMarkerStyle(20);
		Ncoll_pAu_S_sys->SetMarkerSize(1.);                        Ncoll_pAu_S->SetMarkerSize(1.);
		Ncoll_pAu_S_sys->SetFillStyle(0);

		for(int i = 0; i < 6; i++)
		{
			Ncoll_pAu_S->SetPoint(i,ncoll_pAu[i],pAu_S_ptsq[i]);  
			Ncoll_pAu_S->SetPointError(i,0,pAu_S_stat[i]); 

			if(i < 3)
			{
				Ncoll_HeAu_S->SetPoint(i,ncoll_HeAu[i],HeAu_S_ptsq[i]); 
				Ncoll_HeAu_S->SetPointError(i,0,HeAu_S_stat[i]); 
				Ncoll_pAl_S->SetPoint(i,ncoll_pAl[i],pAl_S_ptsq[i]); 
				Ncoll_pAl_S->SetPointError(i,0,pAl_S_stat[i]); 
			}
		}

		TGraphAsymmErrors *Ncoll_HeAu_N_sys = new TGraphAsymmErrors(3,ncoll_HeAu,HeAu_N_ptsq,HeAu_xlow,HeAu_xhigh,HeAu_N_sys_down,HeAu_N_sys_up);     TGraphErrors *Ncoll_HeAu_N = new TGraphErrors(3);

		Ncoll_HeAu_N_sys->SetMarkerColor(kBlue);                         Ncoll_HeAu_N->SetMarkerColor(kBlue);
		Ncoll_HeAu_N_sys->SetLineColor(kBlue);                           Ncoll_HeAu_N->SetLineColor(kBlue);
		Ncoll_HeAu_N_sys->SetFillColor(kBlue);                           Ncoll_HeAu_N->SetFillColor(kBlue); 
		Ncoll_HeAu_N_sys->SetLineWidth(1);                           Ncoll_HeAu_N->SetLineWidth(1);
		Ncoll_HeAu_N_sys->SetMarkerStyle(20);                        Ncoll_HeAu_N->SetMarkerStyle(20);
		Ncoll_HeAu_N_sys->SetMarkerSize(1.);                        Ncoll_HeAu_N->SetMarkerSize(1.);
		Ncoll_HeAu_N_sys->SetFillStyle(0);

		TGraphAsymmErrors *Ncoll_pAl_N_sys = new TGraphAsymmErrors(3,ncoll_pAl,pAl_N_ptsq,pAl_xlow,pAl_xhigh,pAl_N_sys_down,pAl_N_sys_up);     TGraphErrors *Ncoll_pAl_N = new TGraphErrors(3);

		Ncoll_pAl_N_sys->SetMarkerColor(kBlack);                         Ncoll_pAl_N->SetMarkerColor(kBlack);
		Ncoll_pAl_N_sys->SetLineColor(kBlack);                           Ncoll_pAl_N->SetLineColor(kBlack);
		Ncoll_pAl_N_sys->SetFillColor(kBlack);                           Ncoll_pAl_N->SetFillColor(kBlack); 
		Ncoll_pAl_N_sys->SetLineWidth(1);                           Ncoll_pAl_N->SetLineWidth(1);
		Ncoll_pAl_N_sys->SetMarkerStyle(20);                        Ncoll_pAl_N->SetMarkerStyle(20);
		Ncoll_pAl_N_sys->SetMarkerSize(1.);                        Ncoll_pAl_N->SetMarkerSize(1.);
		Ncoll_pAl_N_sys->SetFillStyle(0);

		TGraphAsymmErrors *Ncoll_pAu_N_sys = new TGraphAsymmErrors(6,ncoll_pAu,pAu_N_ptsq,pAu_xlow,pAu_xhigh,pAu_N_sys_down,pAu_N_sys_up);     TGraphErrors *Ncoll_pAu_N = new TGraphErrors(6);

		Ncoll_pAu_N_sys->SetMarkerColor(kRed);                         Ncoll_pAu_N->SetMarkerColor(kRed);
		Ncoll_pAu_N_sys->SetLineColor(kRed);                           Ncoll_pAu_N->SetLineColor(kRed);
		Ncoll_pAu_N_sys->SetFillColor(kRed);                           Ncoll_pAu_N->SetFillColor(kRed); 
		Ncoll_pAu_N_sys->SetLineWidth(1);                           Ncoll_pAu_N->SetLineWidth(1);
		Ncoll_pAu_N_sys->SetMarkerStyle(20);                        Ncoll_pAu_N->SetMarkerStyle(20);
		Ncoll_pAu_N_sys->SetMarkerSize(1.);                        Ncoll_pAu_N->SetMarkerSize(1.);
		Ncoll_pAu_N_sys->SetFillStyle(0);

		for(int i = 0; i < 6; i++)
		{
			Ncoll_pAu_N->SetPoint(i,ncoll_pAu[i],pAu_N_ptsq[i]);  // this appears correct.
			Ncoll_pAu_N->SetPointError(i,0,pAu_N_stat[i]); // TGraph with symemtric Type A errors
			if(i < 3)
			{
				Ncoll_HeAu_N->SetPoint(i,ncoll_HeAu[i],HeAu_N_ptsq[i]); 
				Ncoll_HeAu_N->SetPointError(i,0,HeAu_N_stat[i]); 
				Ncoll_pAl_N->SetPoint(i,ncoll_pAl[i],pAl_N_ptsq[i]); 
				Ncoll_pAl_N->SetPointError(i,0,pAl_N_stat[i]); 
			}
		}

		TCanvas *c00[2];
		for (int iarm=0; iarm<2; iarm++){
			c00[iarm] = new TCanvas(Form("c00_arm%d",iarm),Form("c00_arm%d",iarm),1.3*400,400);
			c00[iarm]->cd(iarm+1);
			SetPadStyle();
			gPad->SetRightMargin(0.03);
			gPad->SetLeftMargin(0.12);
			gPad->SetBottomMargin(0.16);
			htmp = (TH1F*)gPad->DrawFrame(0,0,25,6.5);
			SetHistoStyle("#LTN_{coll}#GT","#LTp_{T}^{2}#GT","",24,20);
			htmp->GetYaxis()->SetTitleOffset(0.8);
			htmp->GetXaxis()->SetTitleOffset(1.0);

			if ( iarm==0 ){
				Ncoll_HeAu_S_sys->Draw("2");
				Ncoll_pAu_S_sys->Draw("2");
				Ncoll_pAl_S_sys->Draw("2");
				Ncoll_HeAu_S->Draw("P");
				Ncoll_pAu_S->Draw("P");
				Ncoll_pAl_S->Draw("P");
			}else{
				Ncoll_HeAu_N_sys->Draw("2");
				Ncoll_pAu_N_sys->Draw("2");
				Ncoll_pAl_N_sys->Draw("2");
				Ncoll_HeAu_N->Draw("P");
				Ncoll_pAu_N->Draw("P");
				Ncoll_pAl_N->Draw("P");
			}

			{
				TLegend *leg = new TLegend(0.2,0.68,0.6,0.93);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				le = leg->AddEntry("","Inclusive J/#psi","h");
				le = leg->AddEntry("","p_{T}<5 GeV/c","h");
				if ( iarm==0 ){
					le = leg->AddEntry("","-2.2<y<-1.2, Au/Al-going","h");
				}else{
					le = leg->AddEntry("","1.2<y<2.2, p/^{3}He-going","h");
				}
				le = leg->AddEntry("","#sqrt{s_{NN}}=200 GeV","h");
				leg->Draw();
			}

			{
				TLegend *leg = new TLegend(0.65,0.2,0.9,0.40);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->AddEntry(Ncoll_pAl_S,"p+Al","P");
				leg->AddEntry(Ncoll_pAu_S,"p+Au","P");
				leg->AddEntry(Ncoll_HeAu_S,"^{3}He+Au","P");
				leg->Draw();
			}
		}//iarm
	}//pt_squared

}

