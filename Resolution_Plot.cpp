#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TSystem.h"

Double_t Resolution(std::string particle = "e-", Double_t energy = 1.0);

void Resolution_Plot() {

  gSystem->Load("Resolution_cpp.so");

    Double_t g4_beamE[7] = {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0};
	Double_t g4_beamE_piplus[14] = {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70., 80., 90., 100.};
	Double_t g4_res_elec_5deg[7] = {0.0889507, 0.0633381, 0.0514269, 0.0405578, 0.0289683, 0.0206064, 0.0170256};
	Double_t g4_res_piplus_5deg[14] = {0.25023, 0.213926, 0.191676, 0.161576, 0.117062, 0.0926133, 0.0813332, 0.0739425, 0.0688675, 0.0639871, 0.0601296,
	0.0598187, 0.0577296, 0.0575132};

	Double_t g4_beamE_piplus_newbirk[13] = {1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70., 80., 90., 100.};
	Double_t g4_res_piplus_5deg_newbirk[13] = {0.248197, 0.205701, 0.15331, 0.110433, 0.0870408, 0.0753077, 0.0679485, 0.0630292, 0.0602589, 0.057515, 0.0557724,
	0.0540185, 0.052538};

	const Int_t num_QGSP = 13;
	Double_t g4_beamE_QGSP[num_QGSP] = {1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
	Double_t g4_res_piplus_5deg_QGSP[num_QGSP] = {0.258843, 0.22203, 0.171962, 0.127462, 0.107721, 0.0966723, 0.089173, 0.0845673, 0.0817124, 0.0788665, 0.0763072, 0.0761002, 0.0727263};
	Double_t g4_res_piplus_5deg_QGSP_TC[num_QGSP] = {0.258816, 0.221997, 0.169287, 0.125666, 0.103886, 0.0909401, 0.0840076, 0.0780045, 0.0744322, 0.0707074, 0.0701562, 0.0689857, 0.0651987};
	Double_t g4_res_piplus_5deg_QGSP_5[num_QGSP] = {0.265107, 0.229308, 0.197124, 0.147366, 0.128219, 0.117786, 0.112769, 0.108098, 0.105276, 0.09942, 0.098145, 0.098001, 0.09346}; // birk = 0.5
	Double_t g4_res_piplus_0deg_QGSP[num_QGSP] = {0.29084, 0.252189, 0.189356, 0.135934, 0.10828, 0.098379, 0.091766, 0.0849641, 0.0847324, 0.078582, 0.0767795, 0.0761457, 0.0732845}; //zdg
        Double_t g4_res_piplus_5deg_QGSP_zji[num_QGSP] = {};
        for(int ig=0; ig<num_QGSP; ig++)
          g4_res_piplus_5deg_QGSP_zji[ig] = Resolution("pi+", g4_beamE_QGSP[ig]);

	const Int_t num_CALICE = 7;
	Double_t g4_beamE_CALICE[num_CALICE] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0};
	Double_t g4_res_piplus_CALICE[num_CALICE] = {.1725, .115, .105, .095, .09, .085, .08}; 

	Double_t g4_res_elec_5deg_QGSP[7] = {0.0880496, 0.0640701, 0.0514047, 0.0403201, 0.0283405, 0.0207327, 0.0169998};

	double beam_E[10] = {3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20, 25, 32};
	double x_err[10] = {0,0,0,0,0,0,0,0,0};
	
	double factor[10] = {1.8, 2.0, 2.0, 1.8, 1.6, 1.6, 1.6, 1.2, 1.2, 1.2};
	double factor_err[10] = {0,0,0,0,0,0,0,0,0,0};
	
	double mean[10] = {487.1, 692.3, 1029.6, 1310.6, 1647.4, 1919.8, 2158.2 + 402.5/*2542.1*/, 2709.6 + 402.5 , 3364.5 + 402.5, 4413.8 + 402.5};
	double mean_err[10] = {0,0,0,0,0,0,0,0,0};
	
	//constant prop factor
	double mean_constant_factor[10] = {470.211, 621.5, 977.6, 1282.3, 1652.7, 1956.9, 2222.56 + 402.5, 2825.3 + 402.5, 3570.22 + 402.5, 4617.58 + 402.5};
	double mean_err_constantfactor[10] = {0,0,0,0,0,0,0,0,0};

	double sigma_constant_factor[10]= {154.9, 215.0, 231.8, 277.2, 345.6, 379.9, 433.185, 546.691, 642.301, 757.404};
	double sigma_err_constant_factor[10] = {0,0,0,0,0,0,0,0,0,0};
	double reso_err_constant_factor[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	
	//sigmas
	double sigma[10] = {165.0, 238.1, 247.6, 265.2, 329.5, 349.0, 402.5, 439.4, 536.3, 623.7};
	double sigma_err[10] = {0,0,0,0,0,0,0,0,0,0};
	double reso_err[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double ecal_reso_err[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	
	double ecal_mean[10] = {478.1, 622.4, 895.7, 1194.3, 1479.5, 1798.2, 2344.0, 2871.1, 3555.6, 4544.0};
	double ecal_mean_err[10] = {0,0,0,0,0,0,0,0,0};
	
	double ecal_sigma[10] = {37.0, 41.7, 51.5, 65.4, 68.5, 81.5, 104.5, 118.1, 137.6, 144.3};
	
	double beam_spread[10] = {.027,.023,.023,.023,.023,.023,.023,.023,.023,.023};


	
	double reso[10];
	double reso_constant_factor[10];
	double ecal_reso[10];
	for(int i=0;i<10;i++) {
		reso[i] = TMath::Sqrt(sigma[i]*sigma[i]/(mean[i]*mean[i]) - beam_spread[i]*beam_spread[i]);
		reso_constant_factor[i] = TMath::Sqrt(sigma_constant_factor[i]*sigma_constant_factor[i]/(mean_constant_factor[i]*mean_constant_factor[i])- beam_spread[i]*beam_spread[i]);
		//ecal
		ecal_reso[i] = TMath::Sqrt(ecal_sigma[i]*ecal_sigma[i]/(ecal_mean[i]*ecal_mean[i]) - beam_spread[i]*beam_spread[i]);
	}

	gStyle->SetOptFit(0000);


	TGraphErrors *gECalResolution = new TGraphErrors(10, beam_E, ecal_reso, 0, 0);
	gECalResolution->SetMarkerStyle(29);
	gECalResolution->SetMarkerColor(kGreen+2);
	gECalResolution->SetMarkerSize(2.4);
	gECalResolution->SetLineStyle(2);
	gECalResolution->SetLineColor(kGreen+2);
	
	TF1 *fresofitecal = new TF1("fresofitecal", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 1, 32);
	fresofitecal->SetLineStyle(2);
	fresofitecal->SetLineColor(kGreen+2);

	
	gECalResolution->Fit(fresofitecal, "R");
	gECalResolution->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gECalResolution->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gECalResolution->GetYaxis()->SetTitleOffset(1.4);

	
	TGraphErrors *gHCalResolution = new TGraphErrors(10, beam_E, reso, x_err, reso_err);
	gHCalResolution->SetMarkerStyle(33);
	gHCalResolution->SetMarkerColor(kRed);
	gHCalResolution->SetMarkerSize(2.4);
	gHCalResolution->SetLineStyle(2);
	gHCalResolution->SetLineColor(kRed);

	
	TGraphErrors *gSumFactor = new TGraphErrors(10, beam_E, factor, 0, factor_err);
	gSumFactor->SetMarkerStyle(20);
	gSumFactor->SetMarkerColor(kMagenta+2);
	gSumFactor->SetMarkerSize(2.1);

	TF1 *fresofit = new TF1("fresofit", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 1, 32);
	fresofit->SetLineStyle(2);
	fresofit->SetLineColor(kRed);



    TGraphErrors* gr_g4elec_res = new TGraphErrors(7, g4_beamE, g4_res_elec_5deg, 0, 0);
	gr_g4elec_res->SetMarkerStyle(kOpenStar);
	gr_g4elec_res->SetMarkerColor(kBlue);
	gr_g4elec_res->SetMarkerSize(2);
	gr_g4elec_res->SetLineStyle(2);
	gr_g4elec_res->SetLineColor(kBlue);
	
	TF1* f_g4elec_fit = new TF1("f_g4elec_fit", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 1, 32);
	f_g4elec_fit->SetParameter(0, .038);
	f_g4elec_fit->SetParameter(1, .43);
	f_g4elec_fit->SetLineStyle(2);
	f_g4elec_fit->SetLineColor(kBlue);

	gr_g4elec_res->Fit(f_g4elec_fit, "R");
	gr_g4elec_res->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4elec_res->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4elec_res->GetYaxis()->SetTitleOffset(1.4);

	TGraphErrors* gr_g4piplus_res = new TGraphErrors(14, g4_beamE_piplus, g4_res_piplus_5deg, 0, 0);
	gr_g4piplus_res->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res->SetMarkerColor(kBlue);
	gr_g4piplus_res->SetMarkerSize(2);
	gr_g4piplus_res->SetLineStyle(2);
	gr_g4piplus_res->SetLineColor(kBlue);
	
	TF1* f_g4piplus_fit = new TF1("f_g4piplus_fit", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	f_g4piplus_fit->SetParameter(0, .038);
	f_g4piplus_fit->SetParameter(1, .43);
	f_g4piplus_fit->SetLineStyle(2);
	f_g4piplus_fit->SetLineColor(kBlue);

	gr_g4piplus_res->Fit(f_g4piplus_fit, "R");
	gr_g4piplus_res->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_newbirk = new TGraphErrors(13, g4_beamE_piplus_newbirk, g4_res_piplus_5deg_newbirk, 0, 0);
	gr_g4piplus_res_newbirk->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res_newbirk->SetMarkerColor(kRed);
	gr_g4piplus_res_newbirk->SetMarkerSize(2);
	gr_g4piplus_res_newbirk->SetLineStyle(2);
	gr_g4piplus_res_newbirk->SetLineColor(kRed);
	
	TF1* f_g4piplus_fit_newbirk = new TF1("f_g4piplus_fit_newbirk", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	f_g4piplus_fit_newbirk->SetParameter(0, .038);
	f_g4piplus_fit_newbirk->SetParameter(1, .43);
	f_g4piplus_fit_newbirk->SetLineStyle(2);
	f_g4piplus_fit_newbirk->SetLineColor(kRed);

	gr_g4piplus_res_newbirk->Fit(f_g4piplus_fit_newbirk, "R");
	gr_g4piplus_res_newbirk->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_newbirk->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_newbirk->GetYaxis()->SetTitleOffset(1.4);

	TGraphErrors* gr_g4piplus_res_QGSP = new TGraphErrors(num_QGSP, g4_beamE_QGSP, g4_res_piplus_5deg_QGSP, 0, 0);
	gr_g4piplus_res_QGSP->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res_QGSP->SetMarkerColor(kBlack);
	gr_g4piplus_res_QGSP->SetMarkerSize(2);
	gr_g4piplus_res_QGSP->SetLineStyle(2);
	gr_g4piplus_res_QGSP->SetLineColor(kBlack);
	
	TF1* f_g4piplus_fit_QGSP = new TF1("f_g4piplus_fit_newbirk", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	f_g4piplus_fit_QGSP->SetLineStyle(2);
	f_g4piplus_fit_QGSP->SetLineColor(kBlack);

	gr_g4piplus_res_QGSP->Fit(f_g4piplus_fit_QGSP, "R");
	gr_g4piplus_res_QGSP->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_QGSP->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_QGSP->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_QGSP_5 = new TGraphErrors(num_QGSP, g4_beamE_QGSP, g4_res_piplus_5deg_QGSP_5, 0, 0);
	gr_g4piplus_res_QGSP_5->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res_QGSP_5->SetMarkerColor(kCyan);
	gr_g4piplus_res_QGSP_5->SetMarkerSize(2);
	gr_g4piplus_res_QGSP_5->SetLineStyle(2);
	gr_g4piplus_res_QGSP_5->SetLineColor(kCyan);
	
	TF1* f_g4piplus_fit_QGSP_5 = new TF1("f_g4piplus_fit_QGSP_5", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);

	f_g4piplus_fit_QGSP_5->SetLineStyle(2);
	f_g4piplus_fit_QGSP_5->SetLineColor(kCyan);

	gr_g4piplus_res_QGSP_5->Fit(f_g4piplus_fit_QGSP_5, "R");
	gr_g4piplus_res_QGSP_5->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_QGSP_5->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_QGSP_5->GetYaxis()->SetTitleOffset(1.4);

	// TGraphErrors* gr_g4piplus_res_QGSP_zdg = new TGraphErrors(num_QGSP, g4_beamE_QGSP, g4_res_piplus_0deg_QGSP, 0, 0);
	// gr_g4piplus_res_QGSP_zdg->SetMarkerStyle(kOpenDiamond);
	// gr_g4piplus_res_QGSP_zdg->SetMarkerColor(kMagenta);
	// gr_g4piplus_res_QGSP_zdg->SetMarkerSize(2);
	// gr_g4piplus_res_QGSP_zdg->SetLineStyle(2);
	// gr_g4piplus_res_QGSP_zdg->SetLineColor(kMagenta);
	
	// TF1* f_g4piplus_fit_QGSP_zdg = new TF1("f_g4piplus_fit_QGSP_zdg", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	// // f_g4piplus_fit_QGSP->SetParameter(0, .038);
	// // f_g4piplus_fit_QGSP->SetParameter(1, .43);
	// f_g4piplus_fit_QGSP_zdg->SetLineStyle(2);
	// f_g4piplus_fit_QGSP_zdg->SetLineColor(kMagenta);

	// gr_g4piplus_res_QGSP_zdg->Fit(f_g4piplus_fit_QGSP_zdg, "R");
	// gr_g4piplus_res_QGSP_zdg->GetXaxis()->SetTitle("Beam Energy (GeV)");
	// gr_g4piplus_res_QGSP_zdg->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	// gr_g4piplus_res_QGSP_zdg->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_QGSP_TC = new TGraphErrors(num_QGSP, g4_beamE_QGSP, g4_res_piplus_5deg_QGSP_TC, 0, 0);
	gr_g4piplus_res_QGSP_TC->SetMarkerStyle(kOpenDiamond);
	gr_g4piplus_res_QGSP_TC->SetMarkerColor(kMagenta);
	gr_g4piplus_res_QGSP_TC->SetMarkerSize(2);
	gr_g4piplus_res_QGSP_TC->SetLineStyle(2);
	gr_g4piplus_res_QGSP_TC->SetLineColor(kMagenta);

	TF1* f_g4piplus_fit_QGSP_TC = new TF1("f_g4piplus_fit_QGSP_TC", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	// f_g4piplus_fit_QGSP->SetParameter(0, .038);
	// f_g4piplus_fit_QGSP->SetParameter(1, .43);
	f_g4piplus_fit_QGSP_TC->SetLineStyle(2);
	f_g4piplus_fit_QGSP_TC->SetLineColor(kMagenta);

	gr_g4piplus_res_QGSP_TC->Fit(f_g4piplus_fit_QGSP_TC, "R");
	gr_g4piplus_res_QGSP_TC->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_QGSP_TC->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_QGSP_TC->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_QGSP_zji = new TGraphErrors(num_QGSP, g4_beamE_QGSP, g4_res_piplus_5deg_QGSP_zji, 0, 0);
	gr_g4piplus_res_QGSP_zji->SetMarkerStyle(kOpenCircle);
	gr_g4piplus_res_QGSP_zji->SetMarkerColor(kRed);
	gr_g4piplus_res_QGSP_zji->SetMarkerSize(2);
	gr_g4piplus_res_QGSP_zji->SetLineStyle(2);
	gr_g4piplus_res_QGSP_zji->SetLineColor(kMagenta);

	TF1* f_g4piplus_fit_QGSP_zji = new TF1("f_g4piplus_fit_QGSP_zji", "sqrt(([1]/sqrt(x))**2 + [0]**2)", 10, 102);
	// f_g4piplus_fit_QGSP->SetParameter(0, .038);
	// f_g4piplus_fit_QGSP->SetParameter(1, .43);
	f_g4piplus_fit_QGSP_zji->SetLineStyle(2);
	f_g4piplus_fit_QGSP_zji->SetLineColor(kMagenta);

	gr_g4piplus_res_QGSP_zji->Fit(f_g4piplus_fit_QGSP_zji, "R");
	gr_g4piplus_res_QGSP_zji->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_QGSP_zji->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_QGSP_zji->GetYaxis()->SetTitleOffset(1.4);


	TGraphErrors* gr_g4piplus_res_CALICE = new TGraphErrors(num_CALICE, g4_beamE_CALICE, g4_res_piplus_CALICE, 0, 0);
	gr_g4piplus_res_CALICE->SetMarkerStyle(kFullCircle);
	gr_g4piplus_res_CALICE->SetMarkerColor(kBlack);
	gr_g4piplus_res_CALICE->SetMarkerSize(1.5);
	gr_g4piplus_res_CALICE->SetLineStyle(2);
	gr_g4piplus_res_CALICE->SetLineColor(kBlack);
	
	TF1* f_g4piplus_fit_CALICE = new TF1("f_g4piplus_fit_CALICE", "sqrt(([1]/sqrt(x))**2 + [0]**2 + (.18/x)**2)", 10, 90);
	// f_g4piplus_fit_QGSP->SetParameter(0, .038);
	// f_g4piplus_fit_QGSP->SetParameter(1, .43);
	f_g4piplus_fit_CALICE->SetLineStyle(2);
	f_g4piplus_fit_CALICE->SetLineColor(kBlack);

	gr_g4piplus_res_CALICE->Fit(f_g4piplus_fit_CALICE, "R");
	gr_g4piplus_res_CALICE->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gr_g4piplus_res_CALICE->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gr_g4piplus_res_CALICE->GetYaxis()->SetTitleOffset(1.4);

	
	// TF1* f_g4elec_fit_QGSP = new TF1("f_g4piplus_fit_newbirk", "[1]/sqrt(x) + [0] + .18/x", .5, 32);
	// // f_g4piplus_fit_QGSP->SetParameter(0, .038);
	// // f_g4piplus_fit_QGSP->SetParameter(1, .43);
	// f_g4elec_fit_QGSP->SetLineStyle(2);
	// f_g4elec_fit_QGSP->SetLineColor(kBlack);

	// gr_g4elec_res_QGSP->Fit(f_g4elec_fit_QGSP, "R");
	// gr_g4elec_res_QGSP->GetXaxis()->SetTitle("Beam Energy (GeV)");
	// gr_g4elec_res_QGSP->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	// gr_g4elec_res_QGSP->GetYaxis()->SetTitleOffset(1.4);

	
 		
	//res from simulation
	TF1 *CALICE_data = new TF1("CALICE_data", "sqrt( ([1]/sqrt(x))**2 + [0]**2 + (.18/x)**2)", 8, 90);
	CALICE_data->SetParameter(0, .04);
	CALICE_data->SetParameter(1, .518);
	CALICE_data->SetLineColor(kBlack);
	CALICE_data->SetLineWidth(2);
	
	TCanvas *cHCalResolution = new TCanvas("cHCalResolution", "HCal Resolution", 900, 700);
	cHCalResolution->cd();
	// gHCalResolution->Fit(fresofit, "R");

	//gECalResolution->Draw("AP");
	gr_g4piplus_res->SetTitle("");

	gr_g4piplus_res->Draw("AP");
	gr_g4piplus_res->GetXaxis()->SetRange(6, 105);
	
    //gr_g4elec_res->Draw("sameP");
	// gHCalResolution->Draw("sameP");
	gr_g4piplus_res_newbirk->Draw("sameP");
	gr_g4piplus_res_QGSP->Draw("sameP");
	// gr_g4piplus_res_QGSP_5->Draw("sameP");
	//gr_g4elec_res_QGSP->Draw("sameP");
	gr_g4piplus_res_QGSP_TC->Draw("sameP");
	gr_g4piplus_res_QGSP_zji->Draw("sameP");
	// gr_g4piplus_res_CALICE->Draw("sameP");
	CALICE_data->Draw("same");


	gECalResolution->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gECalResolution->GetXaxis()->SetLimits(0.0, 102.0);
	gECalResolution->GetYaxis()->SetTitle("Resolution  #frac{#DeltaE}{E}");
	gECalResolution->GetYaxis()->SetTitleOffset(1.4);
	gECalResolution->GetYaxis()->SetRangeUser(0.0, 0.5);
	gECalResolution->SetTitle("");
	
	// MC_hadron_resolution->Draw("sameL");
	
	TLegend *leg = new TLegend(.3, .4, .85, .85);
	leg->SetHeader("Resolution of ECAL + HCAL");
	// leg->AddEntry(gHCalResolution, "hadron,  #frac{#DeltaE}{E} = #frac{0.58 #pm 0.04}{#sqrt{E}} + 0.015 #pm 0.015", "LP");

	// TString piplus_text_QGSP_5;
	// piplus_text_QGSP_5.Form("G4 #pi^{+}, QGSP_BERT kB=.5, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit_QGSP_5->GetParameter(1), f_g4piplus_fit_QGSP_5->GetParameter(0));
	// leg->AddEntry(gr_g4piplus_res_QGSP_5, piplus_text_QGSP_5, "LP");

	TString piplus_text_QGSP_TC;
	piplus_text_QGSP_TC.Form("G4 #pi^{+}, QGSP_BERT kB=.126 TC, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit_QGSP_TC->GetParameter(1), f_g4piplus_fit_QGSP_TC->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res_QGSP_TC, piplus_text_QGSP_TC, "LP");

	TString piplus_text_QGSP_zji;
	piplus_text_QGSP_zji.Form("G4 #pi^{+}, QGSP_BERT kB=.126 zji, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit_QGSP_zji->GetParameter(1), f_g4piplus_fit_QGSP_zji->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res_QGSP_zji, piplus_text_QGSP_zji, "LP");

	TString piplus_text_QGSP;
	piplus_text_QGSP.Form("G4 #pi^{+}, QGSP_BERT kB=.126, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit_QGSP->GetParameter(1), f_g4piplus_fit_QGSP->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res_QGSP, piplus_text_QGSP, "LP");

	TString piplus_text;
	piplus_text.Form("G4 #pi^{+}, FTFP_BERT_HP kB=.126, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit->GetParameter(1), f_g4piplus_fit->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res, piplus_text, "LP");

	TString piplus_text_newbirk;
	piplus_text_newbirk.Form("G4 #pi^{+}, FTFP_BERT_HP kB = .078, #frac{#DeltaE}{E} =  #frac{%0.3f}{#sqrt{E}} #oplus %0.3f", f_g4piplus_fit_newbirk->GetParameter(1), f_g4piplus_fit_newbirk->GetParameter(0));
	leg->AddEntry(gr_g4piplus_res_newbirk, piplus_text_newbirk, "LP");
	
	TString CALICE;
	CALICE.Form("CALICE QGSP_BERT, #frac{#DeltaE}{E} = #frac{0.518}{#sqrt{E}} #oplus 0.04 #oplus #frac{0.18}{E}");
	leg->AddEntry(CALICE_data, CALICE, "LP");

	// leg->AddEntry(gECalResolution, "electron, FNAL, #frac{#DeltaE}{E} =  #frac{0.11 #pm 0.007}{#sqrt{E}} + 0.007 #pm 0.002", "LP");

	// TString electron_text;
	// electron_text.Form("G4 electron, FTFP_BERT_HP kB=.126, #frac{#DeltaE}{E} =  #frac{%0.2f}{#sqrt{E}} #oplus %0.3f ", f_g4elec_fit->GetParameter(1), f_g4elec_fit->GetParameter(0));
	// leg->AddEntry(gr_g4elec_res, electron_text, "LP");

	// TString electron_text_QGSP;
	// electron_text_QGSP.Form("G4 electron, QGSP_BERT, #frac{#DeltaE}{E} =  #frac{%0.2f}{#sqrt{E}} + %0.3f ", f_g4elec_fit_QGSP->GetParameter(1), f_g4elec_fit_QGSP->GetParameter(0));
	// leg->AddEntry(gr_g4elec_res_QGSP, electron_text_QGSP, "LP");

	

	leg->SetTextSize(.028);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->Draw();
	
	cHCalResolution->SetRightMargin(.02);
	gStyle->SetStatX(.8);
	gStyle->SetStatY(.85);
	gStyle->SetStatW(.15);
	gStyle->SetStatH(.13);
	cHCalResolution->Update();

}
