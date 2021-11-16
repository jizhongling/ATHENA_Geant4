// used Masha's basic code
// use two passes, first pass to find optimal weight for eca; which is enrgy dependent
#include<fstream>
#include<iostream>
#include<cstdlib>
#include <string>
#include <sstream>
#include "TH2.h"

using namespace std;
//hcal we want towers 22 32 23 33
/*
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^MAP of TOWERS in HCAL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                                    [00][10][20][30][40][50]
                                    [01][11][21][31][41][51]
                                    [02][12][22][32][42][52]
                                    [03][13][23][33][43][53]
                                    [04][14][24][34][44][54]
                                    [05][15][25][35][45][55]
 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^MAP of TOWERS in HCAL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*/




double threshold = 6.1;

double reco(double E0)
{
double E = E0;
if(E0 < threshold) return 0;
else{ 
double a = 0.1;
double b = 0.0015;
double sigma =  E* TMath::Sqrt(a*a/E +b*b) ;
double random = gRandom->Gaus(E*0.03,sigma); 
E = random;
return E;
}
}

const int num_of_tiles=52;//num of tiles in hadronic calorimeter //36
const int num_of_towers_ecal=576;//num of towers in ecal
const int num_of_towers=36;//num of towers in had cal //36
const int nentries=10000;//number of events;
const float w0 = 0.7; // initial weight for ecal
const float wopt = 1.4; // optimal weight for ecal after first pass

double filename[6] ={0.5,1,2,5,10,20};
const char *energylevel[9] = {"MeV","MeV","MeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV"};
double minenergy[9]= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};//{80.5,160.5,400.5,800.5,1600.5, 4000.5, 8000.5, 16000.5, 30000.5};
double maxenergy[9]= {50.5,100.5,200.5,500.5,1000.5, 2000.5, 3000.5,4000.5,5000.5};

float energyname[9]= {0.5,1,2,5,10,20,30,40,50};

int weight_layers_reco(){

int cen = 5;
    ifstream myfile;
    char buffer[100];

	char fname_log[200];
	sprintf(fname_log,"20mm-mixture-52layer-pi+/20.log");
     myfile.open(fname_log);
    size_t pos;
    
    if(myfile.fail())
    {
        cout << "Input file opening failed.\n";
        exit(1);
    }
    //works up to here
    
 	char fname_out[200];
	sprintf(fname_out,"energy%3.1f_pion.mixture.reco.root", energyname[cen]);
	TFile fout(fname_out,"RECREATE");
    //********************Histograms*************************************
   
    TH1D* emcal_tot_hist=new TH1D("emcal_tot_hist","emcal_tot_hist",100 ,0.5,maxenergy[cen]);
    TH1D* hcal1_tot_hist=new TH1D("hcal1_tot_hist","hcal1_tot_hist",100,0.5,maxenergy[cen]);
         TH1D* emcal_tower_hist =new TH1D("emcal_tower_hist","emcal_tower_hist",100,0.5,maxenergy[cen]);
         TH2F *Hist_ec_hc1 = new TH2F("Hist_ec_hc1","Hist_ec_hc1",200,0.5,0.25+maxenergy[cen]/2,200,0.5,0.25+maxenergy[cen]/2); 

    char  fname[200];

        //histos for different weights in ecal
    TH1D* EW[20];
    char buffer1[100];
    char buffer2[100];

    for(int i=0; i<20; i++) {
        sprintf(buffer1, "EW%d", i);
        sprintf(buffer2, "Weighted FCS weight = %d + 5 ", i);
        EW[i] = new TH1D(buffer1, buffer2, 200, minenergy[cen], maxenergy[cen]);
    }
    
    float w[20];  //weights for ecal
    for(int i=0;i<20;i++){
        w[i] = w0 + 0.05*float(i+1);
    }

    //********************Histograms DONE*************************************
    
    //works up to here
    string search="EM: ";//ti find the right data in the text file for em
    string search_had="HAD: ";//to find the right data in the text file for had
    
    string wordline;//to read text file line by line
    
    double had1_en=0;//had1 energy
    
    int count_em=0;//em tower identifier
    int count_had=0;//had tower identifier
    
    bool em=false;//flag to read em data
    bool had=false;//flag to read had data
    int num_of_events=0;//counter for event number
    int tot_had_count=0;//counter for total energy
    
    
    //works up to here
    double* edep_em_tot= new double [nentries]; //total energy per event in ECAL
    double* edep_had_tot=new double [nentries];//total energy for had1 in all towers for each event
    TH1D* edep_ec_hc1=new TH1D("edep_ec_hc1","edep_ec_hc1",600,minenergy[cen], maxenergy[cen]);

    for(int i=0;i<nentries;i++){
                edep_em_tot[num_of_events-1]=0;
                edep_had_tot[num_of_events-1]=0;
    }    
    
    double*** edep_had_tile=new double**[nentries];//energy in each tower in had1 in each tile
        double** edep_em_tower=new double*[nentries];//energy in each sector in ECAL

    for(int i=0;i<nentries;i++){
        edep_had_tile[i]=new double*[num_of_towers];
        edep_em_tower[i]=new double[576];

        for(int j=0;j<num_of_towers;j++){
            edep_had_tile[i][j]=new double[num_of_tiles];
        }
	for(int j=0;j<576;j++) {
		edep_em_tower[i][j] = 0;
	}

    }
    
    int event_count=0;
    //********************Get DATA*************************************
       
        while(getline(myfile, wordline)){
            
            if(em==true) {
                double em_en=0.0;//em energy
                stringstream geek(wordline);
                geek>>em_en;
                edep_em_tot[num_of_events-1]+= reco(em_en) ;
		edep_em_tower[num_of_events-1][count_em]+=reco(em_en) ;

                if(count_em==(num_of_towers_ecal-1)) em=false, count_em=0;
                if(em==true) count_em++;
            }
            if(wordline==search) em=true, num_of_events++;
            
            if(had==true) {
                
                //tiles for had1
                if(tot_had_count==1){
                    stringstream ss;
                    ss<<wordline;
                    string temp;
                    double found;
                    double h2_ti_en=0.0;
                    int tile_num=0;
                    while(!ss.eof()){
                        ss >> temp;
                        if(stringstream(temp) >> found) h2_ti_en=found;
                        temp="";
                        edep_had_tile[num_of_events-1][count_had][tile_num]+=found;
                        tile_num++;
                        
                    }
                    
                    tot_had_count++;
                }
                
                //total for had 1
                if(tot_had_count==0){
                    stringstream geek(wordline);
                    geek>>had1_en;
                    edep_had_tot[num_of_events-1]+=had1_en;
                    tot_had_count++;
                }
                if(tot_had_count>1) tot_had_count=0, count_had++;
                if(count_had== num_of_towers) had=false, count_had=0;
               
            }
            
            
            if(wordline==search_had) had=true;
            
        }
    //------------------------------------got the data------------------------------------
    
    cout<<"total events: "<<num_of_events<<endl;  
    //fill histograms
    
    for(int i=0;i<nentries;i++){
        //for total energy gotten from the file
        emcal_tot_hist->Fill(edep_em_tot[i]);
        hcal1_tot_hist->Fill(edep_had_tot[i]);
	for(int j =0;j<576;j++) {
		emcal_tower_hist->Fill(edep_em_tower[i][j]);
	}
	Hist_ec_hc1->Fill(edep_em_tot[i] , edep_had_tot[i] , 1);


    }

    for(int i=0;i<nentries;i++){
        edep_ec_hc1->Fill(edep_em_tot[i]/wopt + edep_had_tot[i]);
 
       for(int j=0;j<20;j++){   //filling total emnergy with diff. weight in ecal
            EW[j]->Fill(edep_em_tot[i]/w[j]+edep_had_tot[i]);
        }
    }
     


//-------------------DRAW---------------

   gStyle->SetOptFit(0001);
    gStyle->SetPalette(1);
    
    
//final resolution numbers
    float mean = 0.;
    float sigma = 0.;
    float resolution = 0.;
    
    TCanvas* c = new TCanvas("With Optimal Weight");
    edep_ec_hc1->Fit("gaus");
    mean = edep_ec_hc1->GetFunction("gaus")->GetParameter(1);
    sigma = edep_ec_hc1->GetFunction("gaus")->GetParameter(2);
    if (mean != 0) resolution = sigma/mean; 
    edep_ec_hc1->Draw();
    c->Update();
    	char fname3[200];
	sprintf(fname3,"plots-reco/e.%3.1fgev.pdf", energyname[cen]);
    c->Print(fname3);
    
    cout<<"******************* Energy Resolution: "<< resolution <<endl;



    float P1[64];
    float P2[64];
    float Rw[64];
    
    for(int i=0;i<20;i++){
        EW[i]->Fit("gaus");
        P1[i] = EW[i]->GetFunction("gaus")->GetParameter(1);
        P2[i] = EW[i]->GetFunction("gaus")->GetParameter(2);
        if(P1[i] != 0) Rw[i]=P2[i]/P1[i];
    }
    
    TGraph *BestWeight = new TGraph(20, w, Rw);
    BestWeight->SetMarkerStyle(20);
    BestWeight->SetMarkerColor(kBlue+2);
    BestWeight->SetMarkerSize(2.4);
    TCanvas *c6 = new TCanvas();
       c6->Divide(4,5);
       for(int i=0;i<20;i++){
           c6->cd(i+1);
           EW[i]->SetLineColor(kRed +2);
           EW[i]->Draw();
       }
    c6->Update();//this draws the best weight.
    	sprintf(fname3,"plots-reco/e.%3.1fgev.2.pdf", energyname[cen]);
    c6->Print(fname3);
  
    TCanvas *c7 = new TCanvas();
    BestWeight->Draw("AP"); //this draws the ratio of best weight
    c7->Update();
    	sprintf(fname3,"plots-reco/e.%3.1fgev.3.pdf", energyname[cen]);
    c7->Print(fname3);




//-------------------
	fout.Write();


    delete[] edep_em_tot;
    delete[] edep_had_tot;
    for(int i=0;i<nentries;i++){
	delete[] edep_em_tower[i];

        for(int j=0;j<num_of_towers;j++){
            delete[] edep_had_tile[i][j];
        }
        delete[] edep_had_tile[i] ;
    }
    delete[] edep_em_tower;
    delete[] edep_had_tile;
    
    return 0;
    
}
