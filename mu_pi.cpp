/* 
Analysis code for muon/pion separation in hadron endcap
Available angles: 5 deg, 20 deg
20 deg only has energies up to 20 GeV

Sections: 
1 - ECal
2 - HCal layers 1-9
3 - HCal layers 10-18
4 - HCal layers 19-51
*/


#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include <map>
#include <vector>
#include <tuple>

// Max energies for histogram upper bounds
std::map<Double_t, Double_t> muon_max_energy = {
    {1.0, 50},
    {2.0, 50},
    {5.0, 50},
    {10.0, 50},
    {20.0, 50},
    {30.0, 50},
    {40.0, 50},
    {50.0, 50},
    {60., 50},
    {70., 50},
    {80., 50},
    {90., 50},
    {100., 50}};
std::map<Double_t, Double_t> pion_max_energy = {
    {1.0, 50}, 
    {2.0, 50}, 
    {5.0, 50}, 
    {10.0, 50}, 
    {20.0, 50}, 
    {30.0, 50}, 
    {40.0, 50}, 
    {50.0, 50}, 
    {60., 50}, 
    {70., 50}, 
    {80., 50}, 
    {90., 50}, 
    {100., 50}};

std::map<Double_t, Double_t> ener_to_mult = {
    {1.0, 3.}, 
    {2.0, 4.}, 
    {5.0, 5.}, 
    {10.0, 7.}, 
    {20.0, 8.}, 
    {30.0, 9.}, 
    {40.0, 11.}, 
    {50.0, 12.}, 
    {60., 14.}, 
    {70., 15.}, 
    {80., 16.}, 
    {90., 17.}, 
    {100., 18.}};

std::map<Double_t, Double_t> ener_to_ecal_cut = {{1.0, 9.8}, {2.0, 15.}, {5.0, 20.}, {10.0, 22.}, {20.0, 25.}};
std::map<Double_t, Double_t> ener_to_hcal_cut = {{1.0, 0.}, {2.0, 10.}, {5.0, 12.}, {10.0, 14.}, {20.0, 16.}};

const Int_t num_energies = 13;
Double_t energies[13] = {1., 2., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
const Int_t num_total_layers = 51;
const Int_t num_towers = 36;
const Int_t num_blocks = 64;

void mu_pi(Double_t energy = 1.0, Int_t angle = 20)
{
    // Opening muon file
    TString muon_file_name;
    muon_file_name.Form(
        "mu-_QGSP_Default_%ddeg/mu-_%0.0fGeV_%ddeg.root", 
        angle, 
        energy, 
        angle);
    std::cout << "Opening " << muon_file_name << std::endl;
    TFile *muon_data_file = new TFile(muon_file_name);

    // Opening pion file
    TString pion_file_name;
    pion_file_name.Form(
        "pi+_QGSP_Default_%ddeg/pi+_%0.0fGeV_%ddeg.root", 
        angle, 
        energy, 
        angle);
    std::cout << "Opening " << pion_file_name << std::endl;
    TFile *pion_data_file = new TFile(pion_file_name);

    // Getting the trees from the data file

    // Tree that holds total energies in ECal and HCal
    TTree *muon_total_tree = (TTree *)muon_data_file->Get("EdepTotal");
    // Tree that holds individual block information for ECal
    TTree *muon_ecal_block_tree = (TTree *)muon_data_file->Get("ECalBlocks");
    // Tree that holds individual tower information for HCal
    TTree *muon_hcal_tile_tree = (TTree *)muon_data_file->Get("HCalLayers");

    // Setting up the branch readouts
    // ECal readout variables
    Double_t muon_ecal_block_edep;
    Int_t muon_ecal_block_eventID, muon_ecal_XblockID, muon_ecal_YblockID;
    Double_t muon_hcal_tile_edep;
    // HCal readout variables
    Int_t muon_hcal_layerID, muon_hcal_tile_eventID, muon_hcal_XtowerID, muon_hcal_YtowerID;

    muon_ecal_block_tree->SetBranchAddress("ECal_Edep_Block", &muon_ecal_block_edep);
    muon_ecal_block_tree->SetBranchAddress("ECal_BlockXid", &muon_ecal_XblockID);
    muon_ecal_block_tree->SetBranchAddress("ECal_BlockY", &muon_ecal_YblockID);
    muon_ecal_block_tree->SetBranchAddress("eventID", &muon_ecal_block_eventID);

    muon_hcal_tile_tree->SetBranchAddress("HCal_Edep_Tile", &muon_hcal_tile_edep);
    muon_hcal_tile_tree->SetBranchAddress("eventID", &muon_hcal_tile_eventID);
    muon_hcal_tile_tree->SetBranchAddress("HCal_TowerXid", &muon_hcal_XtowerID);
    muon_hcal_tile_tree->SetBranchAddress("HCal_TowerYid", &muon_hcal_YtowerID);
    muon_hcal_tile_tree->SetBranchAddress("HCal_Layerid", &muon_hcal_layerID);


    // Tree that holds total energies in ECal and HCal
    TTree *pion_total_tree = (TTree *)pion_data_file->Get("EdepTotal");
    // Tree that holds individual block information for ECal
    TTree *pion_ecal_block_tree = (TTree *)pion_data_file->Get("ECalBlocks");
    // Tree that holds individual tower information for HCal
    TTree *pion_hcal_tile_tree = (TTree *)pion_data_file->Get("HCalLayers");

    // Setting up the branch readouts
    // ECal readout variables
    Double_t pion_ecal_block_edep;
    Int_t pion_ecal_block_eventID, pion_ecal_XblockID, pion_ecal_YblockID;
    Double_t pion_hcal_tile_edep;
    // HCal readout variables
    Int_t pion_hcal_layerID, pion_hcal_tile_eventID, pion_hcal_XtowerID, pion_hcal_YtowerID, pion_hcal_numhits;

    pion_ecal_block_tree->SetBranchAddress("ECal_Edep_Block", &pion_ecal_block_edep);
    pion_ecal_block_tree->SetBranchAddress("ECal_BlockXid", &pion_ecal_XblockID);
    pion_ecal_block_tree->SetBranchAddress("ECal_BlockY", &pion_ecal_YblockID);
    pion_ecal_block_tree->SetBranchAddress("eventID", &pion_ecal_block_eventID);

    pion_hcal_tile_tree->SetBranchAddress("HCal_Edep_Tile", &pion_hcal_tile_edep);
    pion_hcal_tile_tree->SetBranchAddress("eventID", &pion_hcal_tile_eventID);
    pion_hcal_tile_tree->SetBranchAddress("HCal_TowerXid", &pion_hcal_XtowerID);
    pion_hcal_tile_tree->SetBranchAddress("HCal_TowerYid", &pion_hcal_YtowerID);
    pion_hcal_tile_tree->SetBranchAddress("HCal_Layerid", &pion_hcal_layerID);
    pion_hcal_tile_tree->SetBranchAddress("HCal_NumHits_Tile", &pion_hcal_numhits);

    const Int_t muon_num_events = (Int_t)muon_total_tree->GetEntries();
    std::cout << "Number of muon events: " << muon_num_events << std::endl;

    // Vector that holds edep in each ecal block in each event
    std::vector<std::vector<Double_t>> muon_ecal_event_edep(muon_num_events);
    std::generate(
        muon_ecal_event_edep.begin(),
        muon_ecal_event_edep.end(),
        []
        {
            std::vector<Double_t> block_event_edep(num_blocks);
            std::fill(block_event_edep.begin(), block_event_edep.end(), 0.);
            return block_event_edep;
        }
    );
    // Vector that holds edep in each hcal section in each tower in each event
    // i.e. select event, then tower number, then section
    std::vector<std::vector<std::vector<Double_t>>> muon_hcal_event_edep(muon_num_events);
    std::generate(
        muon_hcal_event_edep.begin(),
        muon_hcal_event_edep.end(),
        []
        {
            std::vector<std::vector<Double_t>> tower_event_edep(num_towers);
            std::generate(
                tower_event_edep.begin(),
                tower_event_edep.end(),
                []
                {
                    std::vector<Double_t> section_event_edep(3);
                    std::fill(section_event_edep.begin(), section_event_edep.end(), 0.);
                    return section_event_edep;
                });
            return tower_event_edep;
        }
    );

    // Vector that holds edep in each detector section in each event
    std::vector<std::vector<Double_t>> muon_section_event_edep(muon_num_events);
    std::generate(
        muon_section_event_edep.begin(),
        muon_section_event_edep.end(),
        []
        {
            std::vector<Double_t> section_event_edep(4);
            std::fill(section_event_edep.begin(), section_event_edep.end(), 0.);
            return section_event_edep;
    });

    // Vector that holds number of hits in each detector section in each event
    std::vector<std::vector<Int_t>> muon_section_event_hits(muon_num_events);
    std::generate(
        muon_section_event_hits.begin(),
        muon_section_event_hits.end(),
        []
        {
            std::vector<Int_t> section_hits(4);
            std::fill(section_hits.begin(), section_hits.end(), 0.);
            return section_hits;
        }
    );

    const Int_t pion_num_events = (Int_t)pion_total_tree->GetEntries();
    std::cout << "Number of pion events: " << pion_num_events << std::endl;

    // Vector that holds edep in each ecal block in each event
    std::vector<std::vector<Double_t>> pion_ecal_event_edep(pion_num_events);
    std::generate(
        pion_ecal_event_edep.begin(),
        pion_ecal_event_edep.end(),
        []
        {
            std::vector<Double_t> block_event_edep(num_blocks);
            std::fill(block_event_edep.begin(), block_event_edep.end(), 0.);
            return block_event_edep;
        }
    );

    // Vector that holds edep in each hcal section in each tower in each event
    std::vector<std::vector<std::vector<Double_t>>> pion_hcal_event_edep(pion_num_events);
    std::generate(
        pion_hcal_event_edep.begin(),
        pion_hcal_event_edep.end(),
        []
        {
            std::vector<std::vector<Double_t>> tower_event_edep(num_towers);
            std::generate(
                tower_event_edep.begin(),
                tower_event_edep.end(),
                []
                {
                    std::vector<Double_t> section_event_edep(3);
                    std::fill(section_event_edep.begin(), section_event_edep.end(), 0.);
                    return section_event_edep;
                });
            return tower_event_edep;
        }
    );

    // Vector that holds edep in each detector section in each event
    std::vector<std::vector<Double_t>> pion_section_event_edep(pion_num_events);
    std::generate(
        pion_section_event_edep.begin(),
        pion_section_event_edep.end(),
        []
        {
            std::vector<Double_t> section_event_edep(4);
            std::fill(section_event_edep.begin(), section_event_edep.end(), 0.);
            return section_event_edep;
    });

    // Vector that holds number of hits in each detector section in each event
    std::vector<std::vector<Int_t>> pion_section_event_hits(pion_num_events);
    std::generate(
        pion_section_event_hits.begin(),
        pion_section_event_hits.end(),
        []
        {
            std::vector<Int_t> section_hits(4);
            std::fill(section_hits.begin(), section_hits.end(), 0.);
            return section_hits;
        }
    );


    // Different energy cuts for towers/blocks
    Double_t energy_cuts[5] = {1., 1.5, 2., 2.5, 3.};
    // Corresponding avg. hits for each cut
    Double_t muon_ecal_hits[5] = {0., 0., 0., 0., 0.};
    Double_t muon_hcal_sec1_hits[5] = {0., 0., 0., 0., 0.};
    Double_t muon_hcal_sec2_hits[5] = {0., 0., 0., 0., 0.};
    Double_t muon_hcal_sec3_hits[5] = {0., 0., 0., 0., 0.};
    Double_t muon_section_total_hits[4] = {0., 0., 0., 0.};

    // Main Energy cut used throughout the code
    Double_t main_energy_cut = 2.;
    
    
    for (Int_t ievent = 0; ievent < muon_num_events; ievent++)
    {   
        // Each event has num_X entries, so iX + ievent*num_X gets all entries
        for (Int_t iblock = 0; iblock < num_blocks; iblock++)
        {
            muon_ecal_block_tree->GetEntry(iblock + ievent * num_blocks);

            if(muon_ecal_block_edep > energy_cuts[0])
                muon_ecal_hits[0] += 1.;
            if(muon_ecal_block_edep > energy_cuts[1])
                muon_ecal_hits[1] += 1.;
            if(muon_ecal_block_edep > energy_cuts[2])
                muon_ecal_hits[2] += 1.;
            if(muon_ecal_block_edep > energy_cuts[3])
                muon_ecal_hits[3] += 1.;
            if(muon_ecal_block_edep > energy_cuts[4])
                muon_ecal_hits[4] += 1.;
            if(muon_ecal_block_edep > main_energy_cut)
            {
                muon_section_total_hits[0] += 1.;
                muon_section_event_hits[muon_ecal_block_eventID][0] += 1;
            }

            Int_t block_number = 6 * (muon_ecal_XblockID) + (muon_ecal_YblockID);
            muon_ecal_event_edep[muon_ecal_block_eventID][block_number] += muon_ecal_block_edep;
            muon_section_event_edep[muon_ecal_block_eventID][0] += muon_ecal_block_edep;
        }
        
        for (Int_t itile = 0; itile < num_towers*num_total_layers; itile++)
        {
            muon_hcal_tile_tree->GetEntry(itile + ievent * num_towers*num_total_layers);
            Int_t tower_number = 6 * muon_hcal_XtowerID + muon_hcal_YtowerID;

            // HCal Section 1
            if(muon_hcal_layerID < 9){
                muon_hcal_event_edep[muon_hcal_tile_eventID][tower_number][0] += muon_hcal_tile_edep;
                muon_section_event_edep[muon_ecal_block_eventID][1] += muon_hcal_tile_edep;
            }
            // HCal Section 2
            else if (muon_hcal_layerID >= 9 && muon_hcal_layerID < 18){
                muon_hcal_event_edep[muon_hcal_tile_eventID][tower_number][1] += muon_hcal_tile_edep;
                muon_section_event_edep[muon_ecal_block_eventID][2] += muon_hcal_tile_edep;
            }
            // HCal Section 3
            else if(muon_hcal_layerID >= 18 && muon_hcal_layerID < 51)
            {
                muon_hcal_event_edep[muon_hcal_tile_eventID][tower_number][2] += muon_hcal_tile_edep;
                muon_section_event_edep[muon_ecal_block_eventID][3] += muon_hcal_tile_edep;
            }
        }
    }
    
    // Cuts use total tower energy so need to go through HCal data after tile data
    for (Int_t ievent = 0; ievent < muon_num_events; ievent++)
    {
        for (Int_t itower = 0; itower < num_towers; itower++)
        {
            if (muon_hcal_event_edep[ievent][itower][0] > energy_cuts[0])
                muon_hcal_sec1_hits[0] += 1.;
            if (muon_hcal_event_edep[ievent][itower][0] > energy_cuts[1])
                muon_hcal_sec1_hits[1] += 1.;
            if (muon_hcal_event_edep[ievent][itower][0] > energy_cuts[2])
                muon_hcal_sec1_hits[2] += 1.;
            if (muon_hcal_event_edep[ievent][itower][0] > energy_cuts[3])
                muon_hcal_sec1_hits[3] += 1.;
            if (muon_hcal_event_edep[ievent][itower][0] > energy_cuts[4])
                muon_hcal_sec1_hits[4] += 1.;

            if (muon_hcal_event_edep[ievent][itower][1] > energy_cuts[0])
                muon_hcal_sec2_hits[0] += 1.;
            if (muon_hcal_event_edep[ievent][itower][1] > energy_cuts[1])
                muon_hcal_sec2_hits[1] += 1.;
            if (muon_hcal_event_edep[ievent][itower][1] > energy_cuts[2])
                muon_hcal_sec2_hits[2] += 1.;
            if (muon_hcal_event_edep[ievent][itower][1] > energy_cuts[3])
                muon_hcal_sec2_hits[3] += 1.;
            if (muon_hcal_event_edep[ievent][itower][1] > energy_cuts[4])
                muon_hcal_sec2_hits[4] += 1.;

            if (muon_hcal_event_edep[ievent][itower][2] > energy_cuts[0])
                muon_hcal_sec3_hits[0] += 1.;
            if (muon_hcal_event_edep[ievent][itower][2] > energy_cuts[1])
                muon_hcal_sec3_hits[1] += 1.;
            if (muon_hcal_event_edep[ievent][itower][2] > energy_cuts[2])
                muon_hcal_sec3_hits[2] += 1.;
            if (muon_hcal_event_edep[ievent][itower][2] > energy_cuts[3])
                muon_hcal_sec3_hits[3] += 1.;
            if (muon_hcal_event_edep[ievent][itower][2] > energy_cuts[4])
                muon_hcal_sec3_hits[4] += 1.;

            if (muon_hcal_event_edep[ievent][itower][0] > main_energy_cut)
            {   
                muon_section_total_hits[1] += 1.;
                muon_section_event_hits[ievent][1] += 1;
            }                
            if (muon_hcal_event_edep[ievent][itower][1] > main_energy_cut)
            {
                muon_section_total_hits[2] += 1.;
                muon_section_event_hits[ievent][2] += 1;
            }
            if (muon_hcal_event_edep[ievent][itower][2] > main_energy_cut)
            {
                muon_section_total_hits[3] += 1.;
                muon_section_event_hits[ievent][3] += 1;
            }
        }
    }

    // Muon hits in each section 
    TH1D *h_muon_sec1_hits = new TH1D("h_muon_sec1_hits", "", 20, 0, 20);
    TH1D *h_muon_sec2_hits = new TH1D("h_muon_sec2_hits", "", 20, 0, 20);
    TH1D *h_muon_sec3_hits = new TH1D("h_muon_sec3_hits", "", 20, 0, 20);
    TH1D *h_muon_sec4_hits = new TH1D("h_muon_sec4_hits", "", 20, 0, 20);

    // Muon edep in each section
    TH1D *muon_h_ecal_spectra = new TH1D("muon_h_ecal_spectra", "", 100, 0, muon_max_energy[energy]);
    TH1D *muon_h_hcal_spectra_sec1 = new TH1D("muon_h_hcal_spectra_sec1", "", 80, 0, muon_max_energy[energy]);
    TH1D *muon_h_hcal_spectra_sec2 = new TH1D("muon_h_hcal_spectra_sec2", "", 80, 0, muon_max_energy[energy]);
    TH1D *muon_h_hcal_spectra_sec3 = new TH1D("muon_h_hcal_spectra_sec3", "", 80, 0, muon_max_energy[energy]);

    // Edep ratios in hcal sections
    TH1D *muon_h_hcal_sec13_ratio = new TH1D("muon_h_hcal_sec13_ratio", "", 80, 0, 2);
    TH1D *muon_h_hcal_sec23_ratio = new TH1D("muon_h_hcal_sec23_ratio", "", 80, 0, 2);
    
    Double_t muon_efficiencies[3] = {0., 0., 0.};
    Double_t muon_20deg_efficiency = 0.;
    Double_t muon_ecal_acceptance = 0.;

    for (Int_t ievent = 0; ievent < muon_num_events; ievent++)
    {
        // MIP behavior in different sections
        if(muon_section_event_hits[ievent][0] == 1 && 
        muon_section_event_hits[ievent][1] == 1 &&
        muon_section_event_hits[ievent][2] == 1 &&
        muon_section_event_hits[ievent][3] == 1)
        {
            muon_efficiencies[0] += 1.;
        }
        if(muon_section_event_hits[ievent][0] == 1 && 
        muon_section_event_hits[ievent][1] == 1 &&
        muon_section_event_hits[ievent][2] == 1)
        {
            muon_efficiencies[1] += 1.;
        }
        if(muon_section_event_hits[ievent][0] == 1 && 
        muon_section_event_hits[ievent][1] == 1)
        {            
            muon_efficiencies[2] += 1.;
        }

        h_muon_sec1_hits->Fill(muon_section_event_hits[ievent][0]);
        h_muon_sec2_hits->Fill(muon_section_event_hits[ievent][1]);
        h_muon_sec3_hits->Fill(muon_section_event_hits[ievent][2]);
        h_muon_sec4_hits->Fill(muon_section_event_hits[ievent][3]);

        if(angle == 20)
        {
            if(muon_section_event_edep[ievent][0] < ener_to_ecal_cut[energy] &&
            muon_section_event_edep[ievent][1] < 17. &&
            muon_section_event_edep[ievent][2] < 17. &&
            muon_section_event_edep[ievent][3] > ener_to_hcal_cut[energy])
            {
                muon_20deg_efficiency += 1.;
            }
            muon_h_hcal_sec13_ratio->Fill(muon_section_event_edep[ievent][1] / muon_section_event_edep[ievent][3]);
            muon_h_hcal_sec23_ratio->Fill(muon_section_event_edep[ievent][2] / muon_section_event_edep[ievent][3]);

            if( muon_section_event_edep[ievent][0] < ener_to_ecal_cut[energy] )
            {
                muon_ecal_acceptance += 1.;
            }
            
        }
        muon_h_ecal_spectra->Fill(muon_section_event_edep[ievent][0]);
        muon_h_hcal_spectra_sec1->Fill(muon_section_event_edep[ievent][1]);
        muon_h_hcal_spectra_sec2->Fill(muon_section_event_edep[ievent][2]);
        muon_h_hcal_spectra_sec3->Fill(muon_section_event_edep[ievent][3]);

    }

    muon_20deg_efficiency /= muon_num_events;
    muon_20deg_efficiency *= 100;
    muon_ecal_acceptance /= muon_num_events;
    muon_ecal_acceptance *= 100;

    std::cout << "Muon ECal Acceptance: " << muon_ecal_acceptance << std::endl;
    Double_t pseudorapidity = -TMath::Log(TMath::Tan(angle * TMath::DegToRad() / 2.));
    TString muon_title;
    muon_title.Form("#mu- at %0.0f GeV, #eta = %0.2f", energy, pseudorapidity);

    for (Int_t i = 0; i < 5; i++)
    {
        muon_ecal_hits[i] /= muon_num_events;
        muon_hcal_sec1_hits[i] /= muon_num_events;
        muon_hcal_sec2_hits[i] /= muon_num_events;
        muon_hcal_sec3_hits[i] /= muon_num_events;
        if(i<4)
        {
            muon_section_total_hits[i] /= muon_num_events;
        }   
    }
    for (Int_t i = 0; i < 3; i++)
    {
        muon_efficiencies[i] /= muon_num_events;
        muon_efficiencies[i] *= 100;
    }

    Double_t pion_ecal_hits[5] = {0., 0., 0., 0., 0.};
    Double_t pion_hcal_sec1_hits[5] = {0., 0., 0., 0., 0.};
    Double_t pion_hcal_sec2_hits[5] = {0., 0., 0., 0., 0.};
    Double_t pion_hcal_sec3_hits[5] = {0., 0., 0., 0., 0.};
    Double_t pion_section_total_hits[4] = {0., 0., 0., 0.};

    for (Int_t ievent = 0; ievent < pion_num_events; ievent++)
    {   
        // Each event has num_X entries, so iX + ievent*num_X gets all entries
        for (Int_t iblock = 0; iblock < num_blocks; iblock++)
        {
            pion_ecal_block_tree->GetEntry(iblock + ievent * num_blocks);
            if(pion_ecal_block_edep > energy_cuts[0])
                pion_ecal_hits[0] += 1.;
            if(pion_ecal_block_edep > energy_cuts[1])
                pion_ecal_hits[1] += 1.;
            if(pion_ecal_block_edep > energy_cuts[2])
                pion_ecal_hits[2] += 1.;
            if(pion_ecal_block_edep > energy_cuts[3])
                pion_ecal_hits[3] += 1.;
            if(pion_ecal_block_edep > energy_cuts[4])
                pion_ecal_hits[4] += 1.;

            if(pion_ecal_block_edep > main_energy_cut)
            {
                pion_section_total_hits[0] += 1.;
                pion_section_event_hits[pion_ecal_block_eventID][0] += 1;
            }

            Int_t block_number = 6 * (pion_ecal_XblockID) + (pion_ecal_YblockID);
            pion_ecal_event_edep[pion_ecal_block_eventID][block_number] += pion_ecal_block_edep;
            pion_section_event_edep[pion_ecal_block_eventID][0] += pion_ecal_block_edep;
        }
        for (Int_t itile = 0; itile < num_towers*num_total_layers; itile++)
        {
            pion_hcal_tile_tree->GetEntry(itile + ievent * num_towers*num_total_layers);
            Int_t tower_number = 6 * pion_hcal_XtowerID + pion_hcal_YtowerID;

            if(pion_hcal_layerID < 9){
                pion_hcal_event_edep[pion_hcal_tile_eventID][tower_number][0] += pion_hcal_tile_edep;
                pion_section_event_edep[pion_ecal_block_eventID][1] += pion_hcal_tile_edep;
            }
            else if (pion_hcal_layerID >= 9 && pion_hcal_layerID < 18){
                pion_hcal_event_edep[pion_hcal_tile_eventID][tower_number][1] += pion_hcal_tile_edep;
                pion_section_event_edep[pion_ecal_block_eventID][2] += pion_hcal_tile_edep;
            }
            else if(pion_hcal_layerID >= 18 && pion_hcal_layerID < 51)
            {
                pion_hcal_event_edep[pion_hcal_tile_eventID][tower_number][2] += pion_hcal_tile_edep;
                pion_section_event_edep[pion_ecal_block_eventID][3] += pion_hcal_tile_edep;
            }
        }
    }

    for (Int_t ievent = 0; ievent < pion_num_events; ievent++)
    {
        for (Int_t itower = 0; itower < num_towers; itower++)
        {
            if (pion_hcal_event_edep[ievent][itower][0] > energy_cuts[0])
                pion_hcal_sec1_hits[0] += 1.;   
            if(pion_hcal_event_edep[ievent][itower][0] > energy_cuts[1])
                pion_hcal_sec1_hits[1] += 1.;
            if(pion_hcal_event_edep[ievent][itower][0] > energy_cuts[2])
                pion_hcal_sec1_hits[2] += 1.;
            if(pion_hcal_event_edep[ievent][itower][0] > energy_cuts[3])
                pion_hcal_sec1_hits[3] += 1.;
            if(pion_hcal_event_edep[ievent][itower][0] > energy_cuts[4])
                pion_hcal_sec1_hits[4] += 1.;

            if(pion_hcal_event_edep[ievent][itower][1] > energy_cuts[0])
                pion_hcal_sec2_hits[0] += 1.;
            if(pion_hcal_event_edep[ievent][itower][1] > energy_cuts[1])
                pion_hcal_sec2_hits[1] += 1.;
            if(pion_hcal_event_edep[ievent][itower][1] > energy_cuts[2])
                pion_hcal_sec2_hits[2] += 1.;
            if(pion_hcal_event_edep[ievent][itower][1] > energy_cuts[3])
                pion_hcal_sec2_hits[3] += 1.;
            if(pion_hcal_event_edep[ievent][itower][1] > energy_cuts[4])
                pion_hcal_sec2_hits[4] += 1.;

            if(pion_hcal_event_edep[ievent][itower][2] > energy_cuts[0])
                pion_hcal_sec3_hits[0] += 1.;
            if(pion_hcal_event_edep[ievent][itower][2] > energy_cuts[1])
                pion_hcal_sec3_hits[1] += 1.;
            if(pion_hcal_event_edep[ievent][itower][2] > energy_cuts[2])
                pion_hcal_sec3_hits[2] += 1.;
            if(pion_hcal_event_edep[ievent][itower][2] > energy_cuts[3])
                pion_hcal_sec3_hits[3] += 1.;
            if(pion_hcal_event_edep[ievent][itower][2] > energy_cuts[4])
                pion_hcal_sec3_hits[4] += 1.;

            if(pion_hcal_event_edep[ievent][itower][0] > main_energy_cut)
            {
                pion_section_total_hits[1] += 1.;
                pion_section_event_hits[ievent][1] += 1;
            }
            if(pion_hcal_event_edep[ievent][itower][1] > main_energy_cut)
            {
                pion_section_total_hits[2] += 1.;
                pion_section_event_hits[ievent][2] += 1;
            }

            if(pion_hcal_event_edep[ievent][itower][2] > main_energy_cut)
            {
                pion_section_total_hits[3] += 1.;
                pion_section_event_hits[ievent][3] += 1;
            }
        }
    }

    TH1D *h_pion_sec1_hits = new TH1D("h_pion_sec1_hits", "", 20, 0, 20);
    TH1D *h_pion_sec2_hits = new TH1D("h_pion_sec2_hits", "", 20, 0, 20);
    TH1D *h_pion_sec3_hits = new TH1D("h_pion_sec3_hits", "", 20, 0, 20);
    TH1D *h_pion_sec4_hits = new TH1D("h_pion_sec4_hits", "", 20, 0, 20);

    TH1D *pion_h_ecal_spectra = new TH1D("pion_h_ecal_spectra", "", 100, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec1 = new TH1D("pion_h_hcal_spectra_sec1", "", 80, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec2 = new TH1D("pion_h_hcal_spectra_sec2", "", 80, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec3 = new TH1D("pion_h_hcal_spectra_sec3", "", 80, 0, pion_max_energy[energy]);

    // Edep for pions that are identified as muons
    TH1D *pion_h_ecal_spectra_cut = new TH1D("pion_h_ecal_spectra_cut", "", 100, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec1_cut = new TH1D("pion_h_hcal_spectra_sec1_cut", "", 80, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec2_cut = new TH1D("pion_h_hcal_spectra_sec2_cut", "", 80, 0, pion_max_energy[energy]);
    TH1D *pion_h_hcal_spectra_sec3_cut = new TH1D("pion_h_hcal_spectra_sec3_cut", "", 80, 0, pion_max_energy[energy]);

    TH1D *pion_h_hcal_sec13_ratio = new TH1D("pion_h_hcal_sec13_ratio", "", 80, 0, 5);
    TH1D *pion_h_hcal_sec23_ratio = new TH1D("pion_h_hcal_sec23_ratio", "", 80, 0, 5);

    Double_t pion_efficiencies[3] = {0., 0., 0.};
    Double_t pion_20deg_efficiency = 0.;
    Double_t pion_20deg_misidentify = 0.;

    for (Int_t ievent = 0; ievent < pion_num_events; ievent++)
    {
        if(pion_section_event_hits[ievent][0] == 1 && 
        pion_section_event_hits[ievent][1] == 1 &&
        pion_section_event_hits[ievent][2] == 1 &&
        pion_section_event_hits[ievent][3] == 1)
        {
            pion_efficiencies[0] += 1.;
        }
        if(pion_section_event_hits[ievent][0] == 1 && 
        pion_section_event_hits[ievent][1] == 1 &&
        pion_section_event_hits[ievent][2] == 1)
        {
            pion_efficiencies[1] += 1.;
        }
        if(pion_section_event_hits[ievent][0] == 1 && 
        pion_section_event_hits[ievent][1] == 1)
        {            
            pion_efficiencies[2] += 1.;
        }

        h_pion_sec1_hits->Fill(pion_section_event_hits[ievent][0]);
        h_pion_sec2_hits->Fill(pion_section_event_hits[ievent][1]);
        h_pion_sec3_hits->Fill(pion_section_event_hits[ievent][2]);
        h_pion_sec4_hits->Fill(pion_section_event_hits[ievent][3]);

        if(angle == 20)
        {
            if(pion_section_event_edep[ievent][0] < ener_to_ecal_cut[energy] &&
            pion_section_event_edep[ievent][1] < 17. &&
            pion_section_event_edep[ievent][2] < 17. &&
            pion_section_event_edep[ievent][3] > ener_to_hcal_cut[energy])
            {
                pion_20deg_misidentify += 1.;
                pion_h_ecal_spectra_cut->Fill(pion_section_event_edep[ievent][0]);
                pion_h_hcal_spectra_sec1_cut->Fill(pion_section_event_edep[ievent][1]);
                pion_h_hcal_spectra_sec2_cut->Fill(pion_section_event_edep[ievent][2]);
                pion_h_hcal_spectra_sec3_cut->Fill(pion_section_event_edep[ievent][3]);
            }
            else
            {
                pion_20deg_efficiency += 1.;
            }

            pion_h_hcal_sec13_ratio->Fill(pion_section_event_edep[ievent][1] / pion_section_event_edep[ievent][3]);
            pion_h_hcal_sec23_ratio->Fill(pion_section_event_edep[ievent][2] / pion_section_event_edep[ievent][3]);
  
        }
        pion_h_ecal_spectra->Fill(pion_section_event_edep[ievent][0]);
        pion_h_hcal_spectra_sec1->Fill(pion_section_event_edep[ievent][1]);
        pion_h_hcal_spectra_sec2->Fill(pion_section_event_edep[ievent][2]);
        pion_h_hcal_spectra_sec3->Fill(pion_section_event_edep[ievent][3]);
    }
    pion_20deg_efficiency /= pion_num_events;
    pion_20deg_efficiency *= 100;
    
    pion_20deg_misidentify /= pion_num_events;
    pion_20deg_misidentify *= 100;
    
    for (Int_t i = 0; i < 5; i++)
    {
        pion_ecal_hits[i] /= pion_num_events;
        pion_hcal_sec1_hits[i] /= pion_num_events;
        pion_hcal_sec2_hits[i] /= pion_num_events;
        pion_hcal_sec3_hits[i] /= pion_num_events;
        if(i<4)
        {
            pion_section_total_hits[i] /= pion_num_events;
        }
    }

    for (Int_t i = 0; i < 3; i++)
    {
        pion_efficiencies[i] /= pion_num_events;
        pion_efficiencies[i] *= 100; 
    }
    std::cout << "Muon identification: " << muon_20deg_efficiency << std::endl;

    std::cout << "Pion identification: " << pion_20deg_efficiency << std::endl;
    std::cout << "Pion misidentification: " << pion_20deg_misidentify << std::endl;

    TGraph *pion_gr_ecal_cut = new TGraph(5, energy_cuts, pion_ecal_hits);
    TGraph *pion_gr_hcal_sec1_cut = new TGraph(5, energy_cuts, pion_hcal_sec1_hits);
    TGraph *pion_gr_hcal_sec2_cut = new TGraph(5, energy_cuts, pion_hcal_sec2_hits);
    TGraph *pion_gr_hcal_sec3_cut = new TGraph(5, energy_cuts, pion_hcal_sec3_hits);

    TGraph *muon_gr_ecal_cut = new TGraph(5, energy_cuts, muon_ecal_hits);
    TGraph *muon_gr_hcal_sec1_cut = new TGraph(5, energy_cuts, muon_hcal_sec1_hits);
    TGraph *muon_gr_hcal_sec2_cut = new TGraph(5, energy_cuts, muon_hcal_sec2_hits);
    TGraph *muon_gr_hcal_sec3_cut = new TGraph(5, energy_cuts, muon_hcal_sec3_hits);


    TCanvas *c_pion_muon_comparison = new TCanvas("c_pion_muon_comparison", "", 2000, 1000);

    TPad *pad1 = new TPad("pad1", "pad1",0,.47,0.5,1);
    TPad *pad2 = new TPad("pad2", "pad2",.5,.47,1,1);
    TPad *pad3 = new TPad("pad3", "pad3",.25,0,.75,.45);

    pad3->Draw();
    pad2->Draw();
    pad1->Draw();

    // c_pion_muon_comparison->Divide(2, 1);
    // c_pion_muon_comparison->cd(1);

    pad1->cd();


    muon_h_ecal_spectra->Draw();
    muon_h_ecal_spectra->GetXaxis()->SetTitle("Energy deposition (MeV)");
    muon_h_ecal_spectra->GetYaxis()->SetTitle("Number of events (Log)");
    muon_h_ecal_spectra->SetTitle(muon_title);
    muon_h_ecal_spectra->SetLineColor(kGreen+2);

    muon_h_hcal_spectra_sec1->Draw("same");
    muon_h_hcal_spectra_sec1->SetLineColor(kBlack);
    muon_h_hcal_spectra_sec2->Draw("same");
    muon_h_hcal_spectra_sec2->SetLineColor(kBlue);
    muon_h_hcal_spectra_sec3->Draw("same");
    muon_h_hcal_spectra_sec3->SetLineColor(kRed);
    gStyle->SetOptStat(0);

    TLegend *muon_leg_hcal = new TLegend(0.65, 0.66, 0.85, 0.86);
    muon_leg_hcal->AddEntry(muon_h_ecal_spectra, "ECal", "l");
    muon_leg_hcal->AddEntry(muon_h_hcal_spectra_sec1, "HCal Layers 1-9", "l");
    muon_leg_hcal->AddEntry(muon_h_hcal_spectra_sec2, "HCal Layers 10-18", "l");
    muon_leg_hcal->AddEntry(muon_h_hcal_spectra_sec3, "HCal Layers 19-51", "l");
    muon_leg_hcal->SetLineWidth(0);
    muon_leg_hcal->Draw();
    muon_h_ecal_spectra->SetAxisRange(1., 10000., "y");
    gPad->SetLogy();

    // TLine *l_ecal_cut = new TLine(ener_to_ecal_cut[energy], 0, ener_to_ecal_cut[energy], 3000);
    // l_ecal_cut->SetLineColor(kGreen+2);
    // l_ecal_cut->SetLineStyle(2);
    // l_ecal_cut->SetLineWidth(1);
    // l_ecal_cut->Draw("same");

    // TLine *l_hcalsec3_cut = new TLine(ener_to_hcal_cut[energy], 0, ener_to_hcal_cut[energy], 3000);
    // l_hcalsec3_cut->SetLineColor(kRed);
    // l_hcalsec3_cut->SetLineStyle(2);
    // l_hcalsec3_cut->SetLineWidth(1);
    // l_hcalsec3_cut->Draw("same");

    // TLine *l_hcalsec12_cut = new TLine(17., 0, 17., 3000);
    // l_hcalsec12_cut->SetLineColor(kBlack);
    // l_hcalsec12_cut->SetLineStyle(2);
    // l_hcalsec12_cut->SetLineWidth(1);
    // l_hcalsec12_cut->Draw("same");

    pad2->cd();
    // c_pion_muon_comparison->cd(2);
    pion_h_ecal_spectra->Draw();
    pion_h_ecal_spectra->SetLineColor(kGreen+2);
    pion_h_ecal_spectra->GetXaxis()->SetTitle("Energy deposition (MeV)");
    pion_h_ecal_spectra->GetYaxis()->SetTitle("Number of events (Log)");
    pion_h_ecal_spectra->SetAxisRange(1., 10000., "y");

    TString pion_title;
    pion_title.Form("#pi^{+} at %0.0f GeV, #eta = %0.2f", energy, pseudorapidity);
    pion_h_ecal_spectra->SetTitle(pion_title);
    pion_h_hcal_spectra_sec1->Draw("same");
    pion_h_hcal_spectra_sec1->SetLineColor(kBlack);
    pion_h_hcal_spectra_sec2->Draw("same");
    pion_h_hcal_spectra_sec2->SetLineColor(kBlue);
    pion_h_hcal_spectra_sec3->Draw("same");
    pion_h_hcal_spectra_sec3->SetLineColor(kRed);
    gPad->SetLogy();

    // pion_h_ecal_spectra_cut->Draw("same");
    // pion_h_ecal_spectra_cut->SetFillColor(kGreen+2);
    // pion_h_ecal_spectra_cut->SetLineColorAlpha(kGreen+2, .5);
    // pion_h_hcal_spectra_sec1_cut->Draw("same");
    // pion_h_hcal_spectra_sec1_cut->SetFillColorAlpha(kBlack,.5);
    // pion_h_hcal_spectra_sec1_cut->SetLineColor(kBlack);
    // pion_h_hcal_spectra_sec1_cut->SetLineWidth(3);
    // pion_h_hcal_spectra_sec2_cut->Draw("same");
    // pion_h_hcal_spectra_sec2_cut->SetFillColorAlpha(kBlue, .3);
    // pion_h_hcal_spectra_sec2_cut->SetLineColor(kBlue);
    // pion_h_hcal_spectra_sec3_cut->Draw("same");
    // pion_h_hcal_spectra_sec3_cut->SetFillColorAlpha(kRed,.5);
    // pion_h_hcal_spectra_sec3_cut->SetLineColor(kRed);

    TLegend *pion_leg = new TLegend(0.65, 0.69, 0.85, 0.89);
    pion_leg->AddEntry(pion_h_ecal_spectra, "ECal", "l");
    pion_leg->AddEntry(pion_h_hcal_spectra_sec1, "HCal Layers 1-9", "l");
    pion_leg->AddEntry(pion_h_hcal_spectra_sec2, "HCal Layers 10-18", "l");
    pion_leg->AddEntry(pion_h_hcal_spectra_sec3, "HCal Layers 19-51", "l");
    pion_leg->SetLineWidth(0);
    pion_leg->Draw();

    pad3->cd();
    // pion_gr_ecal_cut->Draw("ap");
    
    pion_gr_ecal_cut->GetXaxis()->SetTitle("Energy cut (MeV)");
    pion_gr_ecal_cut->GetYaxis()->SetTitle("Avg. Multiplicity");
    pion_gr_ecal_cut->SetMarkerStyle(kFullCircle);
    pion_gr_ecal_cut->SetMarkerColor(kGreen + 2);
    pion_gr_ecal_cut->SetTitle(pion_title);
    pion_gr_ecal_cut->GetYaxis()->SetRangeUser(0, ener_to_mult[energy]);

    pion_gr_hcal_sec1_cut->Draw("sameP");
    pion_gr_hcal_sec1_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec1_cut->SetMarkerColor(kBlack);

    pion_gr_hcal_sec2_cut->Draw("sameP");
    pion_gr_hcal_sec2_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec2_cut->SetMarkerColor(kBlue);

    pion_gr_hcal_sec3_cut->Draw("sameP");
    pion_gr_hcal_sec3_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec3_cut->SetMarkerColor(kRed);

    TLegend *pion_mult_leg = new TLegend(0.67, 0.68, 0.87, 0.88);
    pion_mult_leg->AddEntry(pion_gr_ecal_cut, "ECal", "p");
    pion_mult_leg->AddEntry(pion_gr_hcal_sec1_cut, "HCal Layers 1-9", "p");
    pion_mult_leg->AddEntry(pion_gr_hcal_sec2_cut, "HCal Layers 10-18", "p");
    pion_mult_leg->AddEntry(pion_gr_hcal_sec3_cut, "HCal Layers 19-51", "p");
    pion_mult_leg->SetLineWidth(0);
    pion_mult_leg->Draw();

    pion_gr_ecal_cut->Draw("ap");
    
    pion_gr_hcal_sec1_cut->Draw("sameP");
    pion_gr_hcal_sec1_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec1_cut->SetMarkerColor(kBlack);

    pion_gr_hcal_sec2_cut->Draw("sameP");
    pion_gr_hcal_sec2_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec2_cut->SetMarkerColor(kBlue);

    pion_gr_hcal_sec3_cut->Draw("sameP");
    pion_gr_hcal_sec3_cut->SetMarkerStyle(kFullCircle);
    pion_gr_hcal_sec3_cut->SetMarkerColor(kRed);

    pion_mult_leg->Draw();

    muon_gr_ecal_cut->Draw("sameP");
    muon_gr_ecal_cut->SetMarkerStyle(kFullStar);
    muon_gr_ecal_cut->SetMarkerSize(1.5);
    muon_gr_ecal_cut->SetMarkerColor(kGreen + 2);

    muon_gr_hcal_sec1_cut->Draw("sameP");
    muon_gr_hcal_sec1_cut->SetMarkerStyle(kFullStar);
    muon_gr_hcal_sec1_cut->SetMarkerSize(1.5);
    muon_gr_hcal_sec1_cut->SetMarkerColor(kBlack);

    muon_gr_hcal_sec2_cut->Draw("sameP");
    muon_gr_hcal_sec2_cut->SetMarkerStyle(kFullStar);
    muon_gr_hcal_sec2_cut->SetMarkerSize(1.5);
    muon_gr_hcal_sec2_cut->SetMarkerColor(kBlue);

    muon_gr_hcal_sec3_cut->Draw("sameP");
    muon_gr_hcal_sec3_cut->SetMarkerStyle(kFullStar);
    muon_gr_hcal_sec3_cut->SetMarkerSize(1.5);
    muon_gr_hcal_sec3_cut->SetMarkerColor(kRed);

    TCanvas* c_hits = new TCanvas("c_hits_pion", "", 2000, 1200);
    c_hits->Divide(2, 1);
    c_hits->cd(1);

    h_muon_sec1_hits->Draw();
    h_muon_sec1_hits->SetLineColor(kGreen + 2);
    h_muon_sec1_hits->SetFillColorAlpha(kGreen + 2, 0.5);
    h_muon_sec1_hits->SetLineWidth(3);
    h_muon_sec2_hits->Draw("same");
    h_muon_sec2_hits->SetLineColor(kBlack);
    h_muon_sec2_hits->SetFillColorAlpha(kBlack, 0.3);
    h_muon_sec2_hits->SetLineWidth(3);
    h_muon_sec3_hits->Draw("same");
    h_muon_sec3_hits->SetLineColor(kBlue);
    h_muon_sec3_hits->SetFillColorAlpha(kBlue, 0.3);
    h_muon_sec3_hits->SetLineWidth(3);
    h_muon_sec4_hits->Draw("same");
    h_muon_sec4_hits->SetLineColor(kRed);
    h_muon_sec4_hits->SetFillColorAlpha(kRed, 0.5);
    h_muon_sec4_hits->SetLineWidth(3);


    h_muon_sec1_hits->SetTitle(muon_title);
    h_muon_sec1_hits->GetXaxis()->SetTitle("Number of hits");
    h_muon_sec1_hits->SetAxisRange(0., 10000., "y");
    h_muon_sec1_hits->SetAxisRange(0., 14., "x");

    TLegend *muon_hit_leg = new TLegend(0.65, 0.69, 0.85, 0.89);
    muon_hit_leg->AddEntry(h_muon_sec1_hits, "ECal", "f");
    muon_hit_leg->AddEntry(h_muon_sec2_hits, "HCal Layers 1-9", "f");
    muon_hit_leg->AddEntry(h_muon_sec3_hits, "HCal Layers 10-18", "f");
    muon_hit_leg->AddEntry(h_muon_sec4_hits, "HCal Layers 19-51", "f");
    muon_hit_leg->SetLineWidth(0);
    muon_hit_leg->SetTextSize(0.02);
    muon_hit_leg->Draw("same");

    c_hits->cd(2);
    h_pion_sec1_hits->Draw();
    h_pion_sec1_hits->SetLineColor(kGreen + 2);
    h_pion_sec1_hits->SetFillColorAlpha(kGreen + 2, 0.5);
    h_pion_sec1_hits->SetLineWidth(3);
    h_pion_sec2_hits->Draw("same");
    h_pion_sec2_hits->SetLineColor(kBlack);
    h_pion_sec2_hits->SetFillColorAlpha(kBlack, 0.3);
    h_pion_sec2_hits->SetLineWidth(3);
    h_pion_sec3_hits->Draw("same");
    h_pion_sec3_hits->SetLineColor(kBlue);
    h_pion_sec3_hits->SetFillColorAlpha(kBlue, 0.3);
    h_pion_sec3_hits->SetLineWidth(3);
    h_pion_sec4_hits->Draw("same");
    h_pion_sec4_hits->SetLineColor(kRed);
    h_pion_sec4_hits->SetFillColorAlpha(kRed, 0.5);
    h_pion_sec4_hits->SetLineWidth(3);

    h_pion_sec1_hits->SetTitle(pion_title);
    h_pion_sec1_hits->GetXaxis()->SetTitle("Number of hits");
    h_pion_sec1_hits->SetAxisRange(0., 10000., "y");
    h_pion_sec1_hits->SetAxisRange(0., 14., "x");

    TLegend *pion_hit_leg = new TLegend(0.65, 0.69, 0.85, 0.89);
    pion_hit_leg->AddEntry(h_pion_sec1_hits, "ECal", "f");
    pion_hit_leg->AddEntry(h_pion_sec2_hits, "HCal Layers 1-9", "f");
    pion_hit_leg->AddEntry(h_pion_sec3_hits, "HCal Layers 10-18", "f");
    pion_hit_leg->AddEntry(h_pion_sec4_hits, "HCal Layers 19-51", "f");
    pion_hit_leg->SetTextSize(0.02);
    pion_hit_leg->SetLineWidth(0);
    pion_hit_leg->Draw("same");

    TCanvas* c_hcal_ratio = new TCanvas("c_hcal_ratio", "", 2000, 1200);
    c_hcal_ratio->Divide(2, 2);

    c_hcal_ratio->cd(1);
    muon_h_hcal_sec13_ratio->Draw();
    muon_h_hcal_sec13_ratio->SetTitle(muon_title);
    muon_h_hcal_sec13_ratio->GetXaxis()->SetTitle("HCal Section 1/Section 3");

    c_hcal_ratio->cd(2);
    muon_h_hcal_sec23_ratio->Draw();
    muon_h_hcal_sec23_ratio->SetTitle(muon_title);
    muon_h_hcal_sec23_ratio->GetXaxis()->SetTitle("HCal Section 2/Section 3");

    c_hcal_ratio->cd(3);
    pion_h_hcal_sec13_ratio->Draw();
    pion_h_hcal_sec13_ratio->SetTitle(pion_title);
    pion_h_hcal_sec13_ratio->GetXaxis()->SetTitle("HCal Section 1/Section 3");

    c_hcal_ratio->cd(4);
    pion_h_hcal_sec23_ratio->Draw();
    pion_h_hcal_sec23_ratio->SetTitle(pion_title);
    pion_h_hcal_sec23_ratio->GetXaxis()->SetTitle("HCal Section 2/Section 3");

    TH1D* h_data = new TH1D("h_data", " ", 16, 0.5, 16.5);
    h_data->SetBinContent(1, muon_section_total_hits[0]);
    h_data->SetBinContent(2, muon_section_total_hits[1]);
    h_data->SetBinContent(3, muon_section_total_hits[2]);
    h_data->SetBinContent(4, muon_section_total_hits[3]);
    h_data->SetBinContent(5, muon_efficiencies[0]);
    h_data->SetBinContent(6, muon_efficiencies[1]);
    h_data->SetBinContent(7, muon_efficiencies[2]);

    h_data->SetBinContent(8, pion_section_total_hits[0]);
    h_data->SetBinContent(9, pion_section_total_hits[1]);
    h_data->SetBinContent(10, pion_section_total_hits[2]);
    h_data->SetBinContent(11, pion_section_total_hits[3]);
    h_data->SetBinContent(12, pion_efficiencies[0]);
    h_data->SetBinContent(13, pion_efficiencies[1]);
    h_data->SetBinContent(14, pion_efficiencies[2]);

    h_data->SetBinContent(15, muon_20deg_efficiency);
    h_data->SetBinContent(16, pion_20deg_misidentify);

    TString outname;
    outname.Form("./mupi_%ddeg/mupi_%0.0fGeV_%ddeg.root", angle, energy, angle);
    TFile* outfile = new TFile(outname, "RECREATE");
    h_data->Write();
}

void energy_plots()
{
    TFile *mupi1GeV_5deg_file = new TFile("./mupi_5deg/mupi_1GeV_5deg.root");
    TFile* mupi2GeV_5deg_file = new TFile("./mupi_5deg/mupi_2GeV_5deg.root");
    TFile* mupi5GeV_5deg_file = new TFile("./mupi_5deg/mupi_5GeV_5deg.root");
    TFile* mupi10GeV_5deg_file = new TFile("./mupi_5deg/mupi_10GeV_5deg.root");
    TFile* mupi20GeV_5deg_file = new TFile("./mupi_5deg/mupi_20GeV_5deg.root");
    TFile* mupi30GeV_5deg_file = new TFile("./mupi_5deg/mupi_30GeV_5deg.root");
    TFile* mupi40GeV_5deg_file = new TFile("./mupi_5deg/mupi_40GeV_5deg.root");
    TFile* mupi50GeV_5deg_file = new TFile("./mupi_5deg/mupi_50GeV_5deg.root");
    TFile* mupi60GeV_5deg_file = new TFile("./mupi_5deg/mupi_60GeV_5deg.root");
    TFile* mupi70GeV_5deg_file = new TFile("./mupi_5deg/mupi_70GeV_5deg.root");
    TFile* mupi80GeV_5deg_file = new TFile("./mupi_5deg/mupi_80GeV_5deg.root");
    TFile* mupi90GeV_5deg_file = new TFile("./mupi_5deg/mupi_90GeV_5deg.root");
    TFile* mupi100GeV_5deg_file = new TFile("./mupi_5deg/mupi_100GeV_5deg.root");

    TFile *mupi1GeV_20deg_file = new TFile("./mupi_20deg/mupi_1GeV_20deg.root");
    TFile* mupi2GeV_20deg_file = new TFile("./mupi_20deg/mupi_2GeV_20deg.root");
    TFile* mupi5GeV_20deg_file = new TFile("./mupi_20deg/mupi_5GeV_20deg.root");
    TFile* mupi10GeV_20deg_file = new TFile("./mupi_20deg/mupi_10GeV_20deg.root");
    TFile* mupi20GeV_20deg_file = new TFile("./mupi_20deg/mupi_20GeV_20deg.root");

    TH1D* mupi1GeV_5deg_data = new TH1D;
    TH1D* mupi2GeV_5deg_data = new TH1D;
    TH1D* mupi5GeV_5deg_data = new TH1D;
    TH1D* mupi10GeV_5deg_data = new TH1D;
    TH1D* mupi20GeV_5deg_data = new TH1D;
    TH1D* mupi30GeV_5deg_data = new TH1D;
    TH1D* mupi40GeV_5deg_data = new TH1D;
    TH1D* mupi50GeV_5deg_data = new TH1D;
    TH1D* mupi60GeV_5deg_data = new TH1D;
    TH1D* mupi70GeV_5deg_data = new TH1D;
    TH1D* mupi80GeV_5deg_data = new TH1D;
    TH1D* mupi90GeV_5deg_data = new TH1D;
    TH1D* mupi100GeV_5deg_data = new TH1D;

    TH1D* mupi1GeV_20deg_data = new TH1D;
    TH1D* mupi2GeV_20deg_data = new TH1D;
    TH1D* mupi5GeV_20deg_data = new TH1D;
    TH1D* mupi10GeV_20deg_data = new TH1D;
    TH1D* mupi20GeV_20deg_data = new TH1D;

    mupi1GeV_5deg_file->GetObject("h_data", mupi1GeV_5deg_data);
    mupi2GeV_5deg_file->GetObject("h_data", mupi2GeV_5deg_data);
    mupi5GeV_5deg_file->GetObject("h_data", mupi5GeV_5deg_data);
    mupi10GeV_5deg_file->GetObject("h_data", mupi10GeV_5deg_data);
    mupi20GeV_5deg_file->GetObject("h_data", mupi20GeV_5deg_data);
    mupi30GeV_5deg_file->GetObject("h_data", mupi30GeV_5deg_data);
    mupi40GeV_5deg_file->GetObject("h_data", mupi40GeV_5deg_data);
    mupi50GeV_5deg_file->GetObject("h_data", mupi50GeV_5deg_data);
    mupi60GeV_5deg_file->GetObject("h_data", mupi60GeV_5deg_data);
    mupi70GeV_5deg_file->GetObject("h_data", mupi70GeV_5deg_data);
    mupi80GeV_5deg_file->GetObject("h_data", mupi80GeV_5deg_data);
    mupi90GeV_5deg_file->GetObject("h_data", mupi90GeV_5deg_data);
    mupi100GeV_5deg_file->GetObject("h_data", mupi100GeV_5deg_data);

    mupi1GeV_20deg_file->GetObject("h_data", mupi1GeV_20deg_data);
    mupi2GeV_20deg_file->GetObject("h_data", mupi2GeV_20deg_data);
    mupi5GeV_20deg_file->GetObject("h_data", mupi5GeV_20deg_data);
    mupi10GeV_20deg_file->GetObject("h_data", mupi10GeV_20deg_data);
    mupi20GeV_20deg_file->GetObject("h_data", mupi20GeV_20deg_data);


    Double_t energies[13] = {1., 2., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
    Double_t energies_20deg[5] = {1., 2., 5., 10., 20.};

    std::vector<Double_t> muon_nhits_sec1;
    std::vector<Double_t> muon_nhits_sec2;
    std::vector<Double_t> muon_nhits_sec3;
    std::vector<Double_t> muon_nhits_sec4;
    std::vector<Double_t> pion_nhits_sec1;
    std::vector<Double_t> pion_nhits_sec2;
    std::vector<Double_t> pion_nhits_sec3;
    std::vector<Double_t> pion_nhits_sec4;
    std::vector<Double_t> muon_efficiency_whole;
    std::vector<Double_t> muon_efficiency_wolast;
    std::vector<Double_t> muon_efficiency_wolasttwo;
    std::vector<Double_t> pion_misidentify_whole;
    std::vector<Double_t> pion_misidentify_wolast;
    std::vector<Double_t> pion_misidentify_wolasttwo;

    std::vector<Double_t> muon_efficiency_20deg;
    std::vector<Double_t> pion_misidentify_20deg;

    std::map<Int_t, TH1D*> hist_5deg_map = {{1, mupi1GeV_5deg_data}, {2, mupi2GeV_5deg_data}, {3, mupi5GeV_5deg_data}, {4, mupi10GeV_5deg_data}, {5, mupi20GeV_5deg_data}, {6, mupi30GeV_5deg_data}, {7, mupi40GeV_5deg_data}, 
    {8, mupi50GeV_5deg_data}, {9, mupi60GeV_5deg_data}, {10, mupi70GeV_5deg_data}, {11, mupi80GeV_5deg_data}, {12, mupi90GeV_5deg_data}, {13, mupi100GeV_5deg_data}};

    std::map<Int_t, TH1D *> hist_20deg_map = {{1, mupi1GeV_20deg_data}, {2, mupi2GeV_20deg_data}, {3, mupi5GeV_20deg_data}, {4, mupi10GeV_20deg_data}, {5, mupi20GeV_20deg_data}};

    for(Int_t i = 1; i < 14; i++)
    {
        muon_nhits_sec1.push_back(hist_5deg_map[i]->GetBinContent(1));
        muon_nhits_sec2.push_back(hist_5deg_map[i]->GetBinContent(2));
        muon_nhits_sec3.push_back(hist_5deg_map[i]->GetBinContent(3));
        muon_nhits_sec4.push_back(hist_5deg_map[i]->GetBinContent(4));
        muon_efficiency_whole.push_back(hist_5deg_map[i]->GetBinContent(5));
        muon_efficiency_wolast.push_back(hist_5deg_map[i]->GetBinContent(6));
        muon_efficiency_wolasttwo.push_back(hist_5deg_map[i]->GetBinContent(7));
        pion_nhits_sec1.push_back(hist_5deg_map[i]->GetBinContent(8));
        pion_nhits_sec2.push_back(hist_5deg_map[i]->GetBinContent(9));
        pion_nhits_sec3.push_back(hist_5deg_map[i]->GetBinContent(10));
        pion_nhits_sec4.push_back(hist_5deg_map[i]->GetBinContent(11));
        pion_misidentify_whole.push_back(hist_5deg_map[i]->GetBinContent(12));
        pion_misidentify_wolast.push_back(hist_5deg_map[i]->GetBinContent(13));
        pion_misidentify_wolasttwo.push_back(hist_5deg_map[i]->GetBinContent(14));

        if(i < 6)
        {
            muon_efficiency_20deg.push_back(hist_20deg_map[i]->GetBinContent(15));
            pion_misidentify_20deg.push_back(hist_20deg_map[i]->GetBinContent(16));
        }
    }


    TCanvas *c_hits = new TCanvas("c_hits", "", 1200, 800);
    TGraph* gr_muon_nhits_sec1 = new TGraph(13, energies, muon_nhits_sec1.data());
    TGraph* gr_muon_nhits_sec2 = new TGraph(13, energies, muon_nhits_sec2.data());
    TGraph* gr_muon_nhits_sec3 = new TGraph(13, energies, muon_nhits_sec3.data());
    TGraph* gr_muon_nhits_sec4 = new TGraph(13, energies, muon_nhits_sec4.data());

    TGraph* gr_pion_nhits_sec1 = new TGraph(13, energies, pion_nhits_sec1.data());
    TGraph* gr_pion_nhits_sec2 = new TGraph(13, energies, pion_nhits_sec2.data());
    TGraph* gr_pion_nhits_sec3 = new TGraph(13, energies, pion_nhits_sec3.data());
    TGraph* gr_pion_nhits_sec4 = new TGraph(13, energies, pion_nhits_sec4.data());

    gr_muon_nhits_sec1->Draw("ap");
    gr_muon_nhits_sec1->GetYaxis()->SetRangeUser(0, 15);
    gr_muon_nhits_sec1->GetXaxis()->SetRangeUser(-1, 101);

    gr_muon_nhits_sec1->SetMarkerStyle(kFullStar);
    gr_muon_nhits_sec1->SetMarkerSize(1.5);
    gr_muon_nhits_sec1->SetMarkerColor(kGreen + 2);
    gr_muon_nhits_sec2->Draw("sameP");
    gr_muon_nhits_sec2->SetMarkerStyle(kFullStar);
    gr_muon_nhits_sec2->SetMarkerColor(kBlack);
    gr_muon_nhits_sec2->SetMarkerSize(1.5);
    gr_muon_nhits_sec3->Draw("sameP");
    gr_muon_nhits_sec3->SetMarkerStyle(kFullStar);
    gr_muon_nhits_sec3->SetMarkerColor(kBlue);
    gr_muon_nhits_sec3->SetMarkerSize(1.5);
    gr_muon_nhits_sec4->Draw("sameP");
    gr_muon_nhits_sec4->SetMarkerStyle(kFullStar);
    gr_muon_nhits_sec4->SetMarkerColor(kRed);
    gr_muon_nhits_sec4->SetMarkerSize(1.5);

    gr_pion_nhits_sec1->Draw("sameP");
    gr_pion_nhits_sec1->SetMarkerStyle(kFullCircle);
    gr_pion_nhits_sec1->SetMarkerColor(kGreen + 2);
    gr_pion_nhits_sec2->Draw("sameP");
    gr_pion_nhits_sec2->SetMarkerStyle(kFullCircle);
    gr_pion_nhits_sec2->SetMarkerColor(kBlack);
    gr_pion_nhits_sec3->Draw("sameP");
    gr_pion_nhits_sec3->SetMarkerStyle(kFullCircle);
    gr_pion_nhits_sec3->SetMarkerColor(kBlue);
    gr_pion_nhits_sec4->Draw("sameP");
    gr_pion_nhits_sec4->SetMarkerStyle(kFullCircle);
    gr_pion_nhits_sec4->SetMarkerColor(kRed);

    gr_muon_nhits_sec1->GetYaxis()->SetTitle("Avg. Multiplicity");
    gr_muon_nhits_sec1->GetXaxis()->SetTitle("Energy (GeV)");
    gr_muon_nhits_sec1->SetTitle("2 MeV threshold");

    TLegend *pion_mult_leg = new TLegend(0.13, 0.68, 0.33, 0.88);
    pion_mult_leg->AddEntry(gr_pion_nhits_sec1, "ECal Block", "p");
    pion_mult_leg->AddEntry(gr_pion_nhits_sec2, "HCal Layers 1-9", "p");
    pion_mult_leg->AddEntry(gr_pion_nhits_sec3, "HCal Layers 10-18", "p");
    pion_mult_leg->AddEntry(gr_pion_nhits_sec4, "HCal Layers 19-51", "p");
    pion_mult_leg->SetLineWidth(0);
    pion_mult_leg->Draw();


    TGraph* gr_muon_efficiency_whole = new TGraph(13, energies, muon_efficiency_whole.data());
    TGraph* gr_muon_efficiency_wolast = new TGraph(13, energies, muon_efficiency_wolast.data());
    TGraph* gr_muon_efficiency_wolast2 = new TGraph(13, energies, muon_efficiency_wolasttwo.data());

    TGraph* gr_pion_misidentify_whole = new TGraph(13, energies, pion_misidentify_whole.data());
    TGraph* gr_pion_misidentify_wolast = new TGraph(13, energies, pion_misidentify_wolast.data());
    TGraph* gr_pion_misidentify_wolast2 = new TGraph(13, energies, pion_misidentify_wolasttwo.data());

    TCanvas *c_efficiency = new TCanvas("c_efficiency", "", 1200, 800);


    gr_muon_efficiency_wolast->Draw("ap");
    gr_muon_efficiency_wolast->SetMarkerStyle(kFullStar);
    gr_muon_efficiency_wolast->SetMarkerSize(1.5);
    gr_muon_efficiency_wolast->SetMarkerColor(kRed);
    gr_muon_efficiency_wolast2->Draw("sameP");
    gr_muon_efficiency_wolast2->SetMarkerStyle(kFullStar);
    gr_muon_efficiency_wolast2->SetMarkerColor(kBlue);
    gr_muon_efficiency_wolast2->SetMarkerSize(1.5);
    gr_muon_efficiency_wolast->GetYaxis()->SetRangeUser(-1, 101);
    gr_muon_efficiency_wolast->GetXaxis()->SetRangeUser(-1, 101);

    gr_pion_misidentify_wolast->Draw("sameP");
    gr_pion_misidentify_wolast->SetMarkerStyle(kFullCircle);
    gr_pion_misidentify_wolast->SetMarkerColor(kRed);
    gr_pion_misidentify_wolast2->Draw("sameP");
    gr_pion_misidentify_wolast2->SetMarkerStyle(kFullCircle);
    gr_pion_misidentify_wolast2->SetMarkerColor(kBlue);

    gr_muon_efficiency_wolast->GetYaxis()->SetTitle("\% of events with one hit in each section");
    gr_muon_efficiency_wolast->GetXaxis()->SetTitle("Energy (GeV)");
    gr_muon_efficiency_wolast->SetTitle("");

    TLegend *pion_eff_leg = new TLegend(0.62, 0.48, 0.87, 0.6);
    pion_eff_leg->AddEntry(gr_pion_misidentify_wolast, "Sections 1-3", "p");
    pion_eff_leg->AddEntry(gr_pion_misidentify_wolast2, "Sections 1-2", "p");
    pion_eff_leg->SetLineWidth(0);
    pion_eff_leg->Draw();


    TGraph* gr_muon_20deg_efficiency = new TGraph(5, energies_20deg, muon_efficiency_20deg.data());
    TGraph* gr_pion_misidentify_20deg = new TGraph(5, energies_20deg, pion_misidentify_20deg.data());

    TCanvas *c_20deg_identification = new TCanvas("c_20deg_identification", "", 1200, 800);


    gr_muon_20deg_efficiency->Draw("ap");
    gr_muon_20deg_efficiency->SetMarkerStyle(kFullStar);
    gr_muon_20deg_efficiency->SetMarkerSize(1.5);
    gr_muon_20deg_efficiency->SetMarkerColor(kBlue);
    gr_muon_20deg_efficiency->GetYaxis()->SetRangeUser(-1, 101);
    gr_muon_20deg_efficiency->GetXaxis()->SetRangeUser(0, 20.2);

    gr_muon_efficiency_wolast->Draw("sameP");

    gr_pion_misidentify_20deg->Draw("sameP");
    // gr_pion_misidentify_20deg->SetMarkerSize(1.65);
    gr_pion_misidentify_20deg->SetMarkerStyle(kFullCircle);
    gr_pion_misidentify_20deg->SetMarkerColor(kBlue);

    gr_pion_misidentify_wolast->Draw("sameP");
    // gr_pion_misidentify_wolast->SetMarkerSize(1.65);

    gr_muon_20deg_efficiency->GetYaxis()->SetTitle("\% of events identified as #mu-");
    gr_muon_20deg_efficiency->GetXaxis()->SetTitle("Energy (GeV)");
    gr_muon_20deg_efficiency->SetTitle("");

    TLegend *eff_20deg_leg = new TLegend(0.45, 0.33, 0.87, 0.67);
    eff_20deg_leg->AddEntry(gr_muon_20deg_efficiency, "#mu-, #eta = 1.74, ECal + 3 HCal Section Cuts", "p");
    eff_20deg_leg->AddEntry(gr_muon_efficiency_wolast, "#mu-, #eta = 3.13, MIP in ECal + 2 HCal Sections", "p");
    eff_20deg_leg->AddEntry(gr_pion_misidentify_20deg, "#pi^{+}, #eta = 1.74, ECal + 3 HCal Section Cuts", "p");
    eff_20deg_leg->AddEntry(gr_pion_misidentify_wolast, "#pi^{+}, #eta = 3.13, MIP in ECal + 2 HCal Sections", "p");
    eff_20deg_leg->SetLineWidth(0);
    eff_20deg_leg->Draw();

    TLine *l_90 = new TLine(0., 90., 20.4, 90.);
    l_90->SetLineStyle(2);
    l_90->Draw("same");

}