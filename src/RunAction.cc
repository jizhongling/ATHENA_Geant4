
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories 
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple

  analysisManager->CreateNtuple("EdepTotal", "Edep");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Total");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Total");
  analysisManager->CreateNtupleIColumn("ECal_NumHits_Total");
  analysisManager->CreateNtupleIColumn("HCal_NumHits_Total");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("ECalBlocks", "ECalBlocks");
  analysisManager->CreateNtupleDColumn("ECal_Edep_Block");
  analysisManager->CreateNtupleIColumn("ECal_BlockXid");
  analysisManager->CreateNtupleIColumn("ECal_BlockYid");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("HCalTowers", "HCalTowers");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Tower");
  analysisManager->CreateNtupleIColumn("HCal_TowerXid");
  analysisManager->CreateNtupleIColumn("HCal_TowerYid");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("HCalLayers", "HCalLayers");
  analysisManager->CreateNtupleDColumn("HCal_Edep_Tile");
  analysisManager->CreateNtupleIColumn("HCal_Layerid");
  analysisManager->CreateNtupleIColumn("HCal_NumHits_Tile");
  analysisManager->CreateNtupleIColumn("HCal_TowerXid");
  analysisManager->CreateNtupleIColumn("HCal_TowerYid");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile(analysisManager->GetFileName()); // File name set via macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  // if(isMaster) 
  // {
  //   G4cout << G4endl << " ----> print histograms statistic ";
    
  //     G4cout << "for the entire run " << G4endl << G4endl; 
    
  //   G4cout << " ECal : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    
  //   G4cout << " HCal : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
  // }
  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
