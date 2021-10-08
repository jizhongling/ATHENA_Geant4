
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"
#include "Analysis.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
 : G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<CalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                              G4double ECalEdep, G4double gapEdep) const
{
  // print event statistics
  G4cout
     << "   ECal: total energy: " 
     << std::setw(7) << G4BestUnit(ECalEdep, "Energy")
     << G4endl
     << "        HCal: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  

  const G4int NumHCalTowers = DetectorConstruction::GetHCalTowers();
  const G4int NumECalBlocks = DetectorConstruction::GetECalBlocks();

  char nameHolder[200];
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  G4double HCal_Edep = 0.;
  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalHitsCollection%d%d", i, j);
      G4int HCalHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto HCalHC = GetHitsCollection(HCalHCID, event);
      auto HCalHit = (*HCalHC)[HCalHC->entries()-1];
      HCal_Edep += HCalHit->GetEdep();
      // Ntuple with id 2 holds HCal information
      analysisManager->FillNtupleDColumn(2, 0, HCalHit->GetEdep());
      analysisManager->FillNtupleIColumn(2, 1, i);
      analysisManager->FillNtupleIColumn(2, 2, j);
      analysisManager->AddNtupleRow(2);
    }
  }

  G4double ECal_Edep = 0.;
  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      sprintf(nameHolder, "ECalHitsCollection%d%d", i, j);
      G4int ECalHCID = G4SDManager::GetSDMpointer()->GetCollectionID(nameHolder);
      auto ECalHC = GetHitsCollection(ECalHCID, event);
      auto ECalHit = (*ECalHC)[ECalHC->entries()-1];
      ECal_Edep += ECalHit->GetEdep();

      // Ntuple with id 1 holds ECal information
      analysisManager->FillNtupleDColumn(1, 0, ECalHit->GetEdep());
      analysisManager->FillNtupleIColumn(1, 1, i);
      analysisManager->FillNtupleIColumn(1, 2, j);
      analysisManager->AddNtupleRow(1);
    }
  }

  auto eventID = event->GetEventID();
  if(eventID % 1000 == 0) G4cout << "---> End of event: " << eventID << G4endl;     

  
 
  // fill histograms
  analysisManager->FillH1(0, ECal_Edep);
  analysisManager->FillH1(1, HCal_Edep);
  // fill ntuple
  
  // Ntuple with id 0 holds total information
  analysisManager->FillNtupleDColumn(0, 0, ECal_Edep);
  analysisManager->FillNtupleDColumn(0, 1, HCal_Edep);
  analysisManager->AddNtupleRow(0); 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
