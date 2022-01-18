
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "GlobalValues.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace GlobalValues;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 
 //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(false)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Fe"); // Using iron for steel
  nistManager->FindOrBuildMaterial("G4_POLYSTYRENE"); // Fiber core material
  nistManager->FindOrBuildMaterial("G4_PLEXIGLASS"); // PMMA for fiber cladding
  nistManager->FindOrBuildMaterial("G4_Galactic"); // Vacuum

  // Print materials
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4cout<<"Constructing Geometry..."<<G4endl;

  // HCal Tower Geometry parameters
  G4double AbsorberPlateThickness = 20.*mm; 
  G4double ActivePlateThickness =  3.*mm;
  G4double HCal_X = 100*mm; // x-dimension of each HCal tower
  G4double HCal_Y = 98.897*mm; // y-dimension of each HCal tower
  G4double HCal_WLS_X = 4.0*mm; // Will have a vertical wavelength-shifting plate in between towers
  G4double HCal_Steel_Y = 1.897*mm; // Will have a horizontal steel (approximated as iron) plate in between towers
  G4double HCal_LayerThickness = AbsorberPlateThickness + ActivePlateThickness; // Each layer has an absorber and active component
  G4double HCal_Thickness = NumHCalLayers * HCal_LayerThickness; // z-dimension of each HCal tower

  // ECal Block Geometry Paramters
  G4double ECal_X = 49.85*mm; // x-dimension of each ECal block
  G4double ECal_Y = 49.30*mm; // y-dimension of each ECal block
  G4double ECal_Thickness = 170*mm; // z-dimension of each ECal block
  G4double ECal_Glue_XY = 0.1*mm; // Glue connecting ECal blocks -- see design specs
  G4double Clearance_Gap = 0.1*mm; // Gap between each set of 4 blocks -- see design specs
  G4double ECal_Fiber_r = 0.235*mm; // Radius of each fiber in ECal
  G4int ECal_Fiber_Rows = 60; // Number of fiber rows in each ECal block
  G4int ECal_Fiber_Cols = 52; // Number of fiber columns in each ECal block

  auto worldSizeXY = 10 * HCal_X;
  auto worldSizeZ  = 2. * (HCal_Thickness + ECal_Thickness); // Arbitrary sizes larger than the detector
  //auto worldSizeXY = 10*m;
  //auto worldSizeZ  = 10*m; // Arbitrary sizes larger than the detector

  // Get materials
  auto DefaultMaterial = G4Material::GetMaterial("G4_Galactic");
  auto AbsorberPlateMaterial = G4Material::GetMaterial("G4_Fe");
  auto ActiveMaterial = G4Material::GetMaterial("G4_POLYSTYRENE");
  auto CladdingMaterial = G4Material::GetMaterial("G4_PLEXIGLASS");

  G4Element* elW = new G4Element("Tungsten", "W", 74., 183.85*g/mole);
  G4Material* ECalAbsorberMaterial = new G4Material("ECalAbsorberMaterial", 12.72*g/cm3, 2);
  ECalAbsorberMaterial->AddElement(elW, 97.0*perCent); // Use mass fraction
  ECalAbsorberMaterial->AddMaterial(ActiveMaterial, 3.0*perCent); // Use mass fraction

  //G4Material* ECal_abs_mat = new G4Material("ECal_fiber_mat", 10.15*g/cm3, 2);
  //ECal_abs_mat->AddElement(elW, 94.8*perCent); // Use mass fraction
  //ECal_abs_mat->AddMaterial(ActiveMaterial, 5.2*perCent); // Use mass fraction

  G4MaterialPropertiesTable* MaterialTable = new G4MaterialPropertiesTable();
  ActiveMaterial->SetMaterialPropertiesTable(MaterialTable);
  ActiveMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
  if ( !DefaultMaterial || !AbsorberPlateMaterial || !ActiveMaterial || !CladdingMaterial || !ECalAbsorberMaterial ) 
  {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  } 
     
  // World
  auto WorldS 
    = new G4Box("WorldSolid",           // its name
                 worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.); // its size
                         
  auto WorldLV
    = new G4LogicalVolume(
                 WorldS,           // its solid
                 DefaultMaterial,  // its material
                 "WorldLogical");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 WorldLV,          // its logical volume                         
                 "WorldPhysical",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  char nameHolder[200];

  // HCal
  G4LogicalVolume* HCalLV[NumHCalTowers][NumHCalTowers];
  G4VSolid* HCalS = new G4Box("HCalSolid", HCal_X/2., HCal_Y/2., HCal_Thickness/2.);

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalLogical%d%d", i, j);
      HCalLV[i][j] = new G4LogicalVolume(HCalS, DefaultMaterial, nameHolder);
      sprintf(nameHolder, "HCalPhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector((-2.5 + i)*HCal_X, (2.5 - j)*HCal_Y, ECal_Thickness/2. + HCal_Thickness/2.), HCalLV[i][j], nameHolder, WorldLV, false, i+j, fCheckOverlaps);
    }
  }
                                                     
  // HCal Layers

  // LayerHolder is used to easily replicate the layers along Z. Dimensions must account for WLS plates and steel plates
  G4LogicalVolume* HCalLayerHolderLV[NumHCalTowers][NumHCalTowers]; 
  G4VSolid* HCalLayerHolderS = new G4Box("HCalLayerHolderSolid", (HCal_X - HCal_WLS_X)/2., (HCal_Y - HCal_Steel_Y)/2., HCal_Thickness/2.);

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalLayerHolderLogical%d%d", i, j);
      HCalLayerHolderLV[i][j] = new G4LogicalVolume(HCalLayerHolderS, DefaultMaterial, nameHolder);
      sprintf(nameHolder, "HCalLayerHolderPhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector(HCal_WLS_X/2., -HCal_Steel_Y/2., 0), HCalLayerHolderLV[i][j], nameHolder, HCalLV[i][j], false, 0, fCheckOverlaps);
    }
  }

  G4LogicalVolume* HCalLayerLV[NumHCalTowers][NumHCalTowers];
  G4VSolid* HCalLayerS = new G4Box("HCalLayerSolid", (HCal_X - HCal_WLS_X)/2., (HCal_Y - HCal_Steel_Y)/2., HCal_LayerThickness/2.);

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalLayerLogical%d%d", i, j);
      HCalLayerLV[i][j] = new G4LogicalVolume(HCalLayerS, DefaultMaterial, nameHolder);
      new G4PVReplica("HCalLayerPhysical", HCalLayerLV[i][j], HCalLayerHolderLV[i][j], kZAxis, NumHCalLayers, HCal_LayerThickness);
    }
  }
  
  // Absorber plates in HCal towers
  
  G4LogicalVolume* HCalAbsorberLV[NumHCalTowers][NumHCalTowers];
  G4VSolid* HCalAbsorberS = new G4Box("HCalAbsorberSolid", (HCal_X - HCal_WLS_X)/2., (HCal_Y - HCal_Steel_Y)/2., AbsorberPlateThickness/2.);

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalAbsorberLogical%d%d", i, j);
      HCalAbsorberLV[i][j] = new G4LogicalVolume(HCalAbsorberS, AbsorberPlateMaterial, nameHolder);
      sprintf(nameHolder, "HCalAbsorberPhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector(0., 0., -ActivePlateThickness/2.), HCalAbsorberLV[i][j], nameHolder, HCalLayerLV[i][j], false, 0, fCheckOverlaps);
    }
  }

  // Scintillating plates in HCal towers

  // Behind the absorber plates in each layer
  G4LogicalVolume* HCalActiveLV[NumHCalTowers][NumHCalTowers];
  G4VSolid* HCalActiveS = new G4Box("HCalActiveSolid", (HCal_X - HCal_WLS_X)/2., (HCal_Y - HCal_Steel_Y)/2., ActivePlateThickness/2.);

  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalActiveLogical%d%d", i, j);
      HCalActiveLV[i][j] = new G4LogicalVolume(HCalActiveS, ActiveMaterial, nameHolder);
      sprintf(nameHolder, "HCalActivePhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector(0., 0., AbsorberPlateThickness/2.), HCalActiveLV[i][j], nameHolder, HCalLayerLV[i][j], false, 0, fCheckOverlaps);
    }
  }

  // Wavelength shifting plates in HCal towers 

  // Only in between towers. Implemented as being part of towers rather than separarte. In right side of towers. 
  // Far right tower section does not have WLS plates.
  G4LogicalVolume* HCalWLS_LV[NumHCalTowers-1][NumHCalTowers];
  G4VSolid* HCalWLS_S = new G4Box("HCalWLSSolid", HCal_WLS_X/2., (HCal_Y - HCal_Steel_Y)/2., HCal_Thickness/2.);
  for(G4int i = 1; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalWLSLogical%d%d", i-1, j);
      HCalWLS_LV[i-1][j] = new G4LogicalVolume(HCalWLS_S, ActiveMaterial, nameHolder);
      sprintf(nameHolder, "HCalWLSPhysical%d%d", i-1, j);
      new G4PVPlacement(0, G4ThreeVector(-(HCal_X-HCal_WLS_X)/2., -HCal_Steel_Y/2., 0), HCalWLS_LV[i-1][j], nameHolder, HCalLV[i][j], false, 0, fCheckOverlaps);
    }
  }

  // Steel plates in HCal towers

  // Only in between towers. Implemented as being part of towers rather than separarte. In top of towers.
  // Top tower section does not have steel plates.
  G4LogicalVolume* HCalSteelLV[NumHCalTowers][NumHCalTowers-1];
  G4VSolid* HCalSteelS = new G4Box("HCalSteelSolid", HCal_X/2., HCal_Steel_Y/2., HCal_Thickness/2.); // Plates that stretch across HCal in x and z directions
  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 1; j < NumHCalTowers; j++)
    {
      sprintf(nameHolder, "HCalSteelLogical%d%d", i, j-1);
      HCalSteelLV[i][j-1] = new G4LogicalVolume(HCalSteelS, AbsorberPlateMaterial, nameHolder);
      sprintf(nameHolder, "HCalSteelPhysical%d%d", i, j-1);
      new G4PVPlacement(0, G4ThreeVector(0, (HCal_Y-HCal_Steel_Y)/2. , 0), HCalSteelLV[i][j-1], nameHolder, HCalLV[i][j], false, 0, fCheckOverlaps);
    }
  }


  // ECal mixture tower structure
  //G4VSolid* ECalS = new G4Tubs("ECalSolid", 0.0*mm, 1542*mm, ECal_Thickness/2, 0, CLHEP::twopi);
  //G4LogicalVolume* ECalLV = new G4LogicalVolume(ECalS, ECal_abs_mat, "ECalLogical");
  //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), ECalLV, "ECalPhysical", WorldLV, false, 0, fCheckOverlaps);

  // // ECal Blocks
  // // First ECal block has origin at ((-2.*HCal_X + ECal_X/2. + Clearance_Gap), (2.*HCal_Y - ECal_Y/2. - Clearance_Gap)), which is top right HCal block shifted by clearance gap
  G4LogicalVolume* ECalLV[NumECalBlocks][NumECalBlocks];
  G4VSolid* ECalS = new G4Box("ECalSolid", ECal_X/2., ECal_Y/2., ECal_Thickness/2.);

  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      sprintf(nameHolder, "ECalLogical%d%d", i, j);
      ECalLV[i][j] = new G4LogicalVolume(ECalS, ECalAbsorberMaterial, nameHolder);

      G4double x0 = -2.*HCal_X + ECal_X/2. + Clearance_Gap; // Top right HCal block
      if(i % 2 != 0) x0 += (ECal_X + ECal_Glue_XY);        // Block to the left of the first block
      G4double y0 = 2.*HCal_Y - ECal_Y/2. - Clearance_Gap; // Top right HCal block
      if(j % 2 != 0) y0 -= (ECal_Y + ECal_Glue_XY);       // Block below the first block
      G4int i_factor = i/2;
      G4int j_factor = j/2;
      sprintf(nameHolder, "ECalPhysical%d%d", i, j);
      // Placement implemented by taking first two blocks and then skipping down by HCal lengths
      new G4PVPlacement(0, G4ThreeVector( x0 + i_factor*HCal_X,  y0 - j_factor*HCal_Y, 0),
               ECalLV[i][j], nameHolder, WorldLV, false, 0, fCheckOverlaps);
      G4cout<<"("<<i<<", "<<j<<"): "<<"("<<x0 + i_factor*HCal_X<<", "<<y0 - j_factor*HCal_Y<<")"<<G4endl;

    }
  }

  // Every 2x2 blocks has glue in the middle

  /* Horizontal glue between ECal blocks 
    □ □
    g g
    □ □
  */
  G4LogicalVolume* ECal_HorizGlueLV[NumECalBlocks][NumECalBlocks/2]; 
  G4VSolid* ECal_HorizGlueS = new G4Box("ECal_HorizGlueSolid", ECal_X/2., ECal_Glue_XY/2., ECal_Thickness/2.);

  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks/2; j ++)
    {
      sprintf(nameHolder, "ECal_HorizGlueLogical%d%d", i, j);
      ECal_HorizGlueLV[i][j] = new G4LogicalVolume(ECal_HorizGlueS, ActiveMaterial, nameHolder);
      G4double x0 = -2.*HCal_X + ECal_X/2. + Clearance_Gap;
      if(i % 2 != 0) x0 += (ECal_X + ECal_Glue_XY);
      G4double y0 = 2.*HCal_Y - ECal_Y - Clearance_Gap - ECal_Glue_XY/2.;
      G4int i_factor = i/2;
      sprintf(nameHolder, "ECal_HorizGluePhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector(x0 + i_factor*HCal_X, y0 - j*HCal_Y, 0), ECal_HorizGlueLV[i][j], nameHolder, WorldLV, false, 0, fCheckOverlaps);
    } 
  }

    /* Vertical glue between ECal blocks 
    glue running down middle of 2x2 blocks
    □ g □
      g
    □ g □
  */

  G4LogicalVolume* ECal_VertGlueLV[NumECalBlocks/2][NumECalBlocks/2];
  G4VSolid* ECal_VertGlueS = new G4Box("ECal_VertGlue", ECal_Glue_XY/2., (2*ECal_Y + ECal_Glue_XY)/2., ECal_Thickness/2.);

  for(G4int i = 0; i < NumECalBlocks/2; i++)
  {
    for(G4int j = 0; j < NumECalBlocks/2; j++)
    {
      sprintf(nameHolder, "ECal_VertGlueLogical%d%d", i, j);
      ECal_VertGlueLV[i][j] = new G4LogicalVolume(ECal_VertGlueS, ActiveMaterial, nameHolder);
      G4double x0 = -2.*HCal_X + ECal_X + Clearance_Gap + ECal_Glue_XY/2.;
      G4double y0 = 2.*HCal_Y - ECal_Y - Clearance_Gap - ECal_Glue_XY/2.;
      sprintf(nameHolder, "ECal_VertGluePhysical%d%d", i, j);
      new G4PVPlacement(0, G4ThreeVector(x0 + i*HCal_X, y0 - j*HCal_Y, 0), ECal_VertGlueLV[i][j], nameHolder, WorldLV, false, 0, fCheckOverlaps);
    } 
  }

  // Fibers
  // Non-sensitive cladding around each fiber
  // Cladding is 3% the total fiber thickness, i.e. Thickness = .03*fiber diameter
  // See documentation for exact fiber placement

  G4LogicalVolume* ECal_FiberCladdingLV[NumECalBlocks][NumECalBlocks];
  G4VSolid* ECal_FiberCladdingS = new G4Tubs("ECal_FiberCladdingSolid", ECal_Fiber_r - .03*ECal_Fiber_r, ECal_Fiber_r, ECal_Thickness/2., 0.0, 360.0*deg);

  // Fiber core
  G4LogicalVolume* ECal_FiberLV[NumECalBlocks][NumECalBlocks];
  G4VSolid* ECal_FiberS = new G4Tubs("ECal_FiberSolid", 0.0, ECal_Fiber_r - .03*ECal_Fiber_r, ECal_Thickness/2., 0.0, 360.0*deg);
  
  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      sprintf(nameHolder, "ECal_FiberCladdingLogical%d%d", i, j);
      ECal_FiberCladdingLV[i][j] = new G4LogicalVolume(ECal_FiberCladdingS, ActiveMaterial, nameHolder);
      sprintf(nameHolder, "ECal_FiberLogical%d%d", i, j);
      ECal_FiberLV[i][j] = new G4LogicalVolume(ECal_FiberS, ActiveMaterial, nameHolder);
      G4int num_fibers_block = 0;

      for(G4int fiber_i = 0; fiber_i < ECal_Fiber_Rows; fiber_i++)
      {
        G4double fiber_xspacing = 0.95865;
        G4double x0 = ECal_X/2. - 0.23966;
        if(fiber_i % 2 != 0) x0 -= fiber_xspacing/2.;
        G4double y0 = (ECal_Y/2. - 0.46) - fiber_i*0.820;

        for(G4int fiber_j = 0; fiber_j < ECal_Fiber_Cols; fiber_j++)
        {
          sprintf(nameHolder, "ECal_FiberCladdingPhysical%02d%02d%02d%02d", i, j, fiber_i, fiber_j);
          new G4PVPlacement(0, G4ThreeVector(x0 - fiber_j*fiber_xspacing, y0, 0), ECal_FiberCladdingLV[i][j], nameHolder, ECalLV[i][j], false, 0, false);

          sprintf(nameHolder, "ECal_FiberPhysical%02d%02d%02d%02d", i, j, fiber_i, fiber_j);
          new G4PVPlacement(0, G4ThreeVector(x0 - fiber_j*fiber_xspacing, y0, 0), ECal_FiberLV[i][j], nameHolder, ECalLV[i][j], false, 0, false);
          num_fibers_block++;

          // if(i == 4 && j == 3 && fiber_i == ECal_Fiber_Rows-1) std::cout<<fiber_j<<": ("<<(x0 - fiber_j*fiber_xspacing )<<", "<<y0<<")"<<std::endl;
        }
      }
      G4cout<<"Number of fibers in ECal block ("<<i<<", "<<j<<"): "<<num_fibers_block<<G4endl;
    }
  }

  G4cout<<"Finished Geometry construction."<<G4endl;
            
  //Visualization attributes
  G4VisAttributes invis=G4VisAttributes::Invisible;
  G4VisAttributes* RedVisAtt= new G4VisAttributes(G4Colour(1,0,0));//red
  G4VisAttributes* BlueVisAtt= new G4VisAttributes(G4Colour(0,0,1));//blue  
  G4VisAttributes* GreenVisAtt = new G4VisAttributes(G4Colour::Green());
  G4VisAttributes* GrayVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));//gray
  G4VisAttributes* CyanVisAtt = new G4VisAttributes(G4Colour::Cyan());
  G4VisAttributes* MagentaVisAtt = new G4VisAttributes(G4Colour::Magenta());
  
  RedVisAtt->SetVisibility(true);
  BlueVisAtt->SetVisibility(true);
  GreenVisAtt->SetVisibility(true);
  GrayVisAtt->SetVisibility(true);
  CyanVisAtt->SetVisibility(true);
  MagentaVisAtt->SetVisibility(true);
  MagentaVisAtt->SetForceSolid(true);
  // BlueVisAtt->SetForceSolid(true);
  
  WorldLV->SetVisAttributes(G4VisAttributes::Invisible);

  // GreenVisAtt->SetForceSolid(true);

  // Only HCal towers, ECal blocks, ECal glue, and one block's fibers are drawn
  // Change invis to other colors to see them
  // Warning: Drawing all the fibers slows down the visualization a lot
  for(int i=0;i<NumHCalTowers;i++)
  {
      for(int j=0;j<NumHCalTowers;j++)
      {
        HCalLV[i][j]->SetVisAttributes(RedVisAtt);
        HCalLayerHolderLV[i][j]->SetVisAttributes(GrayVisAtt);
    
        if( i == 0 && j == 0) 
        {
          HCalActiveLV[i][j]->SetVisAttributes(invis);
          HCalAbsorberLV[i][j]->SetVisAttributes(invis);
          HCalLayerLV[i][j]->SetVisAttributes(invis);
        }
        else if( i == 0 && j == 1)
        {
          HCalActiveLV[i][j]->SetVisAttributes(invis);
          HCalAbsorberLV[i][j]->SetVisAttributes(invis);
          HCalLayerLV[i][j]->SetVisAttributes(invis);
        }
        else if (i == 0 && j == 2)
        {
          HCalActiveLV[i][j]->SetVisAttributes(invis);
          HCalAbsorberLV[i][j]->SetVisAttributes(invis);
          HCalLayerLV[i][j]->SetVisAttributes(invis);
        }
        else
        {
          HCalActiveLV[i][j]->SetVisAttributes(invis);
          HCalAbsorberLV[i][j]->SetVisAttributes(invis);
          HCalLayerLV[i][j]->SetVisAttributes(invis);
        }
        if(i != NumHCalTowers-1) HCalWLS_LV[i][j]->SetVisAttributes(invis);
        if(j != NumHCalTowers - 1) HCalSteelLV[i][j]->SetVisAttributes(invis);
      }
      
  }

  //ECalLV->SetVisAttributes(BlueVisAtt);

  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      ECalLV[i][j]->SetVisAttributes(BlueVisAtt);
      if(j <4) ECal_HorizGlueLV[i][j]->SetVisAttributes(invis);
      if(i < 4 && j < 4) ECal_VertGlueLV[i][j]->SetVisAttributes(invis);

      if(i >4 && j > 4)
      {
        ECal_FiberCladdingLV[i][j]->SetVisAttributes(invis);
        ECal_FiberLV[i][j]->SetVisAttributes(invis);
      }
      else
      {
        ECal_FiberCladdingLV[i][j]->SetVisAttributes(invis);
        ECal_FiberLV[i][j]->SetVisAttributes(invis);
      }
    }
  }

  // Always return the physical World
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(0);

  char HitsNameHolder[200];
  char SDNameHolder[200];
  char DetectorNameHolder[200];

  CalorimeterSD* HCalSD[NumHCalTowers][NumHCalTowers];
    
  for(G4int i = 0; i < NumHCalTowers; i++)
  {
    for(G4int j = 0; j < NumHCalTowers; j++)
    { 
      sprintf(SDNameHolder, "HCalSD%d%d", i, j);
      sprintf(HitsNameHolder, "HCalHitsCollection%d%d", i, j);
      sprintf(DetectorNameHolder, "HCalActiveLogical%d%d", i, j);

      HCalSD[i][j] = new CalorimeterSD(SDNameHolder, HitsNameHolder, NumHCalLayers);
      G4SDManager::GetSDMpointer()->AddNewDetector(HCalSD[i][j]);
      SetSensitiveDetector(DetectorNameHolder, HCalSD[i][j]);
    }
  }

  //CalorimeterSD* ECalSD = new CalorimeterSD("ECalSD", "ECalHitsCollection", 1);
  //G4SDManager::GetSDMpointer()->AddNewDetector(ECalSD);
  //SetSensitiveDetector("ECalLogical", ECalSD);

  CalorimeterSD* ECalSD[NumECalBlocks][NumECalBlocks];
  for(G4int i = 0; i < NumECalBlocks; i++)
  {
    for(G4int j = 0; j < NumECalBlocks; j++)
    {
      sprintf(SDNameHolder, "ECalSD%d%d", i, j);
      sprintf(HitsNameHolder, "ECalHitsCollection%d%d", i, j);
      sprintf(DetectorNameHolder, "ECal_FiberLogical%d%d", i, j);

      ECalSD[i][j] = new CalorimeterSD(SDNameHolder, HitsNameHolder, 1);
      G4SDManager::GetSDMpointer()->AddNewDetector(ECalSD[i][j]);
      SetSensitiveDetector(DetectorNameHolder, ECalSD[i][j]);
    }
  }

  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(0);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
