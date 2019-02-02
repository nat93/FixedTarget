//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "ExExChDetectorConstruction.hh"

#include "ExExChUserActionInitialization.hh"

#include "ExExChPrimaryGeneratorAction.hh"
#include "ExExChTrackingAction.hh"
#include "ExExChStackingAction.hh"
#include "ExExChEventAction.hh"
#include "ExExChRunAction.hh"
#include "ExExChSteppingAction.hh"


#include "ExExChPhysicsList.hh"
#include "QGSP_BERT.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
//    "/media/andrii/F492773C92770302/SPS_SimulationData/hardccbar.root"
    if(argc != 5)
    {
        cout<<"ERROR::Not enough input arguments!"<<endl;
        cout<<"--> [0]: script_name"<<endl;
        cout<<"--> [1]: mac_file_name"<<endl;
        cout<<"--> [2]: input_file_name"<<endl;
        cout<<"--> [3]: output_file_name"<<endl;
        cout<<"--> [4]: event_id_initial"<<endl;
        assert(0);
    }

    // Construct the default run manager

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(NULL));

    G4RunManager* runManager = new G4RunManager;
    G4cout << "MT MODE OFF" << G4endl;
    
    // Activate UI-command base scorer
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);
    
    
    // Set mandatory initialization classes
    G4VUserDetectorConstruction* detector = new ExExChDetectorConstruction();
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new ExExChPhysicsList());    
    runManager->SetUserAction(new ExExChStackingAction());
    runManager->SetUserAction(new ExExChTrackingAction());

    ExExChRunAction* runAction = new ExExChRunAction();
    G4String rootFileName = argv[3];
    runAction->SetOutputFileName(rootFileName);
    G4String inputFileName = argv[2];
    ExExChPrimaryGeneratorAction* primary = new ExExChPrimaryGeneratorAction(atol(argv[4]),runAction,inputFileName);
    runManager->SetUserAction(primary);
    ExExChEventAction* eventAction = new ExExChEventAction(runAction);
    ExExChSteppingAction* stepAction = new ExExChSteppingAction(runAction);
    runManager->SetUserAction(eventAction);
    runManager->SetUserAction(stepAction);
    runManager->SetUserAction(runAction);

    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();
    G4String scriptName = argv[1];

    if(scriptName.contains("vis"))
    {
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();

        G4String command = "/control/execute ";
        UI->ApplyCommand(command+scriptName);

        delete visManager;
    }
    else
    {
        // Batch mode
        G4String command = "/control/execute ";
        UI->ApplyCommand(command+scriptName);
    }
    // Job termination
    delete runManager;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
