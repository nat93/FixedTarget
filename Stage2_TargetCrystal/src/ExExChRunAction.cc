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

#include "ExExChRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "ExExChAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChRunAction::ExExChRunAction(): G4UserRunAction()
{
    _outputFileName = "targetcrystal.root";
    _timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChRunAction::~ExExChRunAction()
{
    delete _timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    _timer->Start();

    tree = new TTree("Tree", "A Tree with Initial particles distributions + LambdaC new zprod, After Target and Crystal");

    tree->Branch("Nparticle",    &EventArray.Nparticle,  "Nparticle/I");
    tree->Branch("Nproj",        &EventArray.Nproj,      "Nproj/I");
    tree->Branch("Ntarg",        &EventArray.Ntarg,      "Ntarg/I");
    tree->Branch("Ncoll",        &EventArray.Ncoll,      "Ncoll/I");
    tree->Branch("index",        EventArray.index,       "index[Nparticle]/F");
    tree->Branch("status",       EventArray.status,      "status[Nparticle]/F");
    tree->Branch("IsFinal",      EventArray.IsFinal,     "IsFinal[Nparticle]/O");
    tree->Branch("Id",           EventArray.ID,          "Id[Nparticle]/F");
    tree->Branch("firstDau",     EventArray.firstDau,    "firstDau[Nparticle]/F");
    tree->Branch("lastDau",      EventArray.lastDau,     "lastDau[Nparticle]/F");
    tree->Branch("Mindex",       EventArray.Mindex,      "Mindex[Nparticle]/F");
    tree->Branch("xprod",        EventArray.xprod,       "xprod[Nparticle]/F");
    tree->Branch("yprod",        EventArray.yprod,       "yprod[Nparticle]/F");
    tree->Branch("zprod",        EventArray.zprod,       "zprod[Nparticle]/F");
    tree->Branch("xdecay",       EventArray.xdecay,      "xdecay[Nparticle]/F");
    tree->Branch("ydecay",       EventArray.ydecay,      "ydecay[Nparticle]/F");
    tree->Branch("zdecay",       EventArray.zdecay,      "zdecay[Nparticle]/F");
    tree->Branch("Px",           EventArray.Px,          "Px[Nparticle]/F");
    tree->Branch("Py",           EventArray.Py,          "Py[Nparticle]/F");
    tree->Branch("Pz",           EventArray.Pz,          "Pz[Nparticle]/F");
    tree->Branch("isGenerated",  &EventArray.isGenerated,"isGenerated/O");
    tree->Branch("Px_Lc_0",      &EventArray.Px_Lc_0,    "Px_Lc_0/F");
    tree->Branch("Py_Lc_0",      &EventArray.Py_Lc_0,    "Py_Lc_0/F");
    tree->Branch("Pz_Lc_0",      &EventArray.Pz_Lc_0,    "Pz_Lc_0/F");
    tree->Branch("isTarget",     &EventArray.isTarget,   "isTarget/O");
    tree->Branch("Px_Lc_1",      &EventArray.Px_Lc_1,    "Px_Lc_1/F");
    tree->Branch("Py_Lc_1",      &EventArray.Py_Lc_1,    "Py_Lc_1/F");
    tree->Branch("Pz_Lc_1",      &EventArray.Pz_Lc_1,    "Pz_Lc_1/F");
    tree->Branch("isCrystal",    &EventArray.isCrystal,  "isCrystal/O");
    tree->Branch("Px_Lc_2",      &EventArray.Px_Lc_2,    "Px_Lc_2/F");
    tree->Branch("Py_Lc_2",      &EventArray.Py_Lc_2,    "Py_Lc_2/F");
    tree->Branch("Pz_Lc_2",      &EventArray.Pz_Lc_2,    "Pz_Lc_2/F");
    tree->Branch("E",            EventArray.E,           "E[Nparticle]/F");
    tree->Branch("mass",         EventArray.mass,        "mass[Nparticle]/F");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChRunAction::EndOfRunAction(const G4Run* /*run*/)
{
    G4cout<<G4endl<<"--> Output file name: "<<_outputFileName<<G4endl;
    TFile* _hfile = new TFile(_outputFileName, "RECREATE");
    if(_hfile->IsZombie()) exit(-1);
    tree->Write();
    _hfile->Close();

    _timer->Stop();

    delete tree;

    delete _hfile;
    G4cout<<G4endl<<G4endl<<"Time: "<<*_timer<<G4endl<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
