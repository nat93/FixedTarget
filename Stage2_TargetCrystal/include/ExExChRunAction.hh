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

#ifndef ExExChRunAction_h
#define ExExChRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4Timer.hh"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

using namespace std;

class G4Run;

class ExExChRunAction : public G4UserRunAction
{
public:
    ExExChRunAction();
    virtual ~ExExChRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void SetOutputFileName(G4String fileName) {_outputFileName = fileName;}
    G4String GetOutputFileName() { return _outputFileName;}

    TTree* tree;

    static const int arraySize = 20000;
    typedef struct
    {
        Float_t index[arraySize];
        Float_t status[arraySize];
        Bool_t IsFinal[arraySize];
        Float_t ID[arraySize];
        Float_t firstDau[arraySize];
        Float_t lastDau[arraySize];
        Float_t Mindex[arraySize];
        Float_t xprod[arraySize];
        Float_t yprod[arraySize];
        Float_t zprod[arraySize];
        Float_t xdecay[arraySize];
        Float_t ydecay[arraySize];
        Float_t zdecay[arraySize];
        Float_t Px[arraySize];
        Float_t Py[arraySize];
        Float_t Pz[arraySize];
        Bool_t isGenerated;
        Float_t Px_Lc_0;
        Float_t Py_Lc_0;
        Float_t Pz_Lc_0;
        Bool_t isTarget;
        Float_t Px_Lc_1;
        Float_t Py_Lc_1;
        Float_t Pz_Lc_1;
        Bool_t isCrystal;
        Float_t Px_Lc_2;
        Float_t Py_Lc_2;
        Float_t Pz_Lc_2;
        Float_t E[arraySize];
        Float_t mass[arraySize];
        Int_t Nparticle;
        Int_t Nproj;
        Int_t Ntarg;
        Int_t Ncoll;
    } OUTPUTEVENT;

    OUTPUTEVENT EventArray;

private:
    G4Timer* _timer;
    G4String _outputFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

