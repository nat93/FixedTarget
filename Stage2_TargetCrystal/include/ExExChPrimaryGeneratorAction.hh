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
// --------------------------------------------------------------
//

#ifndef ExExChPrimaryGeneratorAction_h
#define ExExChPrimaryGeneratorAction_h 1

// C++
#include "iostream"
#include "string"
#include "fstream"
#include "vector"
#include "stdlib.h"
#include "time.h"
#include "iomanip"
#include "assert.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "unistd.h"
#include "sys/time.h"
#include "sstream"

// ROOT
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph2DErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TRandom3.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TVector3.h"

// GEANT4
#include "G4VUserPrimaryGeneratorAction.hh"
#include "ExExChRunAction.hh"
#include "globals.hh"

using namespace std;

class G4GeneralParticleSource;
class G4Event;
class ExExChDetectorConstruction;
class G4ParticleGun;
class ExExChRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExExChPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    ExExChPrimaryGeneratorAction(Long64_t eventIDini, ExExChRunAction*, TString input_file_name);
    virtual ~ExExChPrimaryGeneratorAction();
    
    void GeneratePrimaries(G4Event*);
    void extractIntegerWords(string str);
    G4bool isNotDecayInTarget(G4double l, G4double  tau, G4double beta, G4double gamma, G4double &lInTarget);
    G4bool isNotDecayInCrystal(G4double l, G4double  tau, G4double beta, G4double gamma);
    G4double notDecayRate(G4double l, G4double N0, G4double tau, G4double beta, G4double gamma);
    G4double getParticleFlightTime(G4double l, G4double beta);
    void getThetaAndPhi(G4float px, G4float py, G4float pz, G4double &theta, G4double &phi, G4double &mag);
    G4double getGamma(G4double mass, G4double pmag);
    G4double getBeta(G4double gamma);

private:
    G4ParticleGun* fParticleGun;
    ExExChRunAction* runAction;
    TRandom3* _rnd;

    TChain* fChain;

    static const int arraySize = 20000;
    G4float _index[arraySize];
    G4float _status[arraySize];
    G4bool _IsFinal[arraySize];
    G4float _ID[arraySize];
    G4float _firstDau[arraySize];
    G4float _lastDau[arraySize];
    G4float _Mindex[arraySize];
    G4float _xprod[arraySize];
    G4float _yprod[arraySize];
    G4float _zprod[arraySize];
    G4float _xdecay[arraySize];
    G4float _ydecay[arraySize];
    G4float _zdecay[arraySize];
    G4float _Px[arraySize];
    G4float _Py[arraySize];
    G4float _Pz[arraySize];
    G4float _E[arraySize];
    G4float _mass[arraySize];
    G4int _Nparticle;
    G4int _Nproj;
    G4int _Ntarg;
    G4int _Ncoll;
    Long64_t nEntries;
    Long64_t iEntry;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


