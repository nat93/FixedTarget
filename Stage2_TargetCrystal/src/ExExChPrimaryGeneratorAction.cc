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

#include "ExExChPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ParticleGun.hh"
#include "TMath.h"
#include "G4GeneralParticleSource.hh"

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

#include "../include/constants.hh"

using namespace constants;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChPrimaryGeneratorAction::ExExChPrimaryGeneratorAction(Long64_t eventIDini, ExExChRunAction* runAct, TString input_file_name) : runAction(runAct)
{
    fParticleGun = new G4ParticleGun(1);
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
    fParticleGun->SetParticleMomentum(400.0*CLHEP::GeV);

    fChain = new TChain("Tree");
    fChain->Add(input_file_name);

    fChain->SetBranchAddress("Nparticle",    &_Nparticle);
    fChain->SetBranchAddress("Nproj",        &_Nproj);
    fChain->SetBranchAddress("Ntarg",        &_Ntarg);
    fChain->SetBranchAddress("Ncoll",        &_Ncoll);
    fChain->SetBranchAddress("index",        _index);
    fChain->SetBranchAddress("status",       _status);
    fChain->SetBranchAddress("IsFinal",      _IsFinal);
    fChain->SetBranchAddress("Id",           _ID);
    fChain->SetBranchAddress("firstDau",     _firstDau);
    fChain->SetBranchAddress("lastDau",      _lastDau);
    fChain->SetBranchAddress("Mindex",       _Mindex);
    fChain->SetBranchAddress("xprod",        _xprod);
    fChain->SetBranchAddress("yprod",        _yprod);
    fChain->SetBranchAddress("zprod",        _zprod);
    fChain->SetBranchAddress("xdecay",       _xdecay);
    fChain->SetBranchAddress("ydecay",       _ydecay);
    fChain->SetBranchAddress("zdecay",       _zdecay);
    fChain->SetBranchAddress("Px",           _Px);
    fChain->SetBranchAddress("Py",           _Py);
    fChain->SetBranchAddress("Pz",           _Pz);
    fChain->SetBranchAddress("E",            _E);
    fChain->SetBranchAddress("mass",         _mass);

    G4cout<<"--> Input file: "<<input_file_name<<G4endl;
    nEntries = fChain->GetEntries();
    G4cout<<"--> nEntries: "<<nEntries<<G4endl;

    iEntry = eventIDini;

    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    _rnd = new TRandom3(seed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChPrimaryGeneratorAction::~ExExChPrimaryGeneratorAction(){
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    runAction->EventArray.isGenerated  = false;

    if(iEntry < nEntries)
    {
        fChain->GetEntry(iEntry);

        for(Int_t i = 0; i < _Nparticle; i++)
        {
            if(_ID[i] == 4122) // Lambda_c+ = 4122
            {
                G4double theta, phi, pmag;
                getThetaAndPhi(_Px[i],_Py[i],_Pz[i],theta,phi,pmag);
                G4double gamma = getGamma(_mass[i],pmag);
                G4double beta = getBeta(gamma);
                G4double lInTarget = -999.999;

                if(isNotDecayInTarget(targetLength,tau_LambdaC,beta,gamma,lInTarget))
                {
                    if(isNotDecayInCrystal(crystalLength,tau_LambdaC,beta,gamma))
                    {
                        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
                        G4ParticleDefinition* particle;

                        //** Proton **//
//                        particle = particleTable->FindParticle(2212);

                        //** Lambda_c+ **//
                        particle = particleTable->FindParticle(4122);
                        particle->SetPDGStable(true);

                        if(particle)
                        {
                            _zprod[i] = 0.0 - lInTarget - constants::targetCrystalGap - constants::crystalLength/2.0;
                            G4double mass = particle->GetPDGMass()/CLHEP::GeV;
                            G4double Ekin   = (TMath::Sqrt(pmag*pmag + mass*mass) - mass);

                            fParticleGun->SetParticleDefinition(particle);
                            fParticleGun->SetParticleEnergy(Ekin*CLHEP::GeV);
                            fParticleGun->SetParticleMomentumDirection(G4ThreeVector(_Px[i]/pmag,_Py[i]/pmag,_Pz[i]/pmag));
                            fParticleGun->SetParticlePosition(G4ThreeVector(
                                                                  _xprod[i]*CLHEP::mm,
                                                                  _yprod[i]*CLHEP::mm,
                                                                  _zprod[i]*CLHEP::mm));

                            runAction->EventArray.Nparticle   = _Nparticle;
                            runAction->EventArray.Nproj       = _Nproj;
                            runAction->EventArray.Ntarg       = _Ntarg;
                            runAction->EventArray.Ncoll       = _Ncoll;

                            for(Int_t j = 0; j < _Nparticle; j++)
                            {
                                runAction->EventArray.index[j]    = _index[j];
                                runAction->EventArray.status[j]   = _status[j];
                                runAction->EventArray.IsFinal[j]  = _IsFinal[j];
                                runAction->EventArray.ID[j]       = _ID[j];
                                runAction->EventArray.firstDau[j] = _firstDau[j];
                                runAction->EventArray.lastDau[j]  = _lastDau[j];
                                runAction->EventArray.Mindex[j]   = _Mindex[j];
                                runAction->EventArray.xprod[j]    = _xprod[j];
                                runAction->EventArray.yprod[j]    = _yprod[j];
                                runAction->EventArray.zprod[j]    = _zprod[j];
                                runAction->EventArray.xdecay[j]   = _xdecay[j];
                                runAction->EventArray.ydecay[j]   = _ydecay[j];
                                runAction->EventArray.zdecay[j]   = _zdecay[j];
                                runAction->EventArray.Px[j]       = _Px[j];
                                runAction->EventArray.Py[j]       = _Py[j];
                                runAction->EventArray.Pz[j]       = _Pz[j];
                                runAction->EventArray.E[j]        = _E[j];
                                runAction->EventArray.mass[j]     = _mass[j];
                            }

                            runAction->EventArray.isGenerated  = true;
                            runAction->EventArray.Px_Lc_0  = _Px[i];
                            runAction->EventArray.Py_Lc_0  = _Py[i];
                            runAction->EventArray.Pz_Lc_0  = _Pz[i];

                            runAction->EventArray.isTarget  = false;
                            runAction->EventArray.Px_Lc_1  = -999;
                            runAction->EventArray.Py_Lc_1  = -999;
                            runAction->EventArray.Pz_Lc_1  = -999;

                            runAction->EventArray.isCrystal  = false;
                            runAction->EventArray.Px_Lc_2  = -999;
                            runAction->EventArray.Py_Lc_2  = -999;
                            runAction->EventArray.Pz_Lc_2  = -999;

                            fParticleGun->GeneratePrimaryVertex(anEvent);
                        }
                    }
                }
            }
        }

        iEntry++;
    }
    else
    {
        G4cout<<" EMPTY EVENT !!! "<<iEntry<<"/"<<nEntries<<G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ExExChPrimaryGeneratorAction::isNotDecayInTarget(G4double l, G4double  tau, G4double beta, G4double gamma, G4double &lInTarget)
{
    G4double N0 = 1.0;
    lInTarget = l - _rnd->Uniform(0.0,l);
    G4double probNotDec = notDecayRate(lInTarget, N0, tau, beta, gamma);

    if(probNotDec > _rnd->Uniform(0.0,1.0))
    {
        return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ExExChPrimaryGeneratorAction::isNotDecayInCrystal(G4double l, G4double  tau, G4double beta, G4double gamma)
{
    G4double N0 = 1.0;

    G4double probNotDec = notDecayRate(l, N0, tau, beta, gamma);

    if(probNotDec > _rnd->Uniform(0.0,1.0))
    {
        return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChPrimaryGeneratorAction::notDecayRate(G4double l, G4double N0, G4double tau, G4double beta, G4double gamma)
{
    G4double t;
    if(l >= 0.0 && N0 >= 0.0 && tau >= 0.0 && beta >= 0.0 && gamma >= 1.0)
    {
        t = getParticleFlightTime(l, beta);
        if(t >= 0.0)
        {
            return N0*TMath::Exp(-t/(tau*gamma));
        }
        else
        {
            cout<<" ERROR ---> t < 0"<<endl
               <<" t = "<<t<<endl;
            assert(0);
        }
    }
    else
    {
        cout<<" ERROR: "<<endl<<" l = "
           <<l<<endl
          <<" beta = "<<beta<<endl
         <<" N0 = "<<N0<<endl
        <<" tau = "<<tau<<endl
        <<" gamma = "<<gamma<<endl;
        assert(0);
    }
    return -999.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChPrimaryGeneratorAction::getParticleFlightTime(G4double l, G4double beta)
{
    if(l >= 0.0 && beta > 0.0)
    {
        return l/(c_lightv*beta);
    }
    else
    {
        cout<<" ERROR ---> l or beta < 0"<<endl
           <<" l = "<<l<<endl
          <<" beta = "<<beta<<endl;
        assert(0);
    }
    return -999.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPrimaryGeneratorAction::getThetaAndPhi(G4float px, G4float py, G4float pz, G4double &theta, G4double &phi, G4double &mag)
{
    TVector3 *pv = new TVector3(px, py, pz);
    theta   = pv->Theta();
    phi     = pv->Phi();
    mag     = pv->Mag();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChPrimaryGeneratorAction::getGamma(G4double mass, G4double pmag)
{
    if(mass > 0.0)
    {
        return TMath::Sqrt(1.0 + pmag*pmag/(mass*mass));
    }
    else
    {
        cout<<" ERROR ---> mass <= 0.0"<<endl
           <<" mass = "<<mass<<endl;
        assert(0);
    }
    return -999.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChPrimaryGeneratorAction::getBeta(G4double gamma)
{
    if(gamma >= 1.0)
    {
        return TMath::Sqrt(1.0 - 1.0/(gamma*gamma));
    }
    else
    {
        cout<<" ERROR ---> gamma < 1"<<endl
           <<" gamma = "<<gamma<<endl;
        assert(0);
    }
    return -999.0;
}
