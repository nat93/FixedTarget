#include "ExExChSteppingAction.hh"
#include "ExExChDetectorConstruction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4Transform3D.hh"
#include "TMath.h"

ExExChSteppingAction::ExExChSteppingAction(ExExChRunAction* runAct) : runAction(runAct)
{}

ExExChSteppingAction::~ExExChSteppingAction()
{}

void ExExChSteppingAction::UserSteppingAction(const G4Step *step)
{
    G4Track* track = step->GetTrack();

    G4StepPoint*            prePoint  = step->GetPreStepPoint();
    G4VPhysicalVolume*      prePV     = prePoint->GetPhysicalVolume();
    G4ThreeVector           p         = track->GetMomentum();

    if(prePoint->GetStepStatus() == fGeomBoundary && prePV->GetName().contains("phantom_0"))
    {
        runAction->EventArray.isTarget  = true;
        runAction->EventArray.Px_Lc_1  = p.x()/GeV;
        runAction->EventArray.Py_Lc_1  = p.y()/GeV;
        runAction->EventArray.Pz_Lc_1  = p.z()/GeV;
    }

    if(prePoint->GetStepStatus() == fGeomBoundary && prePV->GetName().contains("phantom_1"))
    {
        runAction->EventArray.isCrystal  = true;
        runAction->EventArray.Px_Lc_2  = p.x()/GeV;
        runAction->EventArray.Py_Lc_2  = p.y()/GeV;
        runAction->EventArray.Pz_Lc_2  = p.z()/GeV;
    }
}
