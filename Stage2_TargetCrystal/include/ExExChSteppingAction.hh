#ifndef ExExChSteppingAction_h
#define ExExChSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "ExExChRunAction.hh"

class G4LogicalVolume;
class ExExChRunAction;

class ExExChSteppingAction : public G4UserSteppingAction
{
public:
    ExExChSteppingAction(ExExChRunAction*);
    virtual ~ExExChSteppingAction();
    virtual void UserSteppingAction(const G4Step *);

private:
    ExExChRunAction* runAction;
};

#endif
