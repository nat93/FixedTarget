//////////////////////////////////////////////////////////////////////
///                Date: 26/01/2019                                ///
///               Autor: Andrii Natochii                           ///
///       Initial Autor: Leonid Burmistrov 9.09.15                 ///
///       Initial Autor: Laure Massacrier 16.01.14                 ///
/// Program description: It illustrates how to generate and        ///
///                      analyze Lambda_c events, inside root      ///
///                      framework with pythia                     ///
///                      (initial version 8185). It produces a     ///
///                      tree with an entry per event. Each event  ///
///                      contains an array of all the particles    ///
///                      produced in the event                     ///
///                      (including the initial protons).          ///
//////////////////////////////////////////////////////////////////////

//C, C++
#include <iostream>
#include <fstream>
#include <assert.h>

//pythia
#include "Pythia8/Pythia.h"

//root
#include "TROOT.h"
#include "TH1.h"
#include "TH2D.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"

//each entry is an event, particles are stored in an array. Array size is the maximum number of particles an array can have

const int arraySize = 20000;
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
    Float_t E[arraySize];
    Float_t mass[arraySize];
    Int_t Nparticle;
    Int_t Nproj;
    Int_t Ntarg;
    Int_t Ncoll;
} OUTPUTEVENT;
static OUTPUTEVENT EventArray; // Event array for pythia output

using namespace Pythia8;
using namespace std;

int main(int argc, char* argv[])
{
    Int_t rndSeed;
    Float_t protonE;
    Int_t nevents;
    TString outRootFileName;

    if(argc == 5)
    {
        rndSeed         = atoi(argv[1]);
        protonE         = atof(argv[2]);
        nevents         = atoi(argv[3]);
        outRootFileName = argv[4];

        cout<<endl;
        cout<<"--> Start <--"<<endl;
        cout<<"--> Lambda_c generation with p p colision "<<endl;
        cout<<"--> random seed      : "<<rndSeed<<endl;
        cout<<"--> proton energy    : "<<protonE<<" GeV/c"<<endl;
        cout<<"--> number of events : "<<nevents<<endl;
        cout<<"--> Out root file    : "<<outRootFileName<<endl;
        cout<<endl;
    }
    else
    {
        cout<<endl;
        cout<<"--> ERROR:: Wrong number in the input arguments!"<<endl;
        cout<<"--> [1] -- random seed"<<endl;
        cout<<"--> [2] -- proton energy"<<endl;
        cout<<"--> [3] -- number of events"<<endl;
        cout<<"--> [4] -- out root file"<<endl;
        cout<<endl;
        assert(0);
    }

    // Generator. Shorthand for the event.
    Pythia pythia;
    Event& event = pythia.event;

    pythia.readString("Random:setSeed = on");
    TString rndSeedStrRoot = "Random:seed = ";
    rndSeedStrRoot += rndSeed;
    pythia.readString(rndSeedStrRoot.Data());

    TString numberOfEvents = "Main:numberOfEvents = ";
    numberOfEvents += nevents;
    pythia.readString(numberOfEvents.Data());

    // Read in commands from external file.
    pythia.readFile("hardccbar.cmnd");

    // Extract settings to be used in the main program.
    int    nEvent    = pythia.mode("Main:numberOfEvents");
    int    nAbort    = pythia.mode("Main:timesAllowErrors");
    bool   showCS    = pythia.flag("Main:showChangedSettings");
    bool   showAS    = pythia.flag("Main:showAllSettings");
    bool   showCPD   = pythia.flag("Main:showChangedParticleData");
    bool   showAPD   = pythia.flag("Main:showAllParticleData");
    bool   showAStat = pythia.flag("Main:showAllStatistics");
    cout<<"showAStat "<<showAStat<<endl;

    // We are interested in Lc -> p K- pi+
    bool isLambdacgenerated;
    bool isLambdacProtongenerated;
    bool isLambdacKaonMgenerated;
    bool isLambdacPionPgenerated;

    // Initialize. Beam parameters set in .cmnd file.
    //pythia.init(2212,2212,400,0); // 400 GeV/c proton on a proton at rest
    pythia.init(2212,2212,protonE,0);

    // List changed data.
    if (showCS) pythia.settings.listChanged();
    if (showAS) pythia.settings.listAll();

    // List particle data.
    if (showCPD) pythia.particleData.listChanged();
    if (showAPD) pythia.particleData.listAll();

    //create NTuple to store all the particles in the event
    //TNtuple *ntuple_particle = new TNtuple("particles","particles","evt:id:status:charge:px:py:pz:e:m:pT:mT:theta:phi:y:eta");

    //create a TTree
    TTree* outputTTree = new TTree("Tree","A tree with all particle parameters");
    outputTTree->Branch("Nparticle",    &EventArray.Nparticle,  "Nparticle/I");
    outputTTree->Branch("Nproj",        &EventArray.Nproj,      "Nproj/I");
    outputTTree->Branch("Ntarg",        &EventArray.Ntarg,      "Ntarg/I");
    outputTTree->Branch("Ncoll",        &EventArray.Ncoll,      "Ncoll/I");
    outputTTree->Branch("index",        EventArray.index,       "index[Nparticle]/F");
    outputTTree->Branch("status",       EventArray.status,      "status[Nparticle]/F");
    outputTTree->Branch("IsFinal",      EventArray.IsFinal,     "IsFinal[Nparticle]/O");
    outputTTree->Branch("Id",           EventArray.ID,          "Id[Nparticle]/F");
    outputTTree->Branch("firstDau",     EventArray.firstDau,    "firstDau[Nparticle]/F");
    outputTTree->Branch("lastDau",      EventArray.lastDau,     "lastDau[Nparticle]/F");
    outputTTree->Branch("Mindex",       EventArray.Mindex,      "Mindex[Nparticle]/F");
    outputTTree->Branch("xprod",        EventArray.xprod,       "xprod[Nparticle]/F");
    outputTTree->Branch("yprod",        EventArray.yprod,       "yprod[Nparticle]/F");
    outputTTree->Branch("zprod",        EventArray.zprod,       "zprod[Nparticle]/F");
    outputTTree->Branch("xdecay",       EventArray.xdecay,      "xdecay[Nparticle]/F");
    outputTTree->Branch("ydecay",       EventArray.ydecay,      "ydecay[Nparticle]/F");
    outputTTree->Branch("zdecay",       EventArray.zdecay,      "zdecay[Nparticle]/F");
    outputTTree->Branch("Px",           EventArray.Px,          "Px[Nparticle]/F");
    outputTTree->Branch("Py",           EventArray.Py,          "Py[Nparticle]/F");
    outputTTree->Branch("Pz",           EventArray.Pz,          "Pz[Nparticle]/F");
    outputTTree->Branch("E",            EventArray.E,           "E[Nparticle]/F");
    outputTTree->Branch("mass",         EventArray.mass,        "mass[Nparticle]/F");

    // number of fills to the output tree
    Long64_t nFills = 0;
    Int_t iAbort = 0;

    TFile *outFile = new TFile(outRootFileName.Data(),"RECREATE");

    Int_t LambdaC_ind, proton_ind, kaon_ind, pion_ind;

    //for(int iEvent = 0; iEvent < nEvent; iEvent++)

    while(nFills < nEvent)
    {
        if(!pythia.next())
        {
            if(++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
            break;
        }

        isLambdacgenerated          = kFALSE;
        isLambdacProtongenerated    = kFALSE;
        isLambdacKaonMgenerated     = kFALSE;
        isLambdacPionPgenerated     = kFALSE;
        // Nparticle is the total number of particles in the event (including the initial protons), the entry 0 is not counted (system)
        EventArray.Nparticle = event.size()-1;
        EventArray.Nproj = 1; // number of incoming particles
        EventArray.Ntarg = 1; // number of target particles
        EventArray.Ncoll = 1; // number of collitions

        //particle loop  (entry 0 is the all event! counting should start at 1)
        for (Int_t i = 1; i < event.size(); i++)
        {
            //fill variables for the tree :
            EventArray.index[i]     = i;
            EventArray.status[i]    = event[i].status();    // status code. The status code includes information on how a particle was produced
            EventArray.IsFinal[i]   = event[i].isFinal();   // true for a remaining particle, i.e. one with positive status code, else false
            EventArray.ID[i]        = event[i].id();        // the identity of a particle, according to the PDG particle codes
            EventArray.firstDau[i]  = event[i].daughter1(); // the indices in the event record where the first and last daughters are stored, if any
            EventArray.lastDau[i]   = event[i].daughter2(); // the indices in the event record where the first and last daughters are stored, if any
            EventArray.Mindex[i]    = event[i].mother1();   // the indices in the event record where the first and last mothers are stored, if any
            EventArray.xprod[i]     = event[i].xProd();     // the production vertex coordinate
            EventArray.yprod[i]     = event[i].yProd();     // the production vertex coordinate
            EventArray.zprod[i]     = event[i].zProd();     // the production vertex coordinate
            EventArray.xdecay[i]    = event[i].xDec();      // the decay vertex coordinate
            EventArray.ydecay[i]    = event[i].yDec();      // the decay vertex coordinate
            EventArray.zdecay[i]    = event[i].zDec();      // the decay vertex coordinate
            EventArray.Px[i]        = event[i].px();        // the particle four-momentum component
            EventArray.Py[i]        = event[i].py();        // the particle four-momentum component
            EventArray.Pz[i]        = event[i].pz();        // the particle four-momentum component
            EventArray.E[i]         = event[i].e();         // the particle four-momentum component
            EventArray.mass[i]      = event[i].m();         // the particle mass, stored with a minus sign (times the absolute value) for spacelike virtual particles

            //true if a id="4122" lambda_c has been generated
            if(event[i].id() == 4122)
            {
                isLambdacgenerated = kTRUE;
                LambdaC_ind = i;
            }
            //true if a id="2212" proton has been generated
            if(event[i].id() == 2212)
            {
                if(event[i].mother1() == LambdaC_ind)
                {
                    isLambdacProtongenerated = kTRUE;
                    proton_ind = i;
                }
            }
            //true if a id="-321" kaon- ( kaon-* )has been generated
            if(event[i].id() == -321 || event[i].id() == -323)
            {
                if(event[i].mother1() == LambdaC_ind)
                {
                    isLambdacKaonMgenerated = kTRUE;
                    kaon_ind = i;
                }
            }
            //true if a id="211" pion+ has been generated
            if(event[i].id() == 211)
            {
                if(event[i].mother1() == LambdaC_ind)
                {
                    isLambdacPionPgenerated = kTRUE;
                    pion_ind = i;
                }
            }
        }
        //end of particle loop

        //fill each event
        if(isLambdacgenerated && isLambdacProtongenerated && isLambdacKaonMgenerated && isLambdacPionPgenerated)
        {
            if(nFills%100 == 0)
            {
                cout<<endl<<endl;
                cout<<"--> nFills = "<<nFills<<endl;
                cout<<"--> Here is the info:"<<endl;
                cout<<"--> LambdaC_ind = "<<LambdaC_ind<<endl;
                cout<<"--> EventArray.firstDau["<<LambdaC_ind<<"] = "<<EventArray.firstDau[LambdaC_ind]<<endl;
                cout<<"--> EventArray.lastDau["<<LambdaC_ind<<"] = "<<EventArray.lastDau[LambdaC_ind]<<endl;
                cout<<"--> proton_ind = "<<proton_ind<<" | "<<EventArray.ID[proton_ind]<<endl;
                cout<<"--> EventArray.Mindex["<<proton_ind<<"] = "<<EventArray.Mindex[proton_ind]<<endl;
                cout<<"--> kaon_ind = "<<kaon_ind<<" | "<<EventArray.ID[kaon_ind]<<endl;
                cout<<"--> EventArray.Mindex["<<kaon_ind<<"] = "<<EventArray.Mindex[kaon_ind]<<endl;
                cout<<"--> pion_ind = "<<pion_ind<<" | "<<EventArray.ID[pion_ind]<<endl;
                cout<<"--> EventArray.Mindex["<<pion_ind<<"] = "<<EventArray.Mindex[pion_ind]<<endl;
                cout<<endl;
            }

            outputTTree->Fill();
            nFills++;
        }
    }
    //end of event loop

    //save TTree
    outputTTree->Write();
    outFile->Close();

    delete outFile;
    delete outputTTree;
    pythia.statistics();

    cout<<"--> nFills = "<<nFills<<endl;

    cout<<endl;
    cout<<"--> Stop <--"<<endl;
    cout<<"--> Lambda_c generation with p p colision "<<endl;
    cout<<"--> random seed      : "<<rndSeed<<endl;
    cout<<"--> proton energy    : "<<protonE<<" GeV/c"<<endl;
    cout<<"--> number of events : "<<nevents<<endl;
    cout<<"--> Out root file    : "<<outRootFileName<<endl;
    cout<<endl;

    return 0;
}
