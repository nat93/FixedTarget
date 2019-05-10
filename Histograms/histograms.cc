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

using namespace std;

TRandom3* _rnd;

void function_1(TString inputFileName);
void function_2(TString inputFileName);
void function_3();
void function_4();

Double_t calculateGamma(Double_t m, Double_t p);
Double_t calculateBetta(Double_t m, Double_t p);
Bool_t ifNotDecayInTarget(Double_t Ltarget, Double_t  tau, Double_t m, Double_t p);
Double_t notDecayRate(Double_t l, Double_t N0, Double_t tau, Double_t m, Double_t p);
Double_t calculateFlightTime(Double_t l, Double_t m, Double_t p);
Double_t getAngleBetweenLatticeAndtrack(Double_t px, Double_t py,Double_t pz, Double_t nx, Double_t ny,Double_t nz);
Bool_t ifThetaEffIsOk(Double_t angleBetweenLatticeAndtrack, Double_t Rcrystal, Double_t p, Double_t m, Int_t intCrystalTypeCondition);
Double_t dThetaEff_Si293( Double_t p, Double_t m, Double_t R);
Double_t dThetaEff_Ge293( Double_t p, Double_t m, Double_t R);
Double_t dThetaEff_Ge80( Double_t p, Double_t m, Double_t R);
Double_t dThetaEff( Double_t p, Double_t m, Double_t R, Double_t u0, Double_t dup);
Double_t getTotEnergy( Double_t p, Double_t m);
Bool_t ifAcceptanceFactorN0IsOk( Double_t p, Double_t m, Double_t Rcrystal, Int_t intCrystalTypeCondition);
Double_t acceptanceFactorN0_Si293( Double_t p, Double_t m, Double_t R);
Double_t acceptanceFactorN0_Ge293( Double_t p, Double_t m, Double_t R);
Double_t acceptanceFactorN0_Ge80( Double_t p, Double_t m, Double_t R);
Double_t acceptanceFactorN0( Double_t p, Double_t m, Double_t R, Double_t Nstr, Double_t kTheta, Double_t dup);
Bool_t ifNotDecay(Double_t Lcrystal, Double_t tau, Double_t m, Double_t p);
Bool_t ifNotDechannel(Double_t p, Double_t m, Double_t Rcrystal, Double_t Lcrystal, Int_t intCrystalTypeCondition);
Double_t dechannelLength_Si293(Double_t p, Double_t m, Double_t R);
Double_t dechannelLength_Ge293(Double_t p, Double_t m, Double_t R);
Double_t dechannelLength_Ge80(Double_t p, Double_t m, Double_t R);
Double_t dechannelLength(Double_t p, Double_t m, Double_t R, Double_t kdech, Double_t R0, Double_t bdech, Double_t dup_dech);
Bool_t ifNotDecayInTargetCrystal(Double_t Ltarget, Double_t Lcrystal, Double_t  tau, Double_t m, Double_t p);

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        cout<<endl;
        cout<<"--> 1 -- function_1() : histograms"<<endl;
        cout<<"--> 2 -- function_2() : (1/dg)^2"<<endl;
        cout<<"--> 3 -- function_3() : common (1/dg)^2"<<endl;
        cout<<"--> 4 -- function_4() : common function_1()"<<endl;
        cout<<endl;
        return -1;
    }

    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    _rnd = new TRandom3(seed);

    switch ( atoi(argv[1]) )
    {
    case 1:
      function_1(argv[2]);
      break;
    case 2:
      function_2(argv[2]);
      break;
    case 3:
      function_3();
      break;
    case 4:
      function_4();
      break;
    default:
      cout<<"--> Nothing to do =)"<<endl;
      break;
    }

    return 0;
}

void function_1(TString inputFileName)
{
    cout<<endl<<"--> function_1() <--"<<endl;
    cout<<endl;
    TString output_file_name    = inputFileName+"_output_function_1.root";

    // Common
    const int arraySize = 20000;
    Float_t _index[arraySize];
    Float_t _status[arraySize];
    Bool_t _IsFinal[arraySize];
    Float_t _ID[arraySize];
    Float_t _firstDau[arraySize];
    Float_t _lastDau[arraySize];
    Float_t _Mindex[arraySize];
    Float_t _xprod[arraySize];
    Float_t _yprod[arraySize];
    Float_t _zprod[arraySize];
    Float_t _xdecay[arraySize];
    Float_t _ydecay[arraySize];
    Float_t _zdecay[arraySize];
    Float_t _Px[arraySize];
    Float_t _Py[arraySize];
    Float_t _Pz[arraySize];
    Float_t _E[arraySize];
    Float_t _mass[arraySize];
    Int_t _Nparticle;
    Int_t _Nproj;
    Int_t _Ntarg;
    Int_t _Ncoll;

    TChain *fChain1 = new TChain("Tree");
    fChain1->Add(inputFileName.Data());
    cout<<"--> inputFileName = "<<inputFileName<<endl;

    fChain1->SetBranchAddress("Nparticle",    &_Nparticle);
    fChain1->SetBranchAddress("Nproj",        &_Nproj);
    fChain1->SetBranchAddress("Ntarg",        &_Ntarg);
    fChain1->SetBranchAddress("Ncoll",        &_Ncoll);
    fChain1->SetBranchAddress("index",        _index);
    fChain1->SetBranchAddress("status",       _status);
    fChain1->SetBranchAddress("IsFinal",      _IsFinal);
    fChain1->SetBranchAddress("Id",           _ID);
    fChain1->SetBranchAddress("firstDau",     _firstDau);
    fChain1->SetBranchAddress("lastDau",      _lastDau);
    fChain1->SetBranchAddress("Mindex",       _Mindex);
    fChain1->SetBranchAddress("xprod",        _xprod);
    fChain1->SetBranchAddress("yprod",        _yprod);
    fChain1->SetBranchAddress("zprod",        _zprod);
    fChain1->SetBranchAddress("xdecay",       _xdecay);
    fChain1->SetBranchAddress("ydecay",       _ydecay);
    fChain1->SetBranchAddress("zdecay",       _zdecay);
    fChain1->SetBranchAddress("Px",           _Px);
    fChain1->SetBranchAddress("Py",           _Py);
    fChain1->SetBranchAddress("Pz",           _Pz);
    fChain1->SetBranchAddress("E",            _E);
    fChain1->SetBranchAddress("mass",         _mass);

    Double_t nEntries = fChain1->GetEntries();
    cout<<"--> nEntries: "<<nEntries<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//
    TH1D* h_1 = new TH1D("h_1","L_{c} E [GeV/c^{2}] all",3000,0.0,300.0);
    TH1D* h_2 = new TH1D("h_2","L_{c} E [GeV/c^{2}, after target+crystal mm] all",3000,0.0,300.0);
    TH2D* h_3 = new TH2D("h_3","L_{c} E vs #Theta",4000,0,0.4,3000,0.0,300.0);
    TH2D* h_4 = new TH2D("h_4","L_{c} E vs #Theta (channeled)",4000,0,0.4,3000,0.0,300.0);
    TH2D* h_5 = new TH2D("h_5","L_{c} E vs #Theta (channeled & non-channeled)",4000,0,0.4,3000,0.0,300.0);
    TH1D* h_6 = new TH1D("h_6","L_{c} P [GeV/c] background",3000,0.0,300.0);
    TH1D* h_7 = new TH1D("h_7","L_{c} P [GeV/c] signal",3000,0.0,300.0);

    TTree* tree = new TTree("Tree", "Tree with Lc+ Signal and Background");
    Double_t Px, Py, Pz, E, mass;
    Int_t CHBKG;
    tree->Branch("Px",      &Px,    "Px/D");
    tree->Branch("Py",      &Py,    "Py/D");
    tree->Branch("Pz",      &Pz,    "Pz/D");
    tree->Branch("E",       &E,     "E/D");
    tree->Branch("mass",    &mass,  "mass/D");
    tree->Branch("CHBKG",   &CHBKG, "CHBKG/I");
    //--------------------------------------------------------------------------//

    const Double_t tau_LambdaC = 200.0*1.0E-6; //ns
    const Int_t intCrystalTypeCondition_Si293 = 0;
    const Double_t _Ltarget   = 3.0;        // mm
    const Double_t _Rcrystal  = 1.0e3;      // mm
    const Double_t _Lcrystal  = 15.0;       // mm
    const Double_t deltaTheta = 100.0e-6;   // rad

    Bool_t channeled    = kFALSE;
    Bool_t background   = kFALSE;

    for(Int_t eventID = 0; eventID < nEntries; eventID++)
    {
        fChain1->GetEntry(eventID);

        if(eventID%100 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)eventID/nEntries);
            fflush(stdout);
        }

        for(Int_t i = 0; i < _Nparticle; i++)
        {
            if(_ID[i] == 4122) // Lambda_c = 4122
            {
                channeled   = kFALSE;
                background  = kFALSE;

                h_1->Fill(_E[i]);
                h_3->Fill(TMath::Abs(TMath::ATan(_Px[i]/_Pz[i])),_E[i]);

                Double_t angleBetweenLatticeAndtrack = getAngleBetweenLatticeAndtrack( _Px[i], _Py[i], _Pz[i], 1, 0, 0);
                Double_t pmag = TMath::Sqrt(_Px[i]*_Px[i] + _Py[i]*_Py[i] + _Pz[i]*_Pz[i]);


//                if(ifNotDecayInTarget(_Ltarget, tau_LambdaC, _mass[i], pmag))
//                if(ifNotDecay(_Lcrystal, tau_LambdaC, _mass[i], pmag))
                if(ifNotDecayInTargetCrystal(_Ltarget, _Lcrystal, tau_LambdaC, _mass[i], pmag))
                {
                    h_2->Fill(_E[i]);

                    Px      = _Px[i];
                    Py      = _Py[i];
                    Pz      = _Pz[i];
                    E       = _E[i];
                    mass    = _mass[i];
                    CHBKG   = 0;

                    //----------------------------------------------------------//
                    // Background
                    if( TMath::ATan(_Px[i]/_Pz[i]) < (_Lcrystal/_Rcrystal+deltaTheta) && TMath::ATan(_Px[i]/_Pz[i]) > (_Lcrystal/_Rcrystal-deltaTheta) )
                    {
                        background = kTRUE;
                    }
                    //----------------------------------------------------------//
                    // Signal
                    if(ifThetaEffIsOk(angleBetweenLatticeAndtrack, _Rcrystal, pmag, _mass[i], intCrystalTypeCondition_Si293))
                    {
                        if(ifAcceptanceFactorN0IsOk( pmag, _mass[i], _Rcrystal, intCrystalTypeCondition_Si293))
                        {
                            if(ifNotDechannel(pmag, _mass[i], _Rcrystal, _Lcrystal, intCrystalTypeCondition_Si293))
                            {
                                h_4->Fill(TMath::Abs(TMath::ATan(_Px[i]/_Pz[i])),_E[i]);
                                channeled = kTRUE;
                            }
                        }
                    }
                    //----------------------------------------------------------//
                }

                if(channeled)
                {
                    h_5->Fill(TMath::Abs(TMath::ATan(_Px[i]/_Pz[i])) + _Lcrystal/_Rcrystal,_E[i]);
                    h_7->Fill(pmag);
                    CHBKG = 1;
                }
                else
                {
                    h_5->Fill(TMath::Abs(TMath::ATan(_Px[i]/_Pz[i])),_E[i]);
                }

                if(background)
                {
                    h_6->Fill(pmag);
                    CHBKG = 2;
                }

                if(channeled || background)
                {
                    tree->Fill();
                }
            }
        }
    }
    cout<<endl;
    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name,"recreate");

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_7->Write();

    tree->Write();

    file->Write();
    file->Close();
    //--------------------------------------------------------------------------//
}

void function_2(TString inputFileName)
{
    cout<<endl<<"--> function_2() <--"<<endl;
    cout<<endl<<"--> inputFileName = "<<inputFileName<<endl;
    TString output_file_name    = inputFileName+"_output_function_2.root";

    // Common
    const int arraySize = 20000;
    Float_t _index[arraySize];
    Float_t _status[arraySize];
    Bool_t _IsFinal[arraySize];
    Float_t _ID[arraySize];
    Float_t _firstDau[arraySize];
    Float_t _lastDau[arraySize];
    Float_t _Mindex[arraySize];
    Float_t _xprod[arraySize];
    Float_t _yprod[arraySize];
    Float_t _zprod[arraySize];
    Float_t _xdecay[arraySize];
    Float_t _ydecay[arraySize];
    Float_t _zdecay[arraySize];
    Float_t _Px[arraySize];
    Float_t _Py[arraySize];
    Float_t _Pz[arraySize];
    Float_t _E[arraySize];
    Float_t _mass[arraySize];
    Int_t _Nparticle;
    Int_t _Nproj;
    Int_t _Ntarg;
    Int_t _Ncoll;

    TChain *fChain1 = new TChain("Tree");
    fChain1->Add(inputFileName.Data());

    fChain1->SetBranchAddress("Nparticle",    &_Nparticle);
    fChain1->SetBranchAddress("Nproj",        &_Nproj);
    fChain1->SetBranchAddress("Ntarg",        &_Ntarg);
    fChain1->SetBranchAddress("Ncoll",        &_Ncoll);
    fChain1->SetBranchAddress("index",        _index);
    fChain1->SetBranchAddress("status",       _status);
    fChain1->SetBranchAddress("IsFinal",      _IsFinal);
    fChain1->SetBranchAddress("Id",           _ID);
    fChain1->SetBranchAddress("firstDau",     _firstDau);
    fChain1->SetBranchAddress("lastDau",      _lastDau);
    fChain1->SetBranchAddress("Mindex",       _Mindex);
    fChain1->SetBranchAddress("xprod",        _xprod);
    fChain1->SetBranchAddress("yprod",        _yprod);
    fChain1->SetBranchAddress("zprod",        _zprod);
    fChain1->SetBranchAddress("xdecay",       _xdecay);
    fChain1->SetBranchAddress("ydecay",       _ydecay);
    fChain1->SetBranchAddress("zdecay",       _zdecay);
    fChain1->SetBranchAddress("Px",           _Px);
    fChain1->SetBranchAddress("Py",           _Py);
    fChain1->SetBranchAddress("Pz",           _Pz);
    fChain1->SetBranchAddress("E",            _E);
    fChain1->SetBranchAddress("mass",         _mass);

    Double_t nEntries = fChain1->GetEntries();
    cout<<"--> nEntries: "<<nEntries<<endl;

    TGraph *polarization = new TGraph("LcPolarization.dat");
    polarization->SetName("polarization");
    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//
    TH2D* h_1 = new TH2D("h_1","(1/#Deltag(L,R))^{2}",30,0.0,30.0,30,0.0,3000.0);
    TH2D* h_2 = new TH2D("h_2","(1/#Deltag(#Theta,R))^{2}",60,0.0,60.0,30,0.0,3000.0);
    TH1D* h_3 = new TH1D("h_3","p_{T} L_{c}^{+} CH",300000,0.0,300.0);
    TH1D* h_4 = new TH1D("h_4","polarization L_{c}^{+} CH",20000,-10.0,10.0);
    //--------------------------------------------------------------------------//

    const Double_t tau_LambdaC = 200.0*1.0E-6;  //ns
    const Int_t intCrystalTypeCondition_Si293 = 0;
    const Double_t _Ltarget   = 3.0;            // mm

    for(Int_t eventID = 0; eventID < nEntries; eventID++)
    {
        fChain1->GetEntry(eventID);

        if(eventID%100 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)eventID/nEntries);
            fflush(stdout);
        }

        for(Double_t _Rcrystal = 100.0; _Rcrystal < 3000.0; _Rcrystal += 100.0)
        {
            for(Double_t _Lcrystal = 1.0; _Lcrystal < 30.0; _Lcrystal += 1.0)
            {
                for(Int_t i = 0; i < _Nparticle; i++)
                {
                    if(_ID[i] == 4122) // Lambda_c = 4122
                    {
                        Double_t angleBetweenLatticeAndtrack = getAngleBetweenLatticeAndtrack( _Px[i], _Py[i], _Pz[i], 1, 0, 0);
                        Double_t pmag = TMath::Sqrt(_Px[i]*_Px[i] + _Py[i]*_Py[i] + _Pz[i]*_Pz[i]);

//                        if(ifNotDecayInTarget(_Ltarget, tau_LambdaC, _mass[i], pmag))
//                        if(ifNotDecay(_Lcrystal, tau_LambdaC, _mass[i], pmag))
                        if(ifNotDecayInTargetCrystal(_Ltarget, _Lcrystal, tau_LambdaC, _mass[i], pmag))
                        {
                            if(ifThetaEffIsOk(angleBetweenLatticeAndtrack, _Rcrystal, pmag, _mass[i], intCrystalTypeCondition_Si293))
                            {
                                if(ifAcceptanceFactorN0IsOk( pmag, _mass[i], _Rcrystal, intCrystalTypeCondition_Si293))
                                {
                                    if(ifNotDechannel(pmag, _mass[i], _Rcrystal, _Lcrystal, intCrystalTypeCondition_Si293))
                                    {
                                        Double_t pT = TMath::Sqrt(_Px[i]*_Px[i] + _Py[i]*_Py[i]);
                                        Double_t xi = TMath::Abs(polarization->Eval(pT));
                                        Double_t theta = _Lcrystal/_Rcrystal;
                                        Double_t ddg = TMath::Power(calculateGamma(_mass[i],pmag)*theta*xi,2);

                                        h_1->Fill(_Lcrystal,_Rcrystal,ddg);
                                        h_2->Fill(theta*1000,_Rcrystal,ddg);
                                        h_3->Fill(pT);
                                        h_4->Fill(xi);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout<<endl;
    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name,"recreate");

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    polarization->Write();

    file->Write();
    file->Close();
    //--------------------------------------------------------------------------//
}

void function_3()
{
    cout<<endl<<"--> function_3() <--"<<endl;
    const Int_t nFiles = 20;
    TFile *_file[nFiles];
    TH2D* h_1 = new TH2D("h_1","(1/#Deltag(L,R))^{2}",30,0.0,30.0,30,0.0,3000.0);
    TH2D* h_2 = new TH2D("h_2","(1/#Deltag(#Theta,R))^{2}",60,0.0,60.0,30,0.0,3000.0);
    TH1D* h_3 = new TH1D("h_3","p_{T} L_{c}^{+} CH",300000,0.0,300.0);
    TH1D* h_4 = new TH1D("h_4","polarization L_{c}^{+} CH",20000,-10.0,10.0);

    TString name;
    for(Int_t i = 0; i < nFiles; i++)
    {
    name = "../../home2/SPS_Sim_Data/hardccbar_";
        name += i+1;
        name += ".root_output_function_2.root";
        _file[i] = TFile::Open(name.Data());
        if(_file[i]->IsOpen())
            cout<<"--> Filename: "<<name<<endl;

        TH2D* hh_1 = (TH2D*)_file[i]->Get("h_1");
        TH2D* hh_2 = (TH2D*)_file[i]->Get("h_2");
        TH1D* hh_3 = (TH1D*)_file[i]->Get("h_3");
        TH1D* hh_4 = (TH1D*)_file[i]->Get("h_4");

        h_1->Add(hh_1);
        h_2->Add(hh_2);
        h_3->Add(hh_3);
        h_4->Add(hh_4);
    }

    cout<<endl;
    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    TString output_file_name = "common_histo_output_function_3.root";
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name,"recreate");

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();

    file->Write();
    file->Close();
    //--------------------------------------------------------------------------//

}

void function_4()
{
    cout<<endl<<"--> function_4() <--"<<endl;
    const Int_t nFiles = 20;
    TFile *_file[nFiles];
    TH1D* h_1 = new TH1D("h_1","L_{c} E [GeV/c^{2}] all",3000,0.0,300.0);
    TH1D* h_2 = new TH1D("h_2","L_{c} E [GeV/c^{2}, after target+crystal mm] all",3000,0.0,300.0);
    TH2D* h_3 = new TH2D("h_3","L_{c} E vs #Theta",4000,0,0.4,3000,0.0,300.0);
    TH2D* h_4 = new TH2D("h_4","L_{c} E vs #Theta (channeled)",4000,0,0.4,3000,0.0,300.0);
    TH2D* h_5 = new TH2D("h_5","L_{c} E vs #Theta (channeled & non-channeled)",4000,0,0.4,3000,0.0,300.0);
    TH1D* h_6 = new TH1D("h_6","L_{c} P [GeV/c] background",3000,0.0,300.0);
    TH1D* h_7 = new TH1D("h_7","L_{c} P [GeV/c] signal",3000,0.0,300.0);

    TString name;
    for(Int_t i = 0; i < nFiles; i++)
    {
        name = "../../home2/SPS_Sim_Data/hardccbar_";
        name += i+1;
        name += ".root_output_function_1.root";
        _file[i] = TFile::Open(name.Data());
        if(_file[i]->IsOpen())
            cout<<"--> Filename: "<<name<<endl;

        TH1D* hh_1 = (TH1D*)_file[i]->Get("h_1");
        TH1D* hh_2 = (TH1D*)_file[i]->Get("h_2");
        TH2D* hh_3 = (TH2D*)_file[i]->Get("h_3");
        TH2D* hh_4 = (TH2D*)_file[i]->Get("h_4");
        TH2D* hh_5 = (TH2D*)_file[i]->Get("h_5");
        TH1D* hh_6 = (TH1D*)_file[i]->Get("h_6");
        TH1D* hh_7 = (TH1D*)_file[i]->Get("h_7");

        h_1->Add(hh_1);
        h_2->Add(hh_2);
        h_3->Add(hh_3);
        h_4->Add(hh_4);
        h_5->Add(hh_5);
        h_6->Add(hh_6);
        h_7->Add(hh_7);
    }

    cout<<endl;
    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    TString output_file_name = "common_histo_output_function_4.root";
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name,"recreate");

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_7->Write();

    file->Write();
    file->Close();
    //--------------------------------------------------------------------------//

}

Double_t calculateGamma(Double_t m, Double_t p)
{
    if(p>=0.0)
    {
        if(m>0.0)
        {
            return TMath::Sqrt(1.0 + p*p/m/m);
        }
    }
    return -999.0;
}

Double_t calculateBetta(Double_t m, Double_t p)
{
    if(p>=0.0)
    {
        if(m>0.0)
        {
            return 1.0/(TMath::Sqrt(1.0 + m*m/p/p));
        }
        else if( m == 0.0)
        {
            return 1.0;
        }
        else
        {
            return -999.0;
        }
    }
    return -999.0;
}

Bool_t ifNotDecayInTargetCrystal(Double_t Ltarget, Double_t Lcrystal, Double_t  tau, Double_t m, Double_t p)
{
    Double_t N0 = 1.0;
    Double_t lInTarget = Ltarget - _rnd->Uniform(0.0,Ltarget);
    Double_t probNotDec = notDecayRate(lInTarget+Lcrystal, N0, tau, m, p);
    if(probNotDec > _rnd->Uniform(0.0,1.0))
        return true;
    return false;
}

Bool_t ifNotDecayInTarget(Double_t Ltarget, Double_t  tau, Double_t m, Double_t p)
{
    Double_t N0 = 1.0;
    Double_t lInTarget = Ltarget - _rnd->Uniform(0.0,Ltarget);
    Double_t probNotDec = notDecayRate(lInTarget, N0, tau, m, p);
    if(probNotDec > _rnd->Uniform(0.0,1.0))
        return true;
    return false;
}

Double_t notDecayRate(Double_t l, Double_t N0, Double_t tau, Double_t m, Double_t p)
{
    Double_t gamma;
    Double_t t;
    if(l >= 0.0)
    {
        if(N0 >= 0.0)
        {
            if(tau >= 0.0)
            {
                if(m >= 0.0)
                {
                    if(p >= 0.0)
                    {
                        gamma = calculateGamma(m, p);
                        t = calculateFlightTime(l, m, p);
                        if(gamma >= 1.0)
                        {
                            if(t >= 0.0)
                            {
                                return N0*TMath::Exp(-t/tau/gamma);
                            }
                        }
                    }
                }
            }
        }
    }
    return -999.0;
}

Double_t calculateFlightTime(Double_t l, Double_t m, Double_t p)
{
    const Double_t c = 299.7924580000;         //mm/ns
    Double_t betta;
    if(l>0.0)
    {
        if(m>=0.0)
        {
            if(p>0.0)
            {
                betta = calculateBetta(m, p);
                if(betta>0.0)
                {
                    return l/c/betta;
                }
            }
        }
    }
    return -999.0;
}

Double_t getAngleBetweenLatticeAndtrack(Double_t px, Double_t py,Double_t pz, Double_t nx, Double_t ny,Double_t nz)
{
    TVector3 linDir(px,py,pz);
    TVector3 planeN(nx,ny,nz);
    Double_t AA = linDir.Dot(planeN);
    Double_t BB = planeN.Mag()*linDir.Mag();
    if(BB != 0.0)
    {
        return TMath::ASin(AA/BB);
    }
    else
    {
        cout<<"ERROR ---> znam == 0.0"<<endl;
        assert(0);
    }
    return -999.0;
}

Bool_t ifThetaEffIsOk(Double_t angleBetweenLatticeAndtrack, Double_t Rcrystal, Double_t p, Double_t m, Int_t intCrystalTypeCondition)
{
    Double_t thetaEff;
    Double_t thetaEffAbs;
    Double_t angleBetweenLatticeAndtrackAbs;

    const Int_t intCrystalTypeCondition_Si293   = 0;
    const Int_t intCrystalTypeCondition_Ge293   = 1;
    const Int_t intCrystalTypeCondition_Ge80    = 2;

    if(intCrystalTypeCondition == intCrystalTypeCondition_Si293)
    {
        thetaEff = dThetaEff_Si293( p, m, Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge293)
    {
        thetaEff = dThetaEff_Ge293( p, m, Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge80)
    {
        thetaEff = dThetaEff_Ge80( p, m, Rcrystal);
    }
    else
    {
        assert(0);
    }
    thetaEffAbs = TMath::Abs(thetaEff);
    angleBetweenLatticeAndtrackAbs = TMath::Abs(angleBetweenLatticeAndtrack);
    if(angleBetweenLatticeAndtrackAbs <= thetaEffAbs)
        return true;
    return false;
}

Double_t dThetaEff_Si293( Double_t p, Double_t m, Double_t R)
{
    const Double_t u0_Si293  = 22.80*1.0E-9; //GeV
    const Double_t dup_Si293 = 0.525;        //GeV/mm

    Double_t u0  = u0_Si293;
    Double_t dup = dup_Si293;
    return dThetaEff( p, m, R, u0, dup);
}

Double_t dThetaEff_Ge293( Double_t p, Double_t m, Double_t R)
{
    const Double_t u0_Ge293  = 40.00*1.0E-9; //GeV
    const Double_t dup_Ge293 = 0.922;        //GeV/mm

    Double_t u0  = u0_Ge293;
    Double_t dup = dup_Ge293;
    return dThetaEff( p, m, R, u0, dup);
}

Double_t dThetaEff_Ge80( Double_t p, Double_t m, Double_t R)
{
    const Double_t u0_Ge80   = 45.20*1.0E-9; //GeV
    const Double_t dup_Ge80  = 1.200;        //GeV/mm

    Double_t u0  = u0_Ge80;
    Double_t dup = dup_Ge80;
    return dThetaEff( p, m, R, u0, dup);
}

Double_t dThetaEff( Double_t p, Double_t m, Double_t R, Double_t u0, Double_t dup)
{
    Double_t E;
    Double_t val;
    if(p>=0.0)
    {
        if(m>0.0)
        {
            E = getTotEnergy(p,m);
            if(R>0.0)
            {
                if(u0>0.0)
                {
                    if(dup>0.0)
                    {
                        val = TMath::Sqrt(2.0*u0/E)*(1-E/R/dup);
                        if(val<0.0)
                            return 0.0;
                        return val;
                    }
                }
            }
        }
    }
    return -999.0;
}

Double_t getTotEnergy( Double_t p, Double_t m)
{
    if(p>=0.0)
    {
        if(m>=0.0)
        {
            return TMath::Sqrt(p*p + m*m);
        }
    }
    return -999.0;
}

Bool_t ifAcceptanceFactorN0IsOk( Double_t p, Double_t m, Double_t Rcrystal, Int_t intCrystalTypeCondition)
{
    const Int_t intCrystalTypeCondition_Si293   = 0;
    const Int_t intCrystalTypeCondition_Ge293   = 1;
    const Int_t intCrystalTypeCondition_Ge80    = 2;

    Double_t acceptanceFactor;
    if(intCrystalTypeCondition == intCrystalTypeCondition_Si293)
    {
        acceptanceFactor = acceptanceFactorN0_Si293( p, m, Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge293)
    {
        acceptanceFactor = acceptanceFactorN0_Ge293( p, m, Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge80)
    {
        acceptanceFactor = acceptanceFactorN0_Ge80( p, m, Rcrystal);
    }
    else
    {
        assert(0);
    }
    if(acceptanceFactor>_rnd->Uniform(0.0,1.0))
        return true;
    return false;
}

Double_t acceptanceFactorN0_Si293( Double_t p, Double_t m, Double_t R)
{
    const Double_t Nstr_Si293 = 0.775;
    const Double_t kTheta_Si293 = 0.405;
    const Double_t dup_Si293 = 0.525;        //GeV/mm

    Double_t Nstr = Nstr_Si293;
    Double_t kTheta = kTheta_Si293;
    Double_t dup = dup_Si293;
    return acceptanceFactorN0( p, m, R, Nstr, kTheta, dup);
}

Double_t acceptanceFactorN0_Ge293( Double_t p, Double_t m, Double_t R)
{
    const Double_t Nstr_Ge293 = 0.785;
    const Double_t kTheta_Ge293 = 0.396;
    const Double_t dup_Ge293 = 0.922;        //GeV/mm

    Double_t Nstr = Nstr_Ge293;
    Double_t kTheta = kTheta_Ge293;
    Double_t dup = dup_Ge293;
    return acceptanceFactorN0( p, m, R, Nstr, kTheta, dup);
}

Double_t acceptanceFactorN0_Ge80( Double_t p, Double_t m, Double_t R)
{
    const Double_t Nstr_Ge80  = 0.812;
    const Double_t kTheta_Ge80  = 0.353;
    const Double_t dup_Ge80  = 1.200;        //GeV/mm

    Double_t Nstr = Nstr_Ge80;
    Double_t kTheta = kTheta_Ge80;
    Double_t dup = dup_Ge80;
    return acceptanceFactorN0( p, m, R, Nstr, kTheta, dup);
}

Double_t acceptanceFactorN0( Double_t p, Double_t m, Double_t R, Double_t Nstr, Double_t kTheta, Double_t dup)
{
    Double_t E;
    if(p>=0.0)
    {
        if(m>0.0)
        {
            E = getTotEnergy(p,m);
            if(R>0.0)
            {
                if(Nstr>0.0)
                {
                    if( kTheta>0.0)
                    {
                        if(dup>0.0)
                        {
                            return Nstr/(1.0 + E*E/R/R/dup/dup/kTheta/kTheta);
                        }
                    }
                }
            }
        }
    }
    return -999.0;
}

Bool_t ifNotDecay(Double_t Lcrystal, Double_t tau, Double_t m, Double_t p)
{
    Double_t N0 = 1.0;
    Double_t probNotDec = notDecayRate(Lcrystal, N0, tau, m, p);
    if(probNotDec>_rnd->Uniform(0.0,1.0))
        return true;
    return false;
}

Bool_t ifNotDechannel(Double_t p, Double_t m, Double_t Rcrystal, Double_t Lcrystal, Int_t intCrystalTypeCondition)
{
    const Int_t intCrystalTypeCondition_Si293   = 0;
    const Int_t intCrystalTypeCondition_Ge293   = 1;
    const Int_t intCrystalTypeCondition_Ge80    = 2;

    Double_t Ldechannel;
    Double_t chanelProb;
    if(intCrystalTypeCondition == intCrystalTypeCondition_Si293)
    {
        Ldechannel = dechannelLength_Si293(p,m,Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge293)
    {
        Ldechannel = dechannelLength_Ge293(p,m,Rcrystal);
    }
    else if(intCrystalTypeCondition == intCrystalTypeCondition_Ge80)
    {
        Ldechannel = dechannelLength_Ge80(p,m,Rcrystal);
    }
    else
    {
        assert(0);
    }
    if(Ldechannel>0.0)
    {
        chanelProb = TMath::Exp(-TMath::Sqrt(Lcrystal/Ldechannel));
        if(chanelProb>_rnd->Uniform(0.0,1.0))
            return true;
    }
    else
    {
        return false;
    }
    return false;
}

Double_t dechannelLength_Si293(Double_t p, Double_t m, Double_t R)
{
    const Double_t kdech_Si293      = 0.0146;
    const Double_t R0_Si293         = 3700.0;   //mm
    const Double_t bdech_Si293      = 0.1;
    const Double_t dup_dech_Si293   = 0.097;    //GeV/mm

    Double_t kdech = kdech_Si293;
    Double_t R0 = R0_Si293;
    Double_t bdech = bdech_Si293;
    Double_t dup_dech = dup_dech_Si293;

    return dechannelLength(p, m, R, kdech, R0, bdech, dup_dech);
}

Double_t dechannelLength_Ge293(Double_t p, Double_t m, Double_t R)
{
    const Double_t kdech_Ge293      = 0.0130;
    const Double_t R0_Ge293         = 3300.0;   //mm
    const Double_t bdech_Ge293      = 0.2;
    const Double_t dup_dech_Ge293   = 0.137;    //GeV/mm

    Double_t kdech = kdech_Ge293;
    Double_t R0 = R0_Ge293;
    Double_t bdech = bdech_Ge293;
    Double_t dup_dech = dup_dech_Ge293;

    return dechannelLength(p, m, R, kdech, R0, bdech, dup_dech);
}

Double_t dechannelLength_Ge80(Double_t p, Double_t m, Double_t R)
{
    const Double_t kdech_Ge80       = 0.0226;
    const Double_t R0_Ge80          = 2900.0;   //mm
    const Double_t bdech_Ge80       = 0.2;
    const Double_t dup_dech_Ge80    = 0.146;    //GeV/mm

    Double_t kdech = kdech_Ge80;
    Double_t R0 = R0_Ge80;
    Double_t bdech = bdech_Ge80;
    Double_t dup_dech = dup_dech_Ge80;

    return dechannelLength(p, m, R, kdech, R0, bdech, dup_dech);
}

Double_t dechannelLength(Double_t p, Double_t m, Double_t R, Double_t kdech, Double_t R0, Double_t bdech, Double_t dup_dech)
{
    Double_t Lmax;
    Double_t Emax;
    Double_t E;
    if(p>=0.0)
    {
        if(m>0.0)
        {
            if(R>0.0)
            {
                if(kdech>0.0)
                {
                    if(R0>0.0)
                    {
                        if(bdech>0.0)
                        {
                            if(dup_dech>0.0)
                            {
                                Lmax = kdech*R*TMath::Power(R0/R,bdech);
                                Emax = R*dup_dech;
                                E = getTotEnergy(p,m);
                                return Lmax*E/Emax*TMath::Exp(1.0-E/Emax);
                            }
                        }
                    }
                }
            }
        }
    }
    return -999.0;
}
