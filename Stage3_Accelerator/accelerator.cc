//my
#include "./src/MagClass.h"
#include "src/Constants.h"

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
#include "TH2D.h"
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

void function_1(TString output_file_name);
void function_2(TString input_filedir_name, TString output_file_name);
void function_3(TString output_file_name);

Bool_t isChanneled(Double_t particleAngleOut);
void passMagnets(Double_t* coord_x0, Double_t* coord_x, Double_t p, Double_t Charge_C, TGraph* gr_x, TGraph* gr_y, Bool_t print);

Int_t main(int argc, char* argv[])
{
    if(argc == 4)
    {
        switch ( atoi(argv[3]) )
        {
        case 1:
          function_1(argv[2]);
          break;
        case 2:
          function_2(argv[1],argv[2]);
          break;
        case 3:
          function_3(argv[2]);
          break;
        default:
          cout<<"--> Nothing to do =)"<<endl;
          break;
        }
    }
    else
    {
        cout<<endl;
        cout<<"--> 1 -- function_1() : single proton trajectory (-200 urad deflection)"<<endl;
        cout<<"--> 2 -- function_2() : Lc products trajectories"<<endl;
        cout<<"--> 3 -- function_3() : for random parameters"<<endl;
        cout<<endl;
        cout<<"--> ERROR:: Wrong imput parameters number:"<<endl<<
              "--> [0] -- script name"<<endl<<
              "--> [1] -- input filename"<<endl<<
              "--> [2] -- output filename"<<endl<<
              "--> [3] -- function id"<<endl;
        return -1;

    }

    return 0;
}

void function_1(TString output_file_name)
{
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//

    Double_t* coord_x0          = new Double_t[mtrx_size];
    Double_t* coord_x           = new Double_t[mtrx_size];

    Double_t p = 270.0;                         // [GeV/c]
    Double_t Charge_C   = 1;

    /*x*/     coord_x0[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
    /*x'*/    coord_x0[1] = Constants::_beamAngleInitialAtCryPosition+Constants::_crystalAngle;// [rad]
    /*y*/     coord_x0[2] = 0.0;                // [m]
    /*y'*/    coord_x0[3] = 0.0;                // [rad]
    /*l*/     coord_x0[4] = 0.0;
    /*dp*/    coord_x0[5] = 0.0;

    TGraph* gr_x = new TGraph();
    gr_x->SetName("gr_x");
    TGraph* gr_y = new TGraph();
    gr_y->SetName("gr_y");

    passMagnets(coord_x0,coord_x,p,Charge_C,gr_x,gr_y,true);

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");
    gr_x->Write();
    gr_y->Write();
    _file->Close();
    cout<<"--> Output file: "<<output_file_name<<endl;
}

void function_2(TString input_filedir_name, TString output_file_name)
{
    //-----------------------------------------//
    // Read file
    //-----------------------------------------//

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
    Bool_t _isGenerated;
    Float_t _Px_Lc_0;
    Float_t _Py_Lc_0;
    Float_t _Pz_Lc_0;
    Bool_t _isTarget;
    Float_t _Px_Lc_1;
    Float_t _Py_Lc_1;
    Float_t _Pz_Lc_1;
    Bool_t _isCrystal;
    Float_t _Px_Lc_2;
    Float_t _Py_Lc_2;
    Float_t _Pz_Lc_2;
    Float_t _E[arraySize];
    Float_t _mass[arraySize];
    Int_t _Nparticle;
    Int_t _Nproj;
    Int_t _Ntarg;
    Int_t _Ncoll;

    TChain* fChain1 = new TChain("Tree");
    TString input_file_name;
    for(Int_t i = 2; i <= 20; i++)
    {
        input_file_name = input_filedir_name;
        input_file_name += "crystaltarget_";
        input_file_name += i;
        input_file_name += ".root";
        fChain1->Add(input_file_name.Data());
    }

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
    fChain1->SetBranchAddress("isGenerated",  &_isGenerated);
    fChain1->SetBranchAddress("Px_Lc_0",      &_Px_Lc_0);
    fChain1->SetBranchAddress("Py_Lc_0",      &_Py_Lc_0);
    fChain1->SetBranchAddress("Pz_Lc_0",      &_Pz_Lc_0);
    fChain1->SetBranchAddress("isTarget",     &_isTarget);
    fChain1->SetBranchAddress("Px_Lc_1",      &_Px_Lc_1);
    fChain1->SetBranchAddress("Py_Lc_1",      &_Py_Lc_1);
    fChain1->SetBranchAddress("Pz_Lc_1",      &_Pz_Lc_1);
    fChain1->SetBranchAddress("isCrystal",    &_isCrystal);
    fChain1->SetBranchAddress("Px_Lc_2",      &_Px_Lc_2);
    fChain1->SetBranchAddress("Py_Lc_2",      &_Py_Lc_2);
    fChain1->SetBranchAddress("Pz_Lc_2",      &_Pz_Lc_2);
    fChain1->SetBranchAddress("E",            _E);
    fChain1->SetBranchAddress("mass",         _mass);

    cout<<"--> Input file: "<<input_filedir_name<<endl;
    Long64_t nEntries = fChain1->GetEntries();
    cout<<"--> nEntries: "<<nEntries<<endl;

    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//

    Double_t p_proton = -999, charge_proton = -999, theta_proton = -999;
    Double_t* coord_x0_proton          = new Double_t[mtrx_size];
    Double_t* coord_x_proton           = new Double_t[mtrx_size];

    Double_t p_kaon = -999, charge_kaon = -999, theta_kaon = -999;
    Double_t* coord_x0_kaon          = new Double_t[mtrx_size];
    Double_t* coord_x_kaon           = new Double_t[mtrx_size];

    Double_t p_pion = -999, charge_pion = -999, theta_pion = -999;
    Double_t* coord_x0_pion          = new Double_t[mtrx_size];
    Double_t* coord_x_pion           = new Double_t[mtrx_size];

    Int_t LambdaC_ind;//, proton_ind, kaon_ind, pion_ind;
    Bool_t isLambdacgenerated, isLambdacProtongenerated, isLambdacKaonMgenerated, isLambdacPionPgenerated;

    TH1D* h_1 = new TH1D("h_1","L_{c} momentum channeled",330,-30.0,300.0);
    TH1D* h_2 = new TH1D("h_2","L_{c} angle",4000000,-4.0,4.0);
    TH2D* h_3 = new TH2D("h_3","XY proton on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_4 = new TH2D("h_4","XY kaon on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_5 = new TH2D("h_5","XY pion on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_6 = new TH2D("h_6","XY proton on RP3",1000,-10,10,1000,-10,10);
    TH2D* h_7 = new TH2D("h_7","XY kaon on RP3",1000,-10,10,1000,-10,10);
    TH2D* h_8 = new TH2D("h_8","XY pion on RP3",1000,-10,10,1000,-10,10);
    TH1D* h_9 = new TH1D("h_9","proton momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_10 = new TH1D("h_10","kaon momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_11 = new TH1D("h_11","pion momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_12 = new TH1D("h_12","proton angle",4000000,-4.0,4.0);
    TH1D* h_13 = new TH1D("h_13","kaon angle",4000000,-4.0,4.0);
    TH1D* h_14 = new TH1D("h_14","pion angle",4000000,-4.0,4.0);

    Int_t nRuns = 0;

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");

    for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        fChain1->GetEntry(iEntry);

        if(iEntry%1000 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)iEntry/nEntries);
            fflush(stdout);
        }

        isLambdacgenerated          = false;
        isLambdacProtongenerated    = false;
        isLambdacKaonMgenerated     = false;
        isLambdacPionPgenerated     = false;

        if(_isGenerated && _isTarget && _isCrystal)
        {
            Double_t thetaLambdaCout = TMath::ATan(_Px_Lc_2/_Pz_Lc_2) - TMath::ATan(_Px_Lc_1/_Pz_Lc_1);

            h_2->Fill(thetaLambdaCout);

            if(isChanneled(thetaLambdaCout))
            {
                for(Int_t iParticle = 0; iParticle < _Nparticle; iParticle++)
                {
//                    cout<<_ID[iParticle]<<endl;
                    if(_ID[iParticle] == 4122) // lambda_c+
                    {
                        isLambdacgenerated = true;
                        LambdaC_ind = iParticle;

                        h_1->Fill(TMath::Sqrt(_Px[iParticle]*_Px[iParticle]+_Py[iParticle]*_Py[iParticle]+_Pz[iParticle]*_Pz[iParticle]));
                    }

                    if(_ID[iParticle] == 2212) // proton
                    {
                        if(_Mindex[iParticle] == LambdaC_ind)// && _IsFinal[iParticle])
                        {
                            isLambdacProtongenerated = true;
//                            proton_ind = iParticle;

                            p_proton = TMath::Sqrt(_Px[iParticle]*_Px[iParticle] + _Py[iParticle]*_Py[iParticle] + _Pz[iParticle]*_Pz[iParticle]);
                            theta_proton = TMath::ATan(_Px[iParticle]/_Pz[iParticle]);
                            charge_proton = 1.0;
                            /*x*/     coord_x0_proton[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                            /*x'*/    coord_x0_proton[1] =
                                    Constants::_beamAngleInitialAtCryPosition +
                                    Constants::_crystalAngle +
                                    TMath::ATan(_Px[iParticle]/_Pz[iParticle]);// [rad]
                            /*y*/     coord_x0_proton[2] = 0.0;                // [m]
                            /*y'*/    coord_x0_proton[3] = TMath::ATan(_Py[iParticle]/_Pz[iParticle]); // [rad]
                            /*l*/     coord_x0_proton[4] = 0.0;
                            /*dp*/    coord_x0_proton[5] = 0.0;
                        }
                    }

                    if(_ID[iParticle] == -321 || _ID[iParticle] == -323) // kaon- || kaon-*
                    {
                        if(_Mindex[iParticle] == LambdaC_ind)// && _IsFinal[iParticle])
                        {
                            isLambdacKaonMgenerated = true;
//                            kaon_ind = iParticle;

                            p_kaon = TMath::Sqrt(_Px[iParticle]*_Px[iParticle] + _Py[iParticle]*_Py[iParticle] + _Pz[iParticle]*_Pz[iParticle]);
                            theta_kaon = TMath::ATan(_Px[iParticle]/_Pz[iParticle]);
                            charge_kaon = -1.0;
                            /*x*/     coord_x0_kaon[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                            /*x'*/    coord_x0_kaon[1] =
                                    Constants::_beamAngleInitialAtCryPosition +
                                    Constants::_crystalAngle +
                                    TMath::ATan(_Px[iParticle]/_Pz[iParticle]);// [rad]
                            /*y*/     coord_x0_kaon[2] = 0.0;                // [m]
                            /*y'*/    coord_x0_kaon[3] = TMath::ATan(_Py[iParticle]/_Pz[iParticle]); // [rad]
                            /*l*/     coord_x0_kaon[4] = 0.0;
                            /*dp*/    coord_x0_kaon[5] = 0.0;
                        }
                    }

                    if(_ID[iParticle] == 211) // pion+
                    {
                        if(_Mindex[iParticle] == LambdaC_ind)// && _IsFinal[iParticle])
                        {
                            isLambdacPionPgenerated = true;
//                            pion_ind = iParticle;

                            p_pion = TMath::Sqrt(_Px[iParticle]*_Px[iParticle] + _Py[iParticle]*_Py[iParticle] + _Pz[iParticle]*_Pz[iParticle]);
                            theta_pion = TMath::ATan(_Px[iParticle]/_Pz[iParticle]);
                            charge_pion = 1.0;
                            /*x*/     coord_x0_pion[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                            /*x'*/    coord_x0_pion[1] =
                                    Constants::_beamAngleInitialAtCryPosition +
                                    Constants::_crystalAngle +
                                    TMath::ATan(_Px[iParticle]/_Pz[iParticle]);// [rad]
                            /*y*/     coord_x0_pion[2] = 0.0;                // [m]
                            /*y'*/    coord_x0_pion[3] = TMath::ATan(_Py[iParticle]/_Pz[iParticle]); // [rad]
                            /*l*/     coord_x0_pion[4] = 0.0;
                            /*dp*/    coord_x0_pion[5] = 0.0;
                        }
                    }
                }

                if(isLambdacgenerated && isLambdacProtongenerated && isLambdacKaonMgenerated && isLambdacPionPgenerated)
                {
                    h_9->Fill(p_proton);
                    h_10->Fill(p_kaon);
                    h_11->Fill(p_pion);
                    h_12->Fill(theta_proton);
                    h_13->Fill(theta_kaon);
                    h_14->Fill(theta_pion);

                    TGraph* gr_proton_x = new TGraph();
                    TGraph* gr_proton_y = new TGraph();
                    TString gr_proton_x_name = "gr_proton_x_";
                    TString gr_proton_y_name = "gr_proton_y_";
                    gr_proton_x_name += nRuns;
                    gr_proton_y_name += nRuns;
                    gr_proton_x->SetName(gr_proton_x_name.Data());
                    gr_proton_y->SetName(gr_proton_y_name.Data());

                    passMagnets(coord_x0_proton,coord_x_proton,p_proton,charge_proton,gr_proton_x,gr_proton_y,false);

                    TGraph* gr_kaon_x = new TGraph();
                    TGraph* gr_kaon_y = new TGraph();
                    TString gr_kaon_x_name = "gr_kaon_x_";
                    TString gr_kaon_y_name = "gr_kaon_y_";
                    gr_kaon_x_name += nRuns;
                    gr_kaon_y_name += nRuns;
                    gr_kaon_x->SetName(gr_kaon_x_name.Data());
                    gr_kaon_y->SetName(gr_kaon_y_name.Data());

                    passMagnets(coord_x0_kaon,coord_x_kaon,p_kaon,charge_kaon,gr_kaon_x,gr_kaon_y,false);

                    TGraph* gr_pion_x = new TGraph();
                    TGraph* gr_pion_y = new TGraph();
                    TString gr_pion_x_name = "gr_pion_x_";
                    TString gr_pion_y_name = "gr_pion_y_";
                    gr_pion_x_name += nRuns;
                    gr_pion_y_name += nRuns;
                    gr_pion_x->SetName(gr_pion_x_name.Data());
                    gr_pion_y->SetName(gr_pion_y_name.Data());

                    passMagnets(coord_x0_pion,coord_x_pion,p_pion,charge_pion,gr_pion_x,gr_pion_y,false);

                    nRuns++;

                    gr_proton_x->Write();
                    gr_proton_y->Write();
                    gr_kaon_x->Write();
                    gr_kaon_y->Write();
                    gr_pion_x->Write();
                    gr_pion_y->Write();

                    h_3->Fill(gr_proton_x->Eval(Constants::_xrph_51937_ua9_pos),gr_proton_y->Eval(Constants::_xrph_51937_ua9_pos));
                    h_4->Fill(gr_kaon_x->Eval(Constants::_xrph_51937_ua9_pos),gr_kaon_y->Eval(Constants::_xrph_51937_ua9_pos));
                    h_5->Fill(gr_pion_x->Eval(Constants::_xrph_51937_ua9_pos),gr_pion_y->Eval(Constants::_xrph_51937_ua9_pos));
                    h_6->Fill(gr_proton_x->Eval(Constants::_xrph_52202_ua9_pos),gr_proton_y->Eval(Constants::_xrph_52202_ua9_pos));
                    h_7->Fill(gr_kaon_x->Eval(Constants::_xrph_52202_ua9_pos),gr_kaon_y->Eval(Constants::_xrph_52202_ua9_pos));
                    h_8->Fill(gr_pion_x->Eval(Constants::_xrph_52202_ua9_pos),gr_pion_y->Eval(Constants::_xrph_52202_ua9_pos));

                    gr_proton_x->Delete();
                    gr_proton_y->Delete();
                    gr_kaon_x->Delete();
                    gr_kaon_y->Delete();
                    gr_pion_x->Delete();
                    gr_pion_y->Delete();
                }
            }
        }
    }
    cout<<endl;
    cout<<"--> nRuns = "<<nRuns<<endl;

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_7->Write();
    h_8->Write();
    h_9->Write();
    h_10->Write();
    h_11->Write();
    h_12->Write();
    h_13->Write();
    h_14->Write();
    _file->Close();

    cout<<"--> Output file: "<<output_file_name<<endl;
}

void function_3(TString output_file_name)
{
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//

    Double_t* coord_x0          = new Double_t[mtrx_size];
    Double_t* coord_x           = new Double_t[mtrx_size];

    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    TRandom3* rnd = new TRandom3(seed);

    Double_t p, a, Charge_C = 1;

    TH3D* h_rp1_x = new TH3D("h_rp1_x","h_rp1_x",300,0,300,1000,-0.007,0.002,100,-1,1);
    TH3D* h_rp3_x = new TH3D("h_rp3_x","h_rp3_x",300,0,300,1000,-0.007,0.002,100,-1,1);

    const Int_t n_i = 100;
    const Int_t n_j = 100;

    cout<<endl;
    for(Int_t i = 0; i < n_i; i++)
    {
        p = rnd->Uniform(0.0,270.0);

        for(Int_t j = 0; j < n_j; j++)
        {
            if(i%1 == 0)
            {
                printf("\r--> BeginOfEvent: %3.1f %%",(double)100.0*i/n_i);
                fflush(stdout);
            }

            a = rnd->Uniform(-0.006,0.006);
            /*x*/     coord_x0[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
            /*x'*/    coord_x0[1] = Constants::_beamAngleInitialAtCryPosition+(a);// [rad]
            /*y*/     coord_x0[2] = 0.0;                // [m]
            /*y'*/    coord_x0[3] = 0.0;                // [rad]
            /*l*/     coord_x0[4] = 0.0;
            /*dp*/    coord_x0[5] = 0.0;
            TGraph* gr_x = new TGraph();
            TGraph* gr_y = new TGraph();
            passMagnets(coord_x0,coord_x,p,Charge_C,gr_x,gr_y,false);

            h_rp1_x->Fill(p,a,gr_x->Eval(Constants::_xrph_51937_ua9_pos));
            h_rp3_x->Fill(p,a,gr_x->Eval(Constants::_xrph_52202_ua9_pos));

            gr_x->Delete();
            gr_y->Delete();
        }
    }
    cout<<endl;

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");
    h_rp1_x->Write();
    h_rp3_x->Write();
    _file->Close();
    cout<<"--> Output file: "<<output_file_name<<endl;
}

Bool_t isChanneled(Double_t particleAngleOut)
{
    if(TMath::Abs(particleAngleOut) < TMath::Abs(Constants::_crystalAngle) + TMath::Abs(Constants::_channeledAngleRange))
    {
        if(TMath::Abs(particleAngleOut) > TMath::Abs(Constants::_crystalAngle) - TMath::Abs(Constants::_channeledAngleRange))
        {
            return true;
        }
    }

    return false;
}

void passMagnets(Double_t* coord_x0, Double_t* coord_x, Double_t p, Double_t Charge_C, TGraph *gr_x, TGraph *gr_y, Bool_t print)
{
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;
    Int_t order = Constants::_order_transport_matrix;

    Double_t KQ1 = Constants::_grad_quad_1;      // [m-2]
    Double_t KQ2 = Constants::_grad_quad_2;      // [m-2]
    Double_t KQ3 = Constants::_grad_quad_3;      // [m-2]
    Double_t KQ4 = Constants::_grad_quad_4;      // [m-2]
    Double_t KS1 = Constants::_grad_sext_1;      // [m-3]
    Double_t ADL = Constants::_angle_dipole;     // [rad]
    Double_t QL = Constants::_quadrupole_length; // [m]
    Double_t DL = Constants::_dipole_length;     // [m]
    Double_t SL = Constants::_sextupole_length;  // [m]
    Double_t P0 = Constants::_nominal_momentum;  // [GeV/c]
    Double_t APH = Constants::_aph;              // [m]
    Double_t APV = Constants::_apv;              // [m]

    Double_t DRIFTL0 = Constants::_xrph_52202_ua9_pos - Constants::_cry3_51799_ua9_pos;   // [m]
    Double_t DRIFTL1 = Constants::_q1_51810_pos-Constants::_cry3_51799_ua9_pos;           // [m]
    Double_t DRIFTL2 = Constants::_q2_51910_pos-Constants::_q1_51810_pos;                 // [m]
    Double_t DRIFTL3 = Constants::_xrph_51937_ua9_pos-Constants::_q2_51910_pos;           // [m]
    Double_t DRIFTL4 = Constants::_lsf_52005_pos-Constants::_xrph_51937_ua9_pos;          // [m]
    Double_t DRIFTL5 = Constants::_q3_52010_pos-Constants::_lsf_52005_pos;                // [m]
    Double_t DRIFTL6 = Constants::_mba_52030_pos-(DL/2.0)-Constants::_q3_52010_pos;       // [m]
    Double_t DRIFTL7 = Constants::_mba_52050_pos-DL-Constants::_mba_52030_pos;            // [m]
    Double_t DRIFTL8 = Constants::_mbb_52070_pos-DL-Constants::_mba_52050_pos;            // [m]
    Double_t DRIFTL9 = Constants::_mbb_52090_pos-DL-Constants::_mbb_52070_pos;            // [m]
    Double_t DRIFTL10 = Constants::_q4_52110_pos-(DL/2.0)-Constants::_mbb_52090_pos;       // [m]
    Double_t DRIFTL11 = Constants::_mbb_52130_pos-(DL/2.0)-Constants::_q4_52110_pos;      // [m]
    Double_t DRIFTL12 = Constants::_mbb_52150_pos-DL-Constants::_mbb_52130_pos;           // [m]
    Double_t DRIFTL13 = Constants::_xrph_52202_ua9_pos-(DL/2.0)-Constants::_mbb_52150_pos;// [m]

    //-----------------------------------------//

    Double_t* coord_temp_1_x    = new Double_t[mtrx_size];
    Double_t* coord_temp_2_x    = new Double_t[mtrx_size];

    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    Int_t gr_x_ipoint = 0, gr_y_ipoint = 0;
    // INITIAL
    if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos);
    gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[0]);
    gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[2]);
    gr_x_ipoint++;
    gr_y_ipoint++;

    if(Constants::_switch_magnets && Charge_C != 0)
    {
        // DRIFT1
        if(!magnet->GetNewCoordDrift(DRIFTL1,order,coord_x0,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}

        // QUAD1: FOC in X, DEFOC in Y
        if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ1,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_q1_51810_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q1_51810_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT2
        if(!magnet->GetNewCoordDrift(DRIFTL2,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}

        // QUAD2: FOC in Y, DEFOC in X
        if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ2,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_q2_51910_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q2_51910_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT3
        if(!magnet->GetNewCoordDrift(DRIFTL3,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_xrph_51937_ua9_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT4
        if(!magnet->GetNewCoordDrift(DRIFTL4,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}

        // SEXT1
        if(!magnet->GetNewCoordSextupole(Charge_C*KS1,p,P0,SL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_lsf_52005_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_lsf_52005_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT5
        if(!magnet->GetNewCoordDrift(DRIFTL5,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}

        // QUAD3: FOC in X, DEFOC in Y
        if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ3,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_q3_52010_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q3_52010_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT6
        if(!magnet->GetNewCoordDrift(DRIFTL6,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q3_52010_pos+DRIFTL6);

        // DIPOLE1: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mba_52030_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52030_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT7
        if(!magnet->GetNewCoordDrift(DRIFTL7,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mba_52030_pos+DL/2+DRIFTL7);

        // DIPOLE2: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mba_52050_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52050_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT8
        if(!magnet->GetNewCoordDrift(DRIFTL8,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mba_52050_pos+DL/2+DRIFTL8);

        // DIPOLE3: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mbb_52070_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52070_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT9
        if(!magnet->GetNewCoordDrift(DRIFTL9,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52070_pos+DL/2+DRIFTL9);

        // DIPOLE4: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mbb_52090_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52090_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT10
        if(!magnet->GetNewCoordDrift(DRIFTL10,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52090_pos+DL/2+DRIFTL10);

        // QUAD4: FOC in Y, DEFOC in X
        if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ4,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_q4_52110_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q4_52110_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT11
        if(!magnet->GetNewCoordDrift(DRIFTL11,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q4_52110_pos+DRIFTL11);

        // DIPOLE5: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mbb_52130_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52130_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT12
        if(!magnet->GetNewCoordDrift(DRIFTL12,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52130_pos+DL/2+DRIFTL12);

        // DIPOLE6: BEND & DEFOC in X, FOC in Y
        if(!magnet->GetNewCoordDipole(ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_mbb_52150_pos<<" [m] ##"<<endl; /*continue;*/}
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52150_pos+DL/2);
        gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[0]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[2]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        // DRIFT13
        if(!magnet->GetNewCoordDrift(DRIFTL13,order,coord_temp_1_x,coord_x,APH,APV))
        {if(print) cout<<"## At : "<<Constants::_xrph_52202_ua9_pos<<" [m] ##"<<endl; /*continue;*/}        
    }
    else
    {
        // DRIFT0
        if(!magnet->GetNewCoordDrift(DRIFTL0,order,coord_x0,coord_x,APH,APV))
        {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
    }

    // FINAL
    if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x[0],coord_x[1],Constants::_xrph_52202_ua9_pos);
    gr_x->SetPoint(gr_x_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[0]);
    gr_y->SetPoint(gr_y_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[2]);
}
