//my
#include "./src/MagClass.h"
#include "./src/DecayClass.h"
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
void function_4(TString output_file_name);
void function_5(TString input_filedir_name, TString output_file_name);
void function_6(Double_t crystalOrientation);

Bool_t isChanneled(Double_t particleAngleOut);
void passMagnets(Double_t s, Double_t* coord_x0, Double_t* coord_x, Double_t p, Double_t Charge_C, TGraph* gr_x, TGraph* gr_y, TGraph *gr_xp, TGraph *gr_yp, Bool_t print);
Double_t getGamma(Double_t mass, Double_t pmag);
Double_t getBeta(Double_t gamma);

Int_t main(int argc, char* argv[])
{
    if(argc == 5)
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
        case 4:
          function_4(argv[2]);
          break;
        case 5:
          function_5(argv[1],argv[2]);
          break;
        case 6:
          function_6(atof(argv[4]));
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
        cout<<"--> 4 -- function_4() : Lc products trajectories (random three-body decay)"<<endl;
        cout<<"--> 5 -- function_5() : Lc products trajectories (random two-body decay)"<<endl;
        cout<<"--> 6 -- function_6() : Lc products trajectories (random two-body decay, fixed momentum)"<<endl;
        cout<<endl;
        cout<<"--> ERROR:: Wrong imput parameters number:"<<endl<<
              "--> [0] -- script name"<<endl<<
              "--> [1] -- input filename"<<endl<<
              "--> [2] -- output filename"<<endl<<
              "--> [3] -- function id"<<endl<<
              "--> [4] -- crystal orientation angle [rad]"<<endl;
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
    TGraph* gr_xp = new TGraph();
    gr_x->SetName("gr_x");
    gr_xp->SetName("gr_xp");
    TGraph* gr_y = new TGraph();
    TGraph* gr_yp = new TGraph();
    gr_y->SetName("gr_y");
    gr_yp->SetName("gr_yp");

    passMagnets(0,coord_x0,coord_x,p,Charge_C,gr_x,gr_y,gr_xp,gr_yp,true);

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");
    gr_x->Write();
    gr_xp->Write();
    gr_y->Write();
    gr_yp->Write();
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
            Double_t thetaXLambdaCout = TMath::ATan(_Px_Lc_2/_Pz_Lc_2) - TMath::ATan(_Px_Lc_1/_Pz_Lc_1);

            h_2->Fill(thetaXLambdaCout);

            if(isChanneled(thetaXLambdaCout))
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
                    TGraph* gr_proton_xp = new TGraph();
                    TGraph* gr_proton_y = new TGraph();
                    TGraph* gr_proton_yp = new TGraph();
                    TString gr_proton_x_name = "gr_proton_x_";
                    TString gr_proton_xp_name = "gr_proton_xp_";
                    TString gr_proton_y_name = "gr_proton_y_";
                    TString gr_proton_yp_name = "gr_proton_yp_";
                    gr_proton_x_name += nRuns;
                    gr_proton_xp_name += nRuns;
                    gr_proton_y_name += nRuns;
                    gr_proton_yp_name += nRuns;
                    gr_proton_x->SetName(gr_proton_x_name.Data());
                    gr_proton_xp->SetName(gr_proton_xp_name.Data());
                    gr_proton_y->SetName(gr_proton_y_name.Data());
                    gr_proton_yp->SetName(gr_proton_yp_name.Data());

                    passMagnets(0,coord_x0_proton,coord_x_proton,p_proton,charge_proton,gr_proton_x,gr_proton_y,gr_proton_xp,gr_proton_yp,false);

                    TGraph* gr_kaon_x = new TGraph();
                    TGraph* gr_kaon_xp = new TGraph();
                    TGraph* gr_kaon_y = new TGraph();
                    TGraph* gr_kaon_yp = new TGraph();
                    TString gr_kaon_x_name = "gr_kaon_x_";
                    TString gr_kaon_xp_name = "gr_kaon_xp_";
                    TString gr_kaon_y_name = "gr_kaon_y_";
                    TString gr_kaon_yp_name = "gr_kaon_yp_";
                    gr_kaon_x_name += nRuns;
                    gr_kaon_xp_name += nRuns;
                    gr_kaon_y_name += nRuns;
                    gr_kaon_yp_name += nRuns;
                    gr_kaon_x->SetName(gr_kaon_x_name.Data());
                    gr_kaon_xp->SetName(gr_kaon_xp_name.Data());
                    gr_kaon_y->SetName(gr_kaon_y_name.Data());
                    gr_kaon_yp->SetName(gr_kaon_yp_name.Data());

                    passMagnets(0,coord_x0_kaon,coord_x_kaon,p_kaon,charge_kaon,gr_kaon_x,gr_kaon_y,gr_kaon_xp,gr_kaon_yp,false);

                    TGraph* gr_pion_x = new TGraph();
                    TGraph* gr_pion_xp = new TGraph();
                    TGraph* gr_pion_y = new TGraph();
                    TGraph* gr_pion_yp = new TGraph();
                    TString gr_pion_x_name = "gr_pion_x_";
                    TString gr_pion_xp_name = "gr_pion_xp_";
                    TString gr_pion_y_name = "gr_pion_y_";
                    TString gr_pion_yp_name = "gr_pion_yp_";
                    gr_pion_x_name += nRuns;
                    gr_pion_xp_name += nRuns;
                    gr_pion_y_name += nRuns;
                    gr_pion_yp_name += nRuns;
                    gr_pion_x->SetName(gr_pion_x_name.Data());
                    gr_pion_xp->SetName(gr_pion_xp_name.Data());
                    gr_pion_y->SetName(gr_pion_y_name.Data());
                    gr_pion_yp->SetName(gr_pion_yp_name.Data());

                    passMagnets(0,coord_x0_pion,coord_x_pion,p_pion,charge_pion,gr_pion_x,gr_pion_y,gr_pion_xp,gr_pion_yp,false);

                    nRuns++;

                    gr_proton_x->Write();
                    gr_proton_y->Write();
                    gr_kaon_x->Write();
                    gr_kaon_y->Write();
                    gr_pion_x->Write();
                    gr_pion_y->Write();

                    gr_proton_xp->Write();
                    gr_proton_yp->Write();
                    gr_kaon_xp->Write();
                    gr_kaon_yp->Write();
                    gr_pion_xp->Write();
                    gr_pion_yp->Write();

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

                    gr_proton_xp->Delete();
                    gr_proton_yp->Delete();
                    gr_kaon_xp->Delete();
                    gr_kaon_yp->Delete();
                    gr_pion_xp->Delete();
                    gr_pion_yp->Delete();
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
            TGraph* gr_xp = new TGraph();
            TGraph* gr_y = new TGraph();
            TGraph* gr_yp = new TGraph();
            passMagnets(0,coord_x0,coord_x,p,Charge_C,gr_x,gr_y,gr_xp,gr_yp,false);

            h_rp1_x->Fill(p,a,gr_x->Eval(Constants::_xrph_51937_ua9_pos));
            h_rp3_x->Fill(p,a,gr_x->Eval(Constants::_xrph_52202_ua9_pos));

            gr_x->Delete();
            gr_xp->Delete();
            gr_y->Delete();
            gr_yp->Delete();
        }
    }
    cout<<endl;

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");
    h_rp1_x->Write();
    h_rp3_x->Write();
    _file->Close();
    cout<<"--> Output file: "<<output_file_name<<endl;
}

void function_4(TString output_file_name)
{
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//
    // Decay section
    //-----------------------------------------//
    DecayClass* decay = new DecayClass();
    //-----------------------------------------//

    Double_t p_proton = -999, charge_proton = -999, theta_x_proton = -999, theta_y_proton = -999;
    Double_t* coord_x0_proton          = new Double_t[mtrx_size];
    Double_t* coord_x_proton           = new Double_t[mtrx_size];

    Double_t p_kaon = -999, charge_kaon = -999, theta_x_kaon = -999, theta_y_kaon = -999;
    Double_t* coord_x0_kaon          = new Double_t[mtrx_size];
    Double_t* coord_x_kaon           = new Double_t[mtrx_size];

    Double_t p_pion = -999, charge_pion = -999, theta_x_pion = -999, theta_y_pion = -999;
    Double_t* coord_x0_pion          = new Double_t[mtrx_size];
    Double_t* coord_x_pion           = new Double_t[mtrx_size];

    Int_t nRuns = 0;
    const Int_t nLc = 200;// number of decayed Lc

    TLorentzVector pProd[4];
    Double_t mProd[4];
    Double_t s_decay = 0;

    // Mass
    mProd[0] = 2.28646; // Lc       [GeV/c2]
    mProd[1] = 0.93827; // Proton   [GeV/c2]
    mProd[2] = 0.49368; // Kaon     [GeV/c2]
    mProd[3] = 0.13957; // Pion     [GeV/c2]

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");

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

    for(Int_t i = 0; i < nLc; i++)
    {
        Double_t pLc = 250.0; // [GeV/c]
        Double_t eLc = TMath::Sqrt(pLc*pLc + mProd[0]*mProd[0]); // [GeV]

        pProd[0].SetPxPyPzE(0,0,pLc,eLc);

        if(decay->threeBody(pProd,mProd))
        {
            //--------------------//
            // Proton
            //--------------------//

            p_proton = pProd[1].P();
//            cout<<"p_proton = "<<p_proton<<endl;
            theta_x_proton = TMath::ATan(pProd[1].Px()/pProd[1].Pz());
            theta_y_proton = TMath::ATan(pProd[1].Py()/pProd[1].Pz());
            charge_proton = 1.0;
            /*x*/     coord_x0_proton[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
            /*x'*/    coord_x0_proton[1] = Constants::_beamAngleInitialAtCryPosition + Constants::_crystalAngle + theta_x_proton;// [rad]
            /*y*/     coord_x0_proton[2] = 0.0;             // [m]
            /*y'*/    coord_x0_proton[3] = theta_y_proton;  // [rad]
            /*l*/     coord_x0_proton[4] = 0.0;
            /*dp*/    coord_x0_proton[5] = 0.0;

            //--------------------//
            // Kaon
            //--------------------//

            p_kaon = pProd[2].P();
//            cout<<"p_kaon = "<<p_kaon<<endl;
            theta_x_kaon = TMath::ATan(pProd[2].Px()/pProd[2].Pz());
            theta_y_kaon = TMath::ATan(pProd[2].Py()/pProd[2].Pz());
            charge_kaon = -1.0;
            /*x*/     coord_x0_kaon[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
            /*x'*/    coord_x0_kaon[1] = Constants::_beamAngleInitialAtCryPosition + Constants::_crystalAngle + theta_x_kaon;// [rad]
            /*y*/     coord_x0_kaon[2] = 0.0;           // [m]
            /*y'*/    coord_x0_kaon[3] = theta_y_kaon;  // [rad]
            /*l*/     coord_x0_kaon[4] = 0.0;
            /*dp*/    coord_x0_kaon[5] = 0.0;

            //--------------------//
            // Pion
            //--------------------//

            p_pion = pProd[3].P();
//            cout<<"p_pion = "<<p_pion<<endl;
            theta_x_pion = TMath::ATan(pProd[3].Px()/pProd[3].Pz());
            theta_y_pion = TMath::ATan(pProd[3].Py()/pProd[3].Pz());
            charge_pion = 1.0;
            /*x*/     coord_x0_pion[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
            /*x'*/    coord_x0_pion[1] = Constants::_beamAngleInitialAtCryPosition + Constants::_crystalAngle + theta_x_pion;// [rad]
            /*y*/     coord_x0_pion[2] = 0.0;           // [m]
            /*y'*/    coord_x0_pion[3] = theta_y_pion;  // [rad]
            /*l*/     coord_x0_pion[4] = 0.0;
            /*dp*/    coord_x0_pion[5] = 0.0;

            //--------------------//
            // Magnets
            //--------------------//

            h_9->Fill(p_proton);
            h_10->Fill(p_kaon);
            h_11->Fill(p_pion);
            h_12->Fill(theta_x_proton);
            h_13->Fill(theta_x_kaon);
            h_14->Fill(theta_x_pion);

            TGraph* gr_proton_x = new TGraph();
            TGraph* gr_proton_y = new TGraph();
            TString gr_proton_x_name = "gr_proton_x_";
            TString gr_proton_y_name = "gr_proton_y_";
            gr_proton_x_name += nRuns;
            gr_proton_y_name += nRuns;
            gr_proton_x->SetName(gr_proton_x_name.Data());
            gr_proton_y->SetName(gr_proton_y_name.Data());

            TGraph* gr_proton_xp = new TGraph();
            TGraph* gr_proton_yp = new TGraph();
            TString gr_proton_xp_name = "gr_proton_xp_";
            TString gr_proton_yp_name = "gr_proton_yp_";
            gr_proton_xp_name += nRuns;
            gr_proton_yp_name += nRuns;
            gr_proton_xp->SetName(gr_proton_xp_name.Data());
            gr_proton_yp->SetName(gr_proton_yp_name.Data());

            passMagnets(s_decay,coord_x0_proton,coord_x_proton,p_proton,charge_proton,gr_proton_x,gr_proton_y,gr_proton_xp,gr_proton_yp,false);

            TGraph* gr_kaon_x = new TGraph();
            TGraph* gr_kaon_y = new TGraph();
            TString gr_kaon_x_name = "gr_kaon_x_";
            TString gr_kaon_y_name = "gr_kaon_y_";
            gr_kaon_x_name += nRuns;
            gr_kaon_y_name += nRuns;
            gr_kaon_x->SetName(gr_kaon_x_name.Data());
            gr_kaon_y->SetName(gr_kaon_y_name.Data());

            TGraph* gr_kaon_xp = new TGraph();
            TGraph* gr_kaon_yp = new TGraph();
            TString gr_kaon_xp_name = "gr_kaon_xp_";
            TString gr_kaon_yp_name = "gr_kaon_yp_";
            gr_kaon_xp_name += nRuns;
            gr_kaon_yp_name += nRuns;
            gr_kaon_xp->SetName(gr_kaon_xp_name.Data());
            gr_kaon_yp->SetName(gr_kaon_yp_name.Data());

            passMagnets(s_decay,coord_x0_kaon,coord_x_kaon,p_kaon,charge_kaon,gr_kaon_x,gr_kaon_y,gr_kaon_xp,gr_kaon_yp,false);

            TGraph* gr_pion_x = new TGraph();
            TGraph* gr_pion_y = new TGraph();
            TString gr_pion_x_name = "gr_pion_x_";
            TString gr_pion_y_name = "gr_pion_y_";
            gr_pion_x_name += nRuns;
            gr_pion_y_name += nRuns;
            gr_pion_x->SetName(gr_pion_x_name.Data());
            gr_pion_y->SetName(gr_pion_y_name.Data());

            TGraph* gr_pion_xp = new TGraph();
            TGraph* gr_pion_yp = new TGraph();
            TString gr_pion_xp_name = "gr_pion_xp_";
            TString gr_pion_yp_name = "gr_pion_yp_";
            gr_pion_xp_name += nRuns;
            gr_pion_yp_name += nRuns;
            gr_pion_xp->SetName(gr_pion_xp_name.Data());
            gr_pion_yp->SetName(gr_pion_yp_name.Data());

            passMagnets(s_decay,coord_x0_pion,coord_x_pion,p_pion,charge_pion,gr_pion_x,gr_pion_y,gr_pion_xp,gr_pion_yp,false);

            nRuns++;

            gr_proton_x->Write();
            gr_proton_y->Write();
            gr_kaon_x->Write();
            gr_kaon_y->Write();
            gr_pion_x->Write();
            gr_pion_y->Write();

            gr_proton_xp->Write();
            gr_proton_yp->Write();
            gr_kaon_xp->Write();
            gr_kaon_yp->Write();
            gr_pion_xp->Write();
            gr_pion_yp->Write();

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

            gr_proton_xp->Delete();
            gr_proton_yp->Delete();
            gr_kaon_xp->Delete();
            gr_kaon_yp->Delete();
            gr_pion_xp->Delete();
            gr_pion_yp->Delete();
        }
    }

    cout<<endl;
    cout<<"--> nRuns = "<<nRuns<<endl;

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

void function_5(TString input_filedir_name, TString output_file_name)
{
    //-----------------------------------------//
    Bool_t plot_graph   = kFALSE;//kTRUE;               // to plot graphs
    const Int_t nLc     = 1e8;                // number of maximum Lc decays
    Int_t detected = 0;
    Int_t nLc_reconstructed = 0;
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
    for(Int_t i = 1; i <= 20; i++)
    {
        input_file_name = input_filedir_name;
        input_file_name += "crystaltarget_";
        input_file_name += i;
        input_file_name += ".root";
        if(fChain1->Add(input_file_name.Data()))
        {
            cout<<"--> File: "<<input_file_name<<" has been added"<<endl;
        }
        else
        {
            assert(0);
        }
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

    Long64_t nEntries = fChain1->GetEntries();
    cout<<"--> nEntries: "<<nEntries<<endl;

    //-----------------------------------------//
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    TRandom3* rnd = new TRandom3(seed);
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//
    // Decay section
    //-----------------------------------------//
    DecayClass* decay = new DecayClass();
    //-----------------------------------------//

    Double_t p_lambda0 = -999, charge_lambda0 = -999, theta_x_lambda0 = -999, theta_y_lambda0 = -999;
    Double_t* coord_x0_lambda0          = new Double_t[mtrx_size];
    Double_t* coord_x_lambda0           = new Double_t[mtrx_size];

    Double_t p_pion_p = -999, charge_pion_p = -999, theta_x_pion_p = -999, theta_y_pion_p = -999;
    Double_t* coord_x0_pion_p          = new Double_t[mtrx_size];
    Double_t* coord_x_pion_p           = new Double_t[mtrx_size];

    Double_t p_proton = -999, charge_proton = -999, theta_x_proton = -999, theta_y_proton = -999;
    Double_t* coord_x0_proton          = new Double_t[mtrx_size];
    Double_t* coord_x_proton           = new Double_t[mtrx_size];

    Double_t p_pion_m = -999, charge_pion_m = -999, theta_x_pion_m = -999, theta_y_pion_m = -999;
    Double_t* coord_x0_pion_m          = new Double_t[mtrx_size];
    Double_t* coord_x_pion_m           = new Double_t[mtrx_size];

    Int_t nRuns = 0;

    TLorentzVector pProd_1[3];
    TLorentzVector pProd_2[3];
    Double_t mProd_1[3];
    Double_t mProd_2[3];

    Double_t s_decay_1 = 0;
    Double_t s_decay_2 = 0;

    // Mass
    mProd_1[0] = 2.28646; // Lc+    [GeV/c2]
    mProd_1[1] = 1.11568; // L0     [GeV/c2]
    mProd_1[2] = 0.13957; // Pion+  [GeV/c2]

    mProd_2[0] = 1.11568; // L0     [GeV/c2]
    mProd_2[1] = 0.93827; // Proton [GeV/c2]
    mProd_2[2] = 0.13957; // Pion-  [GeV/c2]

    TFile* _file = new TFile(output_file_name.Data(),"RECREATE");

    TH1D* h_1 = new TH1D("h_1","L_{c} momentum channeled [GeV/c]",330,-30.0,300.0);
    TH1D* h_2 = new TH1D("h_2","L_{0} momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_3 = new TH1D("h_3","pion+ momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_4 = new TH1D("h_4","proton momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_5 = new TH1D("h_5","pion- momentum [GeV/c]",330,-30.0,300.0);
    TH1D* h_6 = new TH1D("h_6","L_{0} angle",4000000,-4.0,4.0);
    TH1D* h_7 = new TH1D("h_7","pion+ angle",4000000,-4.0,4.0);
    TH1D* h_8 = new TH1D("h_8","proton angle",4000000,-4.0,4.0);
    TH1D* h_9 = new TH1D("h_9","pion- angle",4000000,-4.0,4.0);
    TH2D* h_10 = new TH2D("h_10","XY L_{0} on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_11 = new TH2D("h_11","XY pion+ on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_12 = new TH2D("h_12","XY proton on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_13 = new TH2D("h_13","XY pion- on RP1",1000,-10,10,1000,-10,10);
    TH2D* h_14 = new TH2D("h_14","XY L_{0} on RP3",1000,-10,10,1000,-10,10);
    TH2D* h_15 = new TH2D("h_15","XY pion+ on RP3",1000,-10,10,1000,-10,10);
    TH2D* h_16 = new TH2D("h_16","XY proton on RP3",1000,-10,10,1000,-10,10);
    TH2D* h_17 = new TH2D("h_17","XY pion- on RP3",1000,-10,10,1000,-10,10);
    TH1D* h_18 = new TH1D("h_18","L_{0} decay length [m]",110000,-10,100);
    TH1D* h_19 = new TH1D("h_19","L_{0} decay length (z-proj) [m]",110000,-10,100);
    TH1D* h_20 = new TH1D("h_20","L_{c} angle out [rad]",4000000,-4.0,4.0);
    TH1D* h_21 = new TH1D("h_21","L_{c} angle in [rad]",4000000,-4.0,4.0);
    TH1D* h_22 = new TH1D("h_22","L_{c} momentum in [GeV/c]",330,-30.0,300.0);
    TH2D* h_23 = new TH2D("h_23","L_{c} momentum [GeV/c] vs angle in [rad]",40000,-4.0,4.0,330,-30.0,300.0);
    TH2D* h_24 = new TH2D("h_24","L_{0} momentum [GeV/c] vs L_{c} momentum [GeV/c]",330,-30.0,300.0,330,-30.0,300.0);
    TH2D* h_25 = new TH2D("h_25","pion+ momentum [GeV/c] vs L_{c} momentum [GeV/c]",330,-30.0,300.0,330,-30.0,300.0);
    TH2D* h_26 = new TH2D("h_26","proton momentum [GeV/c] vs L_{c} momentum [GeV/c]",330,-30.0,300.0,330,-30.0,300.0);
    TH2D* h_27 = new TH2D("h_27","pion- momentum [GeV/c] vs L_{c} momentum [GeV/c]",330,-30.0,300.0,330,-30.0,300.0);
    TH2D* h_28 = new TH2D("h_28","XY L_{0} at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_29 = new TH2D("h_29","XY pion+ at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_30 = new TH2D("h_30","XY proton at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_31 = new TH2D("h_31","XY pion- at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_32 = new TH2D("h_32","XY L_{0} at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_33 = new TH2D("h_33","XY pion+ at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_34 = new TH2D("h_34","XY proton at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_35 = new TH2D("h_35","XY pion- at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_36 = new TH2D("h_36","XpYp L_{0} at TR1",4000,-4,4,4000,-4,4);
    TH2D* h_37 = new TH2D("h_37","XpYp pion+ at TR1",4000,-4,4,4000,-4,4);
    TH2D* h_38 = new TH2D("h_38","XpYp proton at TR1",4000,-4,4,4000,-4,4);
    TH2D* h_39 = new TH2D("h_39","XpYp pion- at TR1",4000,-4,4,4000,-4,4);
    TH2D* h_40 = new TH2D("h_40","XpYp L_{0} at TR2",4000,-4,4,4000,-4,4);
    TH2D* h_41 = new TH2D("h_41","XpYp pion+ at TR2",4000,-4,4,4000,-4,4);
    TH2D* h_42 = new TH2D("h_42","XpYp proton at TR2",4000,-4,4,4000,-4,4);
    TH2D* h_43 = new TH2D("h_43","XpYp pion- at TR2",4000,-4,4,4000,-4,4);
    TH2D* h_44 = new TH2D("h_44","XpX L_{0} at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_45 = new TH2D("h_45","XpX pion+ at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_46 = new TH2D("h_46","XpX proton at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_47 = new TH2D("h_47","XpX pion- at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_48 = new TH2D("h_48","XpX L_{0} at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_49 = new TH2D("h_49","XpX pion+ at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_50 = new TH2D("h_50","XpX proton at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_51 = new TH2D("h_51","XpX pion- at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_52 = new TH2D("h_52","YpY L_{0} at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_53 = new TH2D("h_53","YpY pion+ at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_54 = new TH2D("h_54","YpY proton at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_55 = new TH2D("h_55","YpY pion- at TR1",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_56 = new TH2D("h_56","YpY L_{0} at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_57 = new TH2D("h_57","YpY pion+ at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_58 = new TH2D("h_58","YpY proton at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_59 = new TH2D("h_59","YpY pion- at TR2",4000,-0.2,0.2,6000,-0.3,0.3);
    TH2D* h_60 = new TH2D("h_60","proton angle vs pion- angle at L_{0} decay",6000,-0.3,0.3,6000,-0.3,0.3);
    TH2D* h_61 = new TH2D("h_61","XY pion+ at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_62 = new TH2D("h_62","XY proton at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_63 = new TH2D("h_63","XY pion- at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_64 = new TH2D("h_64","XY pion+ at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_65 = new TH2D("h_65","XY proton at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_66 = new TH2D("h_66","XY pion- at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
    TH2D* h_67 = new TH2D("h_67","L_{0} theta vs phi",800,-4.0,4.0,20000,-0.1,0.1);
    TH2D* h_68 = new TH2D("h_68","pion+ theta vs phi",800,-4.0,4.0,20000,-0.1,0.1);
    TH2D* h_69 = new TH2D("h_69","proton theta vs phi (cut)",800,-4.0,4.0,20000,-0.1,0.1);
    TH2D* h_70 = new TH2D("h_70","pion- theta vs phi",800,-4.0,4.0,20000,-0.1,0.1);
    TH2D* h_71 = new TH2D("h_71","L_{c} momentum [GeV/c] vs angle in Y [rad]",40000,-4.0,4.0,330,-30.0,300.0);

    cout<<endl;
    Bool_t isLambdacgenerated;
    Double_t pLc;
    Double_t thetaXLambdaCout, thetaXLambdaCin, thetaYLambdaCout;

    for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        fChain1->GetEntry(iEntry);

        if(iEntry%10000 == 0)
        {
            printf("\r--> Working: %3.1f | nRuns = %d",100.0*iEntry/nEntries,nRuns);
            fflush(stdout);
        }

        if(_isGenerated && _isTarget && _isCrystal)
        {
            thetaXLambdaCout = TMath::ATan(_Px_Lc_2/_Pz_Lc_2) - TMath::ATan(_Px_Lc_1/_Pz_Lc_1);            

            h_20->Fill(thetaXLambdaCout);            

            isLambdacgenerated = kFALSE;

            for(Int_t iParticle = 0; iParticle < _Nparticle; iParticle++)
            {
                if(_ID[iParticle] == 4122 && !isLambdacgenerated) // lambda_c+
                {
                    isLambdacgenerated = kTRUE;

                    pLc = TMath::Sqrt(_Px[iParticle]*_Px[iParticle]+_Py[iParticle]*_Py[iParticle]+_Pz[iParticle]*_Pz[iParticle]);
                    thetaYLambdaCout    = TMath::ATan(_Py[iParticle]/_Pz[iParticle]);
                    thetaXLambdaCin     = TMath::ATan(_Px[iParticle]/_Pz[iParticle]);

                    h_21->Fill(thetaXLambdaCin);
                    h_22->Fill(pLc);
                    h_23->Fill(thetaXLambdaCin,pLc);
                    h_71->Fill(thetaYLambdaCout,pLc);
                }
            }

            Double_t eLc = TMath::Sqrt(pLc*pLc + mProd_1[0]*mProd_1[0]); // [GeV]

            pProd_1[0].SetPxPyPzE(0,0,pLc,eLc);

            detected = 0;

            if(/*isChanneled(thetaXLambdaCout) &&*/ isLambdacgenerated)
            {
                h_1->Fill(pLc);

                if(nRuns >= nLc) continue;

                if(decay->twoBody(pProd_1,mProd_1))
                {
                    //--------------------//
                    // Lambda0
                    //--------------------//

                    p_lambda0 = pProd_1[1].P();
                    h_2->Fill(p_lambda0);
                    h_24->Fill(pLc,p_lambda0);
                    h_67->Fill(pProd_1[1].Phi(),pProd_1[1].Theta());

                    theta_x_lambda0 =
                            Constants::_beamAngleInitialAtCryPosition +
                            Constants::_crystalAngle +
                            Constants::_crystalOrientation +
                            TMath::ATan(pProd_1[1].Px()/pProd_1[1].Pz());
                    h_6->Fill(theta_x_lambda0);
                    theta_y_lambda0 = thetaYLambdaCout + TMath::ATan(pProd_1[1].Py()/pProd_1[1].Pz());
                    charge_lambda0 = 0.0;
                    /*x*/     coord_x0_lambda0[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                    /*x'*/    coord_x0_lambda0[1] = theta_x_lambda0;// [rad]
                    /*y*/     coord_x0_lambda0[2] = 0.0;              // [m]
                    /*y'*/    coord_x0_lambda0[3] = theta_y_lambda0;  // [rad]
                    /*l*/     coord_x0_lambda0[4] = 0.0;
                    /*dp*/    coord_x0_lambda0[5] = 0.0;

                    // RANDOM //
                    Double_t gamma_L0   = getGamma(mProd_2[0],p_lambda0);
                    Double_t beta_L0    = getBeta(gamma_L0);
                    Double_t prob       = rnd->Uniform(0.0,1.0);
                    Double_t cosTheta = pProd_1[1].CosTheta();

                    if(prob == 0)
                    {
                        while( 1 )
                        {
                            prob = rnd->Uniform(0.0,1.0);
                            if(prob > 0) break;
                        }
                    }
                    s_decay_2 = Constants::_ctau_Lambda0*beta_L0*gamma_L0*TMath::Log(1.0/prob); // [m]
                    h_18->Fill(s_decay_2);
                    s_decay_2 = s_decay_2*cosTheta;
                    h_19->Fill(s_decay_2);

                    //--------------------//
                    // Pion+
                    //--------------------//

                    p_pion_p = pProd_1[2].P();
                    h_3->Fill(p_pion_p);
                    h_25->Fill(pLc,p_pion_p);
                    h_68->Fill(pProd_1[2].Phi(),pProd_1[2].Theta());

                    theta_x_pion_p =
                            Constants::_beamAngleInitialAtCryPosition +
                            Constants::_crystalAngle +
                            Constants::_crystalOrientation +
                            TMath::ATan(pProd_1[2].Px()/pProd_1[2].Pz());
                    h_7->Fill(theta_x_pion_p);
                    theta_y_pion_p = thetaYLambdaCout + TMath::ATan(pProd_1[2].Py()/pProd_1[2].Pz());
                    charge_pion_p = 1.0;
                    /*x*/     coord_x0_pion_p[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                    /*x'*/    coord_x0_pion_p[1] = theta_x_pion_p;// [rad]
                    /*y*/     coord_x0_pion_p[2] = 0.0;           // [m]
                    /*y'*/    coord_x0_pion_p[3] = theta_y_pion_p;  // [rad]
                    /*l*/     coord_x0_pion_p[4] = 0.0;
                    /*dp*/    coord_x0_pion_p[5] = 0.0;

                    //--------------------//
                    // Magnets
                    //--------------------//

                    TGraph* gr_lambda0_x = new TGraph();
                    TGraph* gr_lambda0_xp = new TGraph();
                    TGraph* gr_lambda0_y = new TGraph();
                    TGraph* gr_lambda0_yp = new TGraph();
                    TString gr_lambda0_x_name = "gr_lambda0_x_";
                    TString gr_lambda0_xp_name = "gr_lambda0_xp_";
                    TString gr_lambda0_y_name = "gr_lambda0_y_";
                    TString gr_lambda0_yp_name = "gr_lambda0_yp_";
                    gr_lambda0_x_name += nRuns;
                    gr_lambda0_xp_name += nRuns;
                    gr_lambda0_y_name += nRuns;
                    gr_lambda0_yp_name += nRuns;
                    gr_lambda0_x->SetName(gr_lambda0_x_name.Data());
                    gr_lambda0_xp->SetName(gr_lambda0_xp_name.Data());
                    gr_lambda0_y->SetName(gr_lambda0_y_name.Data());
                    gr_lambda0_yp->SetName(gr_lambda0_yp_name.Data());

                    if(charge_lambda0 != 0)
                    {
                        assert(0);
                    }

                    // INITIAL
                    gr_lambda0_x->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[0]);
                    gr_lambda0_xp->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[1]);
                    gr_lambda0_y->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[2]);
                    gr_lambda0_yp->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[3]);

                    // DRIFT & FINAL
                    magnet->GetNewCoordDrift(s_decay_2,Constants::_order_transport_matrix,coord_x0_lambda0,coord_x_lambda0,Constants::_aph,Constants::_apv);
                    gr_lambda0_x->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[0]);
                    gr_lambda0_xp->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[1]);
                    gr_lambda0_y->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[2]);
                    gr_lambda0_yp->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[3]);

                    TGraph* gr_pion_p_x = new TGraph();
                    TGraph* gr_pion_p_xp = new TGraph();
                    TGraph* gr_pion_p_y = new TGraph();
                    TGraph* gr_pion_p_yp = new TGraph();
                    TString gr_pion_p_x_name = "gr_pion_p_x_";
                    TString gr_pion_p_xp_name = "gr_pion_p_xp_";
                    TString gr_pion_p_y_name = "gr_pion_p_y_";
                    TString gr_pion_p_yp_name = "gr_pion_p_yp_";
                    gr_pion_p_x_name += nRuns;
                    gr_pion_p_xp_name += nRuns;
                    gr_pion_p_y_name += nRuns;
                    gr_pion_p_yp_name += nRuns;
                    gr_pion_p_x->SetName(gr_pion_p_x_name.Data());
                    gr_pion_p_xp->SetName(gr_pion_p_xp_name.Data());
                    gr_pion_p_y->SetName(gr_pion_p_y_name.Data());
                    gr_pion_p_yp->SetName(gr_pion_p_yp_name.Data());

                    passMagnets(s_decay_1,coord_x0_pion_p,coord_x_pion_p,p_pion_p,charge_pion_p,gr_pion_p_x,gr_pion_p_y,gr_pion_p_xp,gr_pion_p_yp,false);

                    if(plot_graph)
                    {
                        gr_lambda0_x->Write();
                        gr_lambda0_xp->Write();
                        gr_lambda0_y->Write();
                        gr_lambda0_yp->Write();
                        gr_pion_p_x->Write();
                        gr_pion_p_xp->Write();
                        gr_pion_p_y->Write();
                        gr_pion_p_yp->Write();
                    }

                    if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 > Constants::_xrph_51937_ua9_pos)
                        h_10->Fill(gr_lambda0_x->Eval(Constants::_xrph_51937_ua9_pos),gr_lambda0_y->Eval(Constants::_xrph_51937_ua9_pos));
                    if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 > Constants::_xrph_52202_ua9_pos)
                        h_14->Fill(gr_lambda0_x->Eval(Constants::_xrph_52202_ua9_pos),gr_lambda0_y->Eval(Constants::_xrph_52202_ua9_pos));

                    if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 > Constants::_tr1_s_pos)
                    {
                        h_28->Fill(gr_lambda0_x->Eval(Constants::_tr1_s_pos),gr_lambda0_y->Eval(Constants::_tr1_s_pos));
                        h_36->Fill(gr_lambda0_xp->Eval(Constants::_tr1_s_pos),gr_lambda0_yp->Eval(Constants::_tr1_s_pos));
                        h_44->Fill(gr_lambda0_x->Eval(Constants::_tr1_s_pos),gr_lambda0_xp->Eval(Constants::_tr1_s_pos));
                        h_52->Fill(gr_lambda0_y->Eval(Constants::_tr1_s_pos),gr_lambda0_yp->Eval(Constants::_tr1_s_pos));
                    }

                    if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 > Constants::_tr2_s_pos)
                    {
                        h_32->Fill(gr_lambda0_x->Eval(Constants::_tr2_s_pos),gr_lambda0_y->Eval(Constants::_tr2_s_pos));
                        h_40->Fill(gr_lambda0_xp->Eval(Constants::_tr2_s_pos),gr_lambda0_yp->Eval(Constants::_tr2_s_pos));
                        h_48->Fill(gr_lambda0_x->Eval(Constants::_tr2_s_pos),gr_lambda0_xp->Eval(Constants::_tr2_s_pos));
                        h_56->Fill(gr_lambda0_y->Eval(Constants::_tr2_s_pos),gr_lambda0_yp->Eval(Constants::_tr2_s_pos));
                    }

                    h_11->Fill(gr_pion_p_x->Eval(Constants::_xrph_51937_ua9_pos),gr_pion_p_y->Eval(Constants::_xrph_51937_ua9_pos));
                    h_15->Fill(gr_pion_p_x->Eval(Constants::_xrph_52202_ua9_pos),gr_pion_p_y->Eval(Constants::_xrph_52202_ua9_pos));
                    h_29->Fill(gr_pion_p_x->Eval(Constants::_tr1_s_pos),gr_pion_p_y->Eval(Constants::_tr1_s_pos));
                    h_37->Fill(gr_pion_p_xp->Eval(Constants::_tr1_s_pos),gr_pion_p_yp->Eval(Constants::_tr1_s_pos));
                    h_45->Fill(gr_pion_p_x->Eval(Constants::_tr1_s_pos),gr_pion_p_xp->Eval(Constants::_tr1_s_pos));
                    h_53->Fill(gr_pion_p_y->Eval(Constants::_tr1_s_pos),gr_pion_p_yp->Eval(Constants::_tr1_s_pos));
                    h_33->Fill(gr_pion_p_x->Eval(Constants::_tr2_s_pos),gr_pion_p_y->Eval(Constants::_tr2_s_pos));
                    h_41->Fill(gr_pion_p_xp->Eval(Constants::_tr2_s_pos),gr_pion_p_yp->Eval(Constants::_tr2_s_pos));
                    h_49->Fill(gr_pion_p_x->Eval(Constants::_tr2_s_pos),gr_pion_p_xp->Eval(Constants::_tr2_s_pos));
                    h_57->Fill(gr_pion_p_y->Eval(Constants::_tr2_s_pos),gr_pion_p_yp->Eval(Constants::_tr2_s_pos));

                    // TRACKER 1
                    if(gr_pion_p_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                    {
                        detected += 1;
                        h_61->Fill(gr_pion_p_x->Eval(Constants::_tr1_s_pos),gr_pion_p_y->Eval(Constants::_tr1_s_pos));

                        // TRACKER 2
                        if(gr_pion_p_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_p_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_p_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_p_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                        {
                            detected += 1;
                            h_64->Fill(gr_pion_p_x->Eval(Constants::_tr2_s_pos),gr_pion_p_y->Eval(Constants::_tr2_s_pos));
                        }
                    }

                    pProd_2[0].SetPxPyPzE(0.0,0.0,pProd_1[1].P(),pProd_1[1].E());

                    if(decay->twoBody(pProd_2,mProd_2))
                    {
                        //--------------------//
                        // Proton
                        //--------------------//

                        p_proton = pProd_2[1].P();
                        h_4->Fill(p_proton);
                        h_26->Fill(pLc,p_proton);
                        h_69->Fill(pProd_2[1].Phi(),pProd_2[1].Theta());

                        theta_x_proton = coord_x_lambda0[1] + TMath::ATan(pProd_2[1].Px()/pProd_2[1].Pz());
                        h_8->Fill(theta_x_proton);
                        theta_y_proton = coord_x_lambda0[3] + TMath::ATan(pProd_2[1].Py()/pProd_2[1].Pz());
                        charge_proton = 1.0;
                        /*x*/     coord_x0_proton[0] = coord_x_lambda0[0];  // [m]
                        /*x'*/    coord_x0_proton[1] = theta_x_proton;      // [rad]
                        /*y*/     coord_x0_proton[2] = coord_x_lambda0[2];  // [m]
                        /*y'*/    coord_x0_proton[3] = theta_y_proton;      // [rad]
                        /*l*/     coord_x0_proton[4] = 0.0;
                        /*dp*/    coord_x0_proton[5] = 0.0;

                        //--------------------//
                        // Pion-
                        //--------------------//

                        p_pion_m = pProd_2[2].P();
                        h_5->Fill(p_pion_m);
                        h_27->Fill(pLc,p_pion_m);
                        h_70->Fill(pProd_2[2].Phi(),pProd_2[2].Theta());

                        theta_x_pion_m = coord_x_lambda0[1] + TMath::ATan(pProd_2[2].Px()/pProd_2[2].Pz());
                        h_9->Fill(theta_x_pion_m);
                        theta_y_pion_m = coord_x_lambda0[3] + TMath::ATan(pProd_2[2].Py()/pProd_2[2].Pz());
                        charge_pion_m = -1.0;
                        /*x*/     coord_x0_pion_m[0] = coord_x_lambda0[0];  // [m]
                        /*x'*/    coord_x0_pion_m[1] = theta_x_pion_m;      // [rad]
                        /*y*/     coord_x0_pion_m[2] = coord_x_lambda0[2];  // [m]
                        /*y'*/    coord_x0_pion_m[3] = theta_y_pion_m;      // [rad]
                        /*l*/     coord_x0_pion_m[4] = 0.0;
                        /*dp*/    coord_x0_pion_m[5] = 0.0;

                        h_60->Fill(theta_x_pion_m,theta_x_proton);

                        //--------------------//
                        // Magnets
                        //--------------------//

                        TGraph* gr_proton_x = new TGraph();
                        TGraph* gr_proton_xp = new TGraph();
                        TGraph* gr_proton_y = new TGraph();
                        TGraph* gr_proton_yp = new TGraph();
                        TString gr_proton_x_name = "gr_proton_x_";
                        TString gr_proton_xp_name = "gr_proton_xp_";
                        TString gr_proton_y_name = "gr_proton_y_";
                        TString gr_proton_yp_name = "gr_proton_yp_";
                        gr_proton_x_name += nRuns;
                        gr_proton_xp_name += nRuns;
                        gr_proton_y_name += nRuns;
                        gr_proton_yp_name += nRuns;
                        gr_proton_x->SetName(gr_proton_x_name.Data());
                        gr_proton_xp->SetName(gr_proton_xp_name.Data());
                        gr_proton_y->SetName(gr_proton_y_name.Data());
                        gr_proton_yp->SetName(gr_proton_yp_name.Data());

                        passMagnets(s_decay_2,coord_x0_proton,coord_x_proton,p_proton,charge_proton,gr_proton_x,gr_proton_y,gr_proton_xp,gr_proton_yp,false);

                        TGraph* gr_pion_m_x = new TGraph();
                        TGraph* gr_pion_m_xp = new TGraph();
                        TGraph* gr_pion_m_y = new TGraph();
                        TGraph* gr_pion_m_yp = new TGraph();
                        TString gr_pion_m_x_name = "gr_pion_m_x_";
                        TString gr_pion_m_xp_name = "gr_pion_m_xp_";
                        TString gr_pion_m_y_name = "gr_pion_m_y_";
                        TString gr_pion_m_yp_name = "gr_pion_m_yp_";
                        gr_pion_m_x_name += nRuns;
                        gr_pion_m_xp_name += nRuns;
                        gr_pion_m_y_name += nRuns;
                        gr_pion_m_yp_name += nRuns;
                        gr_pion_m_x->SetName(gr_pion_m_x_name.Data());
                        gr_pion_m_xp->SetName(gr_pion_m_xp_name.Data());
                        gr_pion_m_y->SetName(gr_pion_m_y_name.Data());
                        gr_pion_m_yp->SetName(gr_pion_m_yp_name.Data());

                        passMagnets(s_decay_2,coord_x0_pion_m,coord_x_pion_m,p_pion_m,charge_pion_m,gr_pion_m_x,gr_pion_m_y,gr_pion_m_xp,gr_pion_m_yp,false);

                        if(plot_graph)
                        {
                            gr_proton_x->Write();
                            gr_proton_xp->Write();
                            gr_proton_y->Write();
                            gr_proton_yp->Write();
                            gr_pion_m_x->Write();
                            gr_pion_m_xp->Write();
                            gr_pion_m_y->Write();
                            gr_pion_m_yp->Write();
                        }

                        if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 < Constants::_xrph_51937_ua9_pos)
                        {
                            h_12->Fill(gr_proton_x->Eval(Constants::_xrph_51937_ua9_pos),gr_proton_y->Eval(Constants::_xrph_51937_ua9_pos));
                            h_13->Fill(gr_pion_m_x->Eval(Constants::_xrph_51937_ua9_pos),gr_pion_m_y->Eval(Constants::_xrph_51937_ua9_pos));
                        }
                        if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 < Constants::_xrph_52202_ua9_pos)
                        {
                            h_16->Fill(gr_proton_x->Eval(Constants::_xrph_52202_ua9_pos),gr_proton_y->Eval(Constants::_xrph_52202_ua9_pos));
                            h_17->Fill(gr_pion_m_x->Eval(Constants::_xrph_52202_ua9_pos),gr_pion_m_y->Eval(Constants::_xrph_52202_ua9_pos));
                        }

                        if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 < Constants::_tr1_s_pos)
                        {
                            h_30->Fill(gr_proton_x->Eval(Constants::_tr1_s_pos),gr_proton_y->Eval(Constants::_tr1_s_pos));
                            h_38->Fill(gr_proton_xp->Eval(Constants::_tr1_s_pos),gr_proton_yp->Eval(Constants::_tr1_s_pos));
                            h_46->Fill(gr_proton_x->Eval(Constants::_tr1_s_pos),gr_proton_xp->Eval(Constants::_tr1_s_pos));
                            h_54->Fill(gr_proton_y->Eval(Constants::_tr1_s_pos),gr_proton_yp->Eval(Constants::_tr1_s_pos));
                            h_31->Fill(gr_pion_m_x->Eval(Constants::_tr1_s_pos),gr_pion_m_y->Eval(Constants::_tr1_s_pos));
                            h_39->Fill(gr_pion_m_xp->Eval(Constants::_tr1_s_pos),gr_pion_m_yp->Eval(Constants::_tr1_s_pos));
                            h_47->Fill(gr_pion_m_x->Eval(Constants::_tr1_s_pos),gr_pion_m_xp->Eval(Constants::_tr1_s_pos));
                            h_55->Fill(gr_pion_m_y->Eval(Constants::_tr1_s_pos),gr_pion_m_yp->Eval(Constants::_tr1_s_pos));

                            // TRACKER 1
                            if(gr_proton_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                            {
                                detected += 1;
                                h_62->Fill(gr_proton_x->Eval(Constants::_tr1_s_pos),gr_proton_y->Eval(Constants::_tr1_s_pos));

                                // TRACKER 2
                                if(gr_proton_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                                   gr_proton_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                                   gr_proton_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                                   gr_proton_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                                {
                                    detected += 1;
                                    h_65->Fill(gr_proton_x->Eval(Constants::_tr2_s_pos),gr_proton_y->Eval(Constants::_tr2_s_pos));
                                }
                            }

                            // TRACKER 1
                            if(gr_pion_m_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                            {
                                detected += 1;
                                h_63->Fill(gr_pion_m_x->Eval(Constants::_tr1_s_pos),gr_pion_m_y->Eval(Constants::_tr1_s_pos));

                                // TRACKER 2
                                if(gr_pion_m_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                                   gr_pion_m_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                                   gr_pion_m_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                                   gr_pion_m_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                                {
                                    detected += 1;
                                    h_66->Fill(gr_pion_m_x->Eval(Constants::_tr2_s_pos),gr_pion_m_y->Eval(Constants::_tr2_s_pos));
                                }
                            }

                            h_34->Fill(gr_proton_x->Eval(Constants::_tr2_s_pos),gr_proton_y->Eval(Constants::_tr2_s_pos));
                            h_42->Fill(gr_proton_xp->Eval(Constants::_tr2_s_pos),gr_proton_yp->Eval(Constants::_tr2_s_pos));
                            h_50->Fill(gr_proton_x->Eval(Constants::_tr2_s_pos),gr_proton_xp->Eval(Constants::_tr2_s_pos));
                            h_58->Fill(gr_proton_y->Eval(Constants::_tr2_s_pos),gr_proton_yp->Eval(Constants::_tr2_s_pos));
                            h_35->Fill(gr_pion_m_x->Eval(Constants::_tr2_s_pos),gr_pion_m_y->Eval(Constants::_tr2_s_pos));
                            h_43->Fill(gr_pion_m_xp->Eval(Constants::_tr2_s_pos),gr_pion_m_yp->Eval(Constants::_tr2_s_pos));
                            h_51->Fill(gr_pion_m_x->Eval(Constants::_tr2_s_pos),gr_pion_m_xp->Eval(Constants::_tr2_s_pos));
                            h_59->Fill(gr_pion_m_y->Eval(Constants::_tr2_s_pos),gr_pion_m_yp->Eval(Constants::_tr2_s_pos));
                        }

                        gr_proton_x->Delete();
                        gr_proton_xp->Delete();
                        gr_proton_y->Delete();
                        gr_proton_yp->Delete();
                        gr_pion_m_x->Delete();
                        gr_pion_m_xp->Delete();
                        gr_pion_m_y->Delete();
                        gr_pion_m_yp->Delete();
                    }

                    gr_lambda0_x->Delete();
                    gr_lambda0_xp->Delete();
                    gr_lambda0_y->Delete();
                    gr_lambda0_yp->Delete();
                    gr_pion_p_x->Delete();
                    gr_pion_p_xp->Delete();
                    gr_pion_p_y->Delete();
                    gr_pion_p_yp->Delete();

                    nRuns++;
                }

                if(detected == 6)
                {
                    nLc_reconstructed += 1;
                }
            }
        }
    }

    cout<<endl;

    cout<<"--> _quadpix_xy_dim = "<<Constants::_quadpix_xy_dim<<endl;
    cout<<"--> nLc_reconstructed = "<<nLc_reconstructed<<endl;

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
    h_15->Write();
    h_16->Write();
    h_17->Write();
    h_18->Write();
    h_19->Write();
    h_20->Write();
    h_21->Write();
    h_22->Write();
    h_23->Write();
    h_24->Write();
    h_25->Write();
    h_26->Write();
    h_27->Write();
    h_28->Write();
    h_29->Write();
    h_30->Write();
    h_31->Write();
    h_32->Write();
    h_33->Write();
    h_34->Write();
    h_35->Write();
    h_36->Write();
    h_37->Write();
    h_38->Write();
    h_39->Write();
    h_40->Write();
    h_41->Write();
    h_42->Write();
    h_43->Write();
    h_44->Write();
    h_45->Write();
    h_46->Write();
    h_47->Write();
    h_48->Write();
    h_49->Write();
    h_50->Write();
    h_51->Write();
    h_52->Write();
    h_53->Write();
    h_54->Write();
    h_55->Write();
    h_56->Write();
    h_57->Write();
    h_58->Write();
    h_59->Write();
    h_60->Write();
    h_61->Write();
    h_62->Write();
    h_63->Write();
    h_64->Write();
    h_65->Write();
    h_66->Write();
    h_67->Write();
    h_68->Write();
    h_69->Write();
    h_70->Write();
    h_71->Write();

    _file->Close();

    cout<<"--> Output file: "<<output_file_name<<endl;
}

void function_6(Double_t crystalOrientation)
{
    //------------------------------//
    const Int_t nLc     = 1e5;      // number of maximum Lc decays
    const Int_t nMom    = 45;       // number of diff. Lc mom.
    //------------------------------//
    //------------------------------//
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    TRandom3* rnd = new TRandom3(seed);
    //-----------------------------------------//
    // Magnets section
    //-----------------------------------------//

    MagClass* magnet = new MagClass();
    const Int_t mtrx_size = magnet->_mtrx_size;

    //-----------------------------------------//
    // Decay section
    //-----------------------------------------//
    DecayClass* decay = new DecayClass();
    //-----------------------------------------//

    Double_t p_lambda0 = -999, charge_lambda0 = -999, theta_x_lambda0 = -999, theta_y_lambda0 = -999;
    Double_t* coord_x0_lambda0          = new Double_t[mtrx_size];
    Double_t* coord_x_lambda0           = new Double_t[mtrx_size];

    Double_t p_pion_p = -999, charge_pion_p = -999, theta_x_pion_p = -999, theta_y_pion_p = -999;
    Double_t* coord_x0_pion_p          = new Double_t[mtrx_size];
    Double_t* coord_x_pion_p           = new Double_t[mtrx_size];

    Double_t p_proton = -999, charge_proton = -999, theta_x_proton = -999, theta_y_proton = -999;
    Double_t* coord_x0_proton          = new Double_t[mtrx_size];
    Double_t* coord_x_proton           = new Double_t[mtrx_size];

    Double_t p_pion_m = -999, charge_pion_m = -999, theta_x_pion_m = -999, theta_y_pion_m = -999;
    Double_t* coord_x0_pion_m          = new Double_t[mtrx_size];
    Double_t* coord_x_pion_m           = new Double_t[mtrx_size];

    TLorentzVector pProd_1[3];
    TLorentzVector pProd_2[3];
    Double_t mProd_1[3];
    Double_t mProd_2[3];

    Double_t s_decay_1 = 0;
    Double_t s_decay_2 = 0;

    // Mass
    mProd_1[0] = 2.28646; // Lc+    [GeV/c2]
    mProd_1[1] = 1.11568; // L0     [GeV/c2]
    mProd_1[2] = 0.13957; // Pion+  [GeV/c2]

    mProd_2[0] = 1.11568; // L0     [GeV/c2]
    mProd_2[1] = 0.93827; // Proton [GeV/c2]
    mProd_2[2] = 0.13957; // Pion-  [GeV/c2]

    cout<<endl;

    Double_t pLc, thetaYLambdaCout;
    TString output_file_name;

    for(Int_t i = 0; i <= nMom; i++)
    {
        pLc = 50.0 + i*10.0;// [GeV/c]
        cout<<endl<<"--> pLc = "<<pLc<<" [GeV/c]"<<endl;
        cout<<"--> crystalOrientation = "<<(Int_t)(crystalOrientation*1e3)<<" [mrad]"<<endl;
        output_file_name = "./output/accelerator_";
        output_file_name += (Int_t)(crystalOrientation*1e3);
        output_file_name += "mrad_";
        output_file_name += (Int_t)pLc;
        output_file_name += "GeV.root";

        TFile* file = new TFile(output_file_name.Data(),"RECREATE");

        TH1D* h_1 = new TH1D("h_1","L_{c} momentum [GeV/c]",550,-50.0,500.0);

        TH1D* h_2 = new TH1D("h_2","L_{0} momentum [GeV/c]",550,-50.0,500.0);
        TH1D* h_3 = new TH1D("h_3","pion+ momentum [GeV/c]",550,-50.0,500.0);
        TH1D* h_4 = new TH1D("h_4","proton momentum [GeV/c]",550,-50.0,500.0);
        TH1D* h_5 = new TH1D("h_5","pion- momentum [GeV/c]",550,-50.0,500.0);

        TH1D* h_6 = new TH1D("h_6","L_{0} thetaX",400000,-4.0,4.0);
        TH1D* h_7 = new TH1D("h_7","pion+ thetaX",400000,-4.0,4.0);
        TH1D* h_8 = new TH1D("h_8","proton thetaX",400000,-4.0,4.0);
        TH1D* h_9 = new TH1D("h_9","pion- thetaX",400000,-4.0,4.0);

        TH1D* h_10 = new TH1D("h_10","L_{c} thetaY",400000,-4.0,4.0);
        TH1D* h_11 = new TH1D("h_11","L_{0} thetaY",400000,-4.0,4.0);
        TH1D* h_12 = new TH1D("h_12","pion+ thetaY",400000,-4.0,4.0);
        TH1D* h_13 = new TH1D("h_13","proton thetaY",400000,-4.0,4.0);
        TH1D* h_14 = new TH1D("h_14","pion- thetaY",400000,-4.0,4.0);

        TH1D* h_18 = new TH1D("h_18","L_{0} decay length [m]",11000,-10,100);
        TH1D* h_19 = new TH1D("h_19","L_{0} decay length (z-proj) [m]",11000,-10,100);

        TH2D* h_24 = new TH2D("h_24","L_{0} momentum [GeV/c] vs L_{c} momentum [GeV/c]",550,-50.0,500.0,550,-50.0,500.0);
        TH2D* h_25 = new TH2D("h_25","pion+ momentum [GeV/c] vs L_{c} momentum [GeV/c]",550,-50.0,500.0,550,-50.0,500.0);
        TH2D* h_26 = new TH2D("h_26","proton momentum [GeV/c] vs L_{c} momentum [GeV/c]",550,-50.0,500.0,550,-50.0,500.0);
        TH2D* h_27 = new TH2D("h_27","pion- momentum [GeV/c] vs L_{c} momentum [GeV/c]",550,-50.0,500.0,550,-50.0,500.0);

        TH2D* h_29 = new TH2D("h_29","XY pion+ at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_30 = new TH2D("h_30","XY proton at TR1",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_31 = new TH2D("h_31","XY pion- at TR1",4000,-0.2,0.2,4000,-0.2,0.2);

        TH2D* h_33 = new TH2D("h_33","XY pion+ at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_34 = new TH2D("h_34","XY proton at TR2",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_35 = new TH2D("h_35","XY pion- at TR2",4000,-0.2,0.2,4000,-0.2,0.2);

        TH2D* h_61 = new TH2D("h_61","XY pion+ at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_62 = new TH2D("h_62","XY proton at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_63 = new TH2D("h_63","XY pion- at TR1 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);

        TH2D* h_64 = new TH2D("h_64","XY pion+ at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_65 = new TH2D("h_65","XY proton at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);
        TH2D* h_66 = new TH2D("h_66","XY pion- at TR2 (detected)",4000,-0.2,0.2,4000,-0.2,0.2);

        for(Int_t j = 0; j < nLc; j++)
        {
            if(j%1000 == 0)
            {
                printf("\r--> Working: %3.1f",100.0*j/nLc);
                fflush(stdout);
            }

            Double_t eLc = TMath::Sqrt(pLc*pLc + mProd_1[0]*mProd_1[0]); // [GeV/c]
            pProd_1[0].SetPxPyPzE(0,0,pLc,eLc);
            thetaYLambdaCout = gRandom->Gaus(0.0,4.0e-3);

            h_1->Fill(pLc);
            h_10->Fill(thetaYLambdaCout);

            if(decay->twoBody(pProd_1,mProd_1))
            {
                //--------------------//
                // Lambda0
                //--------------------//

                p_lambda0 = pProd_1[1].P();
                h_2->Fill(p_lambda0);
                h_24->Fill(pLc,p_lambda0);

                theta_x_lambda0 =
                        Constants::_beamAngleInitialAtCryPosition +
                        Constants::_crystalAngle +
                        crystalOrientation +
                        TMath::ATan(pProd_1[1].Px()/pProd_1[1].Pz());
                h_6->Fill(theta_x_lambda0);
                theta_y_lambda0 = thetaYLambdaCout + TMath::ATan(pProd_1[1].Py()/pProd_1[1].Pz());
                h_11->Fill(theta_y_lambda0);
                charge_lambda0 = 0.0;
                /*x*/     coord_x0_lambda0[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                /*x'*/    coord_x0_lambda0[1] = theta_x_lambda0;// [rad]
                /*y*/     coord_x0_lambda0[2] = 0.0;              // [m]
                /*y'*/    coord_x0_lambda0[3] = theta_y_lambda0;  // [rad]
                /*l*/     coord_x0_lambda0[4] = 0.0;
                /*dp*/    coord_x0_lambda0[5] = 0.0;

                // RANDOM //
                Double_t gamma_L0   = getGamma(mProd_2[0],p_lambda0);
                Double_t beta_L0    = getBeta(gamma_L0);
                Double_t prob       = rnd->Uniform(0.0,1.0);
                Double_t cosTheta = pProd_1[1].CosTheta();

                if(prob == 0)
                {
                    while( 1 )
                    {
                        prob = rnd->Uniform(0.0,1.0);
                        if(prob > 0) break;
                    }
                }
                s_decay_2 = Constants::_ctau_Lambda0*beta_L0*gamma_L0*TMath::Log(1.0/prob); // [m]
                h_18->Fill(s_decay_2);
                s_decay_2 = s_decay_2*cosTheta;
                h_19->Fill(s_decay_2);

                //--------------------//
                // Pion+
                //--------------------//

                p_pion_p = pProd_1[2].P();
                h_3->Fill(p_pion_p);
                h_25->Fill(pLc,p_pion_p);

                theta_x_pion_p =
                        Constants::_beamAngleInitialAtCryPosition +
                        Constants::_crystalAngle +
                        crystalOrientation +
                        TMath::ATan(pProd_1[2].Px()/pProd_1[2].Pz());
                h_7->Fill(theta_x_pion_p);
                theta_y_pion_p = thetaYLambdaCout + TMath::ATan(pProd_1[2].Py()/pProd_1[2].Pz());
                h_12->Fill(theta_y_pion_p);
                charge_pion_p = 1.0;
                /*x*/     coord_x0_pion_p[0] = Constants::_beamPositionInitialAtCryPosition;    // [m]
                /*x'*/    coord_x0_pion_p[1] = theta_x_pion_p;// [rad]
                /*y*/     coord_x0_pion_p[2] = 0.0;           // [m]
                /*y'*/    coord_x0_pion_p[3] = theta_y_pion_p;  // [rad]
                /*l*/     coord_x0_pion_p[4] = 0.0;
                /*dp*/    coord_x0_pion_p[5] = 0.0;

                //--------------------//
                // Magnets
                //--------------------//

                TGraph* gr_lambda0_x = new TGraph();
                TGraph* gr_lambda0_xp = new TGraph();
                TGraph* gr_lambda0_y = new TGraph();
                TGraph* gr_lambda0_yp = new TGraph();
                TString gr_lambda0_x_name = "gr_lambda0_x_";
                TString gr_lambda0_xp_name = "gr_lambda0_xp_";
                TString gr_lambda0_y_name = "gr_lambda0_y_";
                TString gr_lambda0_yp_name = "gr_lambda0_yp_";
                gr_lambda0_x->SetName(gr_lambda0_x_name.Data());
                gr_lambda0_xp->SetName(gr_lambda0_xp_name.Data());
                gr_lambda0_y->SetName(gr_lambda0_y_name.Data());
                gr_lambda0_yp->SetName(gr_lambda0_yp_name.Data());

                if(charge_lambda0 != 0)
                {
                    assert(0);
                }

                // INITIAL
                gr_lambda0_x->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[0]);
                gr_lambda0_xp->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[1]);
                gr_lambda0_y->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[2]);
                gr_lambda0_yp->SetPoint(0,Constants::_cry3_51799_ua9_pos+s_decay_1,coord_x0_lambda0[3]);

                // DRIFT & FINAL
                magnet->GetNewCoordDrift(s_decay_2,Constants::_order_transport_matrix,coord_x0_lambda0,coord_x_lambda0,Constants::_aph,Constants::_apv);
                gr_lambda0_x->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[0]);
                gr_lambda0_xp->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[1]);
                gr_lambda0_y->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[2]);
                gr_lambda0_yp->SetPoint(1,Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2,coord_x_lambda0[3]);

                TGraph* gr_pion_p_x = new TGraph();
                TGraph* gr_pion_p_xp = new TGraph();
                TGraph* gr_pion_p_y = new TGraph();
                TGraph* gr_pion_p_yp = new TGraph();
                TString gr_pion_p_x_name = "gr_pion_p_x_";
                TString gr_pion_p_xp_name = "gr_pion_p_xp_";
                TString gr_pion_p_y_name = "gr_pion_p_y_";
                TString gr_pion_p_yp_name = "gr_pion_p_yp_";
                gr_pion_p_x->SetName(gr_pion_p_x_name.Data());
                gr_pion_p_xp->SetName(gr_pion_p_xp_name.Data());
                gr_pion_p_y->SetName(gr_pion_p_y_name.Data());
                gr_pion_p_yp->SetName(gr_pion_p_yp_name.Data());

                passMagnets(s_decay_1,coord_x0_pion_p,coord_x_pion_p,p_pion_p,charge_pion_p,gr_pion_p_x,gr_pion_p_y,gr_pion_p_xp,gr_pion_p_yp,false);

                h_29->Fill(gr_pion_p_x->Eval(Constants::_tr1_s_pos),gr_pion_p_y->Eval(Constants::_tr1_s_pos));
                h_33->Fill(gr_pion_p_x->Eval(Constants::_tr2_s_pos),gr_pion_p_y->Eval(Constants::_tr2_s_pos));

                // TRACKER 1
                if(gr_pion_p_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                   gr_pion_p_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                   gr_pion_p_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                   gr_pion_p_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                {
                    h_61->Fill(gr_pion_p_x->Eval(Constants::_tr1_s_pos),gr_pion_p_y->Eval(Constants::_tr1_s_pos));

                    // TRACKER 2
                    if(gr_pion_p_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                       gr_pion_p_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                    {
                        h_64->Fill(gr_pion_p_x->Eval(Constants::_tr2_s_pos),gr_pion_p_y->Eval(Constants::_tr2_s_pos));
                    }
                }

                pProd_2[0].SetPxPyPzE(0.0,0.0,pProd_1[1].P(),pProd_1[1].E());

                if(decay->twoBody(pProd_2,mProd_2))
                {
                    //--------------------//
                    // Proton
                    //--------------------//

                    p_proton = pProd_2[1].P();
                    h_4->Fill(p_proton);
                    h_26->Fill(pLc,p_proton);

                    theta_x_proton = coord_x_lambda0[1] + TMath::ATan(pProd_2[1].Px()/pProd_2[1].Pz());
                    h_8->Fill(theta_x_proton);
                    theta_y_proton = coord_x_lambda0[3] + TMath::ATan(pProd_2[1].Py()/pProd_2[1].Pz());
                    h_13->Fill(theta_y_proton);
                    charge_proton = 1.0;
                    /*x*/     coord_x0_proton[0] = coord_x_lambda0[0];  // [m]
                    /*x'*/    coord_x0_proton[1] = theta_x_proton;      // [rad]
                    /*y*/     coord_x0_proton[2] = coord_x_lambda0[2];  // [m]
                    /*y'*/    coord_x0_proton[3] = theta_y_proton;      // [rad]
                    /*l*/     coord_x0_proton[4] = 0.0;
                    /*dp*/    coord_x0_proton[5] = 0.0;

                    //--------------------//
                    // Pion-
                    //--------------------//

                    p_pion_m = pProd_2[2].P();
                    h_5->Fill(p_pion_m);
                    h_27->Fill(pLc,p_pion_m);

                    theta_x_pion_m = coord_x_lambda0[1] + TMath::ATan(pProd_2[2].Px()/pProd_2[2].Pz());
                    h_9->Fill(theta_x_pion_m);
                    theta_y_pion_m = coord_x_lambda0[3] + TMath::ATan(pProd_2[2].Py()/pProd_2[2].Pz());
                    h_14->Fill(theta_y_pion_m);
                    charge_pion_m = -1.0;
                    /*x*/     coord_x0_pion_m[0] = coord_x_lambda0[0];  // [m]
                    /*x'*/    coord_x0_pion_m[1] = theta_x_pion_m;      // [rad]
                    /*y*/     coord_x0_pion_m[2] = coord_x_lambda0[2];  // [m]
                    /*y'*/    coord_x0_pion_m[3] = theta_y_pion_m;      // [rad]
                    /*l*/     coord_x0_pion_m[4] = 0.0;
                    /*dp*/    coord_x0_pion_m[5] = 0.0;

                    //--------------------//
                    // Magnets
                    //--------------------//

                    TGraph* gr_proton_x = new TGraph();
                    TGraph* gr_proton_xp = new TGraph();
                    TGraph* gr_proton_y = new TGraph();
                    TGraph* gr_proton_yp = new TGraph();
                    TString gr_proton_x_name = "gr_proton_x_";
                    TString gr_proton_xp_name = "gr_proton_xp_";
                    TString gr_proton_y_name = "gr_proton_y_";
                    TString gr_proton_yp_name = "gr_proton_yp_";
                    gr_proton_x->SetName(gr_proton_x_name.Data());
                    gr_proton_xp->SetName(gr_proton_xp_name.Data());
                    gr_proton_y->SetName(gr_proton_y_name.Data());
                    gr_proton_yp->SetName(gr_proton_yp_name.Data());

                    passMagnets(s_decay_2,coord_x0_proton,coord_x_proton,p_proton,charge_proton,gr_proton_x,gr_proton_y,gr_proton_xp,gr_proton_yp,false);

                    TGraph* gr_pion_m_x = new TGraph();
                    TGraph* gr_pion_m_xp = new TGraph();
                    TGraph* gr_pion_m_y = new TGraph();
                    TGraph* gr_pion_m_yp = new TGraph();
                    TString gr_pion_m_x_name = "gr_pion_m_x_";
                    TString gr_pion_m_xp_name = "gr_pion_m_xp_";
                    TString gr_pion_m_y_name = "gr_pion_m_y_";
                    TString gr_pion_m_yp_name = "gr_pion_m_yp_";
                    gr_pion_m_x->SetName(gr_pion_m_x_name.Data());
                    gr_pion_m_xp->SetName(gr_pion_m_xp_name.Data());
                    gr_pion_m_y->SetName(gr_pion_m_y_name.Data());
                    gr_pion_m_yp->SetName(gr_pion_m_yp_name.Data());

                    passMagnets(s_decay_2,coord_x0_pion_m,coord_x_pion_m,p_pion_m,charge_pion_m,gr_pion_m_x,gr_pion_m_y,gr_pion_m_xp,gr_pion_m_yp,false);

                    if(Constants::_cry3_51799_ua9_pos+s_decay_1+s_decay_2 < Constants::_tr1_s_pos)
                    {
                        h_30->Fill(gr_proton_x->Eval(Constants::_tr1_s_pos),gr_proton_y->Eval(Constants::_tr1_s_pos));
                        h_31->Fill(gr_pion_m_x->Eval(Constants::_tr1_s_pos),gr_pion_m_y->Eval(Constants::_tr1_s_pos));

                        // TRACKER 1
                        if(gr_proton_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_proton_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                           gr_proton_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_proton_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                        {
                            h_62->Fill(gr_proton_x->Eval(Constants::_tr1_s_pos),gr_proton_y->Eval(Constants::_tr1_s_pos));

                            // TRACKER 2
                            if(gr_proton_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_proton_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                            {
                                h_65->Fill(gr_proton_x->Eval(Constants::_tr2_s_pos),gr_proton_y->Eval(Constants::_tr2_s_pos));
                            }
                        }

                        // TRACKER 1
                        if(gr_pion_m_x->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_m_x->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_m_y->Eval(Constants::_tr1_s_pos) <= Constants::_tr1_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                           gr_pion_m_y->Eval(Constants::_tr1_s_pos) >= Constants::_tr1_y_pos - Constants::_quadpix_xy_dim/2.0)
                        {
                            h_63->Fill(gr_pion_m_x->Eval(Constants::_tr1_s_pos),gr_pion_m_y->Eval(Constants::_tr1_s_pos));

                            // TRACKER 2
                            if(gr_pion_m_x->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_x_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_x->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_x_pos - Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_y->Eval(Constants::_tr2_s_pos) <= Constants::_tr2_y_pos + Constants::_quadpix_xy_dim/2.0 &&
                               gr_pion_m_y->Eval(Constants::_tr2_s_pos) >= Constants::_tr2_y_pos - Constants::_quadpix_xy_dim/2.0)
                            {
                                h_66->Fill(gr_pion_m_x->Eval(Constants::_tr2_s_pos),gr_pion_m_y->Eval(Constants::_tr2_s_pos));
                            }
                        }

                        h_34->Fill(gr_proton_x->Eval(Constants::_tr2_s_pos),gr_proton_y->Eval(Constants::_tr2_s_pos));
                        h_35->Fill(gr_pion_m_x->Eval(Constants::_tr2_s_pos),gr_pion_m_y->Eval(Constants::_tr2_s_pos));
                    }

                    gr_proton_x->Delete();
                    gr_proton_xp->Delete();
                    gr_proton_y->Delete();
                    gr_proton_yp->Delete();
                    gr_pion_m_x->Delete();
                    gr_pion_m_xp->Delete();
                    gr_pion_m_y->Delete();
                    gr_pion_m_yp->Delete();
                }

                gr_lambda0_x->Delete();
                gr_lambda0_xp->Delete();
                gr_lambda0_y->Delete();
                gr_lambda0_yp->Delete();
                gr_pion_p_x->Delete();
                gr_pion_p_xp->Delete();
                gr_pion_p_y->Delete();
                gr_pion_p_yp->Delete();
            }
        }

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
        h_18->Write();
        h_19->Write();
        h_24->Write();
        h_25->Write();
        h_26->Write();
        h_27->Write();
        h_29->Write();
        h_30->Write();
        h_31->Write();
        h_33->Write();
        h_34->Write();
        h_35->Write();
        h_61->Write();
        h_62->Write();
        h_63->Write();
        h_64->Write();
        h_65->Write();
        h_66->Write();

        file->Close();

        cout<<endl<<"--> Output file: "<<output_file_name<<endl;
    }

    cout<<endl;

    cout<<"--> _quadpix_xy_dim = "<<Constants::_quadpix_xy_dim<<endl;
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

void passMagnets(Double_t s, Double_t* coord_x0, Double_t* coord_x, Double_t p, Double_t Charge_C, TGraph *gr_x, TGraph *gr_y, TGraph *gr_xp, TGraph *gr_yp, Bool_t print)
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
    Bool_t s_status = kFALSE;

    if(!s_status && s <= 0)
    {
        // INITIAL
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[0]);
        gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[1]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[2]);
        gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos,coord_x0[3]);
        gr_x_ipoint++;
        gr_y_ipoint++;

        s_status = kTRUE;
    }

    if(Constants::_switch_magnets && Charge_C != 0)
    {
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_q1_51810_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL1 = Constants::_q1_51810_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT1-S
            if(!magnet->GetNewCoordDrift(DRIFTL1,order,coord_x0,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT1
            if(!magnet->GetNewCoordDrift(DRIFTL1,order,coord_x0,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
        if(s_status)
        {
            // QUAD1: FOC in X, DEFOC in Y
            if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ1,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_q1_51810_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q1_51810_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_q1_51810_pos,coord_temp_2_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_q2_51910_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL2 = Constants::_q2_51910_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT2-S
            if(!magnet->GetNewCoordDrift(DRIFTL2,order,coord_x0,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT2
            if(!magnet->GetNewCoordDrift(DRIFTL2,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
        if(s_status)
        {
            // QUAD2: FOC in Y, DEFOC in X
            if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ2,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_q2_51910_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q2_51910_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_q2_51910_pos,coord_temp_2_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_xrph_51937_ua9_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL3 = Constants::_xrph_51937_ua9_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT3-S
            if(!magnet->GetNewCoordDrift(DRIFTL3,order,coord_x0,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT3
            if(!magnet->GetNewCoordDrift(DRIFTL3,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_xrph_51937_ua9_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_xrph_51937_ua9_pos,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_lsf_52005_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL4 = Constants::_lsf_52005_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT4-S
            if(!magnet->GetNewCoordDrift(DRIFTL4,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT4
            if(!magnet->GetNewCoordDrift(DRIFTL4,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
        if(s_status)
        {
            // SEXT1
            if(!magnet->GetNewCoordSextupole(Charge_C*KS1,p,P0,SL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_lsf_52005_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_lsf_52005_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_lsf_52005_pos,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_q3_52010_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL5 = Constants::_q3_52010_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT5-S
            if(!magnet->GetNewCoordDrift(DRIFTL5,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT5
            if(!magnet->GetNewCoordDrift(DRIFTL5,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
        if(s_status)
        {
            // QUAD3: FOC in X, DEFOC in Y
            if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ3,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_q3_52010_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q3_52010_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_q3_52010_pos,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_q3_52010_pos + DRIFTL6 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL6 = Constants::_q3_52010_pos + DRIFTL6 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT6-S
            if(!magnet->GetNewCoordDrift(DRIFTL6,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT6
            if(!magnet->GetNewCoordDrift(DRIFTL6,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
        if(s_status)
        {
            // DIPOLE1: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mba_52030_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52030_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mba_52030_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_mba_52030_pos + DL/2 + DRIFTL7 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL7 = Constants::_mba_52030_pos + DL/2 + DRIFTL7 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT7
            if(!magnet->GetNewCoordDrift(DRIFTL7,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT7
            if(!magnet->GetNewCoordDrift(DRIFTL7,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mba_52030_pos+DL/2+DRIFTL7);
        }
        if(s_status)
        {
            // DIPOLE2: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mba_52050_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52050_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mba_52050_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_mba_52050_pos + DL/2 + DRIFTL8 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL8 = Constants::_mba_52050_pos + DL/2 + DRIFTL8 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT8
            if(!magnet->GetNewCoordDrift(DRIFTL8,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT8
            if(!magnet->GetNewCoordDrift(DRIFTL8,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mba_52050_pos+DL/2+DRIFTL8);
        }
        if(s_status)
        {
            // DIPOLE3: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mbb_52070_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52070_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mbb_52070_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_mbb_52070_pos + DL/2 + DRIFTL9 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL9 = Constants::_mbb_52070_pos + DL/2 + DRIFTL9 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT9-S
            if(!magnet->GetNewCoordDrift(DRIFTL9,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT9
            if(!magnet->GetNewCoordDrift(DRIFTL9,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52070_pos+DL/2+DRIFTL9);
        }
        if(s_status)
        {
            // DIPOLE4: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mbb_52090_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52090_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mbb_52090_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_mbb_52090_pos + DL/2 + DRIFTL10 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL10 = Constants::_mbb_52090_pos+DL/2+DRIFTL10 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT10-S
            if(!magnet->GetNewCoordDrift(DRIFTL10,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT10
            if(!magnet->GetNewCoordDrift(DRIFTL10,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52090_pos+DL/2+DRIFTL10);
        }
        if(s_status)
        {
            // QUAD4: FOC in Y, DEFOC in X
            if(!magnet->GetNewCoordQuadrupole(Charge_C*KQ4,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_q4_52110_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q4_52110_pos);
            gr_x->SetPoint(gr_x_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_q4_52110_pos,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_q4_52110_pos + DRIFTL11 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL11 = Constants::_q4_52110_pos + DRIFTL11 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT11-S
            if(!magnet->GetNewCoordDrift(DRIFTL11,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT11
            if(!magnet->GetNewCoordDrift(DRIFTL11,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q4_52110_pos+DRIFTL11);
        }
        if(s_status)
        {
            // DIPOLE5: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mbb_52130_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52130_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mbb_52130_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_mbb_52130_pos + DL/2 + DRIFTL12 - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL12 = Constants::_mbb_52130_pos + DL/2 + DRIFTL12 - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT12-S
            if(!magnet->GetNewCoordDrift(DRIFTL12,order,coord_x0,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT12
            if(!magnet->GetNewCoordDrift(DRIFTL12,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_mbb_52130_pos+DL/2+DRIFTL12);
        }
        if(s_status)
        {
            // DIPOLE6: BEND & DEFOC in X, FOC in Y
            if(!magnet->GetNewCoordDipole(P0,p,ADL,Charge_C,DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_mbb_52150_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52150_pos+DL/2);
            gr_x->SetPoint(gr_x_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_mbb_52150_pos+DL/2,coord_temp_1_x[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
        if(!s_status && s < Constants::_xrph_52202_ua9_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL13 = Constants::_xrph_52202_ua9_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT13-S
            if(!magnet->GetNewCoordDrift(DRIFTL13,order,coord_x0,coord_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_xrph_52202_ua9_pos<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT13
            if(!magnet->GetNewCoordDrift(DRIFTL13,order,coord_temp_1_x,coord_x,APH,APV))
            {if(print) cout<<"## At : "<<Constants::_xrph_52202_ua9_pos<<" [m] ##"<<endl; /*continue;*/}
        }
        //-----------------------------------------------------------------------------------------------------------------------------//
    }
    else
    {
        if(!s_status && s < Constants::_xrph_52202_ua9_pos - Constants::_cry3_51799_ua9_pos)
        {
            DRIFTL0 = Constants::_xrph_52202_ua9_pos - (Constants::_cry3_51799_ua9_pos + s);
            // DRIFT0-S
            if(!magnet->GetNewCoordDrift(DRIFTL0,order,coord_x0,coord_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
            if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],Constants::_cry3_51799_ua9_pos+s);
            gr_x->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[0]);
            gr_xp->SetPoint(gr_x_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[1]);
            gr_y->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[2]);
            gr_yp->SetPoint(gr_y_ipoint,Constants::_cry3_51799_ua9_pos+s,coord_x0[3]);
            gr_x_ipoint++;
            gr_y_ipoint++;

            s_status = kTRUE;
        }
        else if(s_status)
        {
            // DRIFT0
            if(!magnet->GetNewCoordDrift(DRIFTL0,order,coord_x0,coord_x,APH,APV))
            {if(print) cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; /*continue;*/}
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------//
    if(s_status)
    {
        // FINAL
        if(print) printf("%10.15f , %10.15f , %10.15f\n",coord_x[0],coord_x[1],Constants::_xrph_52202_ua9_pos);
        gr_x->SetPoint(gr_x_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[0]);
        gr_xp->SetPoint(gr_x_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[1]);
        gr_y->SetPoint(gr_y_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[2]);
        gr_yp->SetPoint(gr_y_ipoint,Constants::_xrph_52202_ua9_pos,coord_x[3]);
    }
}

Double_t getGamma(Double_t mass, Double_t pmag)
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

Double_t getBeta(Double_t gamma)
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

