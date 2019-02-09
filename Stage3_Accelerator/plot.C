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

int plot()
{
    //----------------------------------------------------------------//
    // ACCEL 270 GeV/c proton 5mrad
    //----------------------------------------------------------------//
    Double_t s_accl[]={5180.829499999999825,    5185.170399999999972,    5217.168099999999868,    5223.927249999999731,
    5246.616299999999683,    5249.165799999999763,    5251.068299999999908,    5257.328300000000127,    5257.728299999999763,
    5263.988299999999981,    5264.378300000000309,    5270.638300000000527,    5271.018299999999726,    5277.278299999999945,
    5281.163499999999658,    5281.163499999999658,    5283.055999999999585,    5289.315999999999804,    5289.695999999999913,
    5295.956000000000131,    5309.700034999999843};

    Double_t x_accl[]={-0.022728400000000,    -0.089680588035702,    -0.480689728726413,    -0.679365876581929,
    -1.346280067345753,    -1.421219200070676,    -1.380550322565201,    -1.246735534695320,    -1.238184916560144,
    -1.104370128690272,    -1.096033276008461,    -0.962218488138598,    -0.954095400910191,    -0.820280613040335,
    -0.737228459093313,    -0.737228459093313,    -0.746620245387641,    -0.777685785116631,    -0.779571586327646,
    -0.810637126056642,    -0.878843751970417};

    TGraph* gr_accl = new TGraph(21,s_accl,x_accl);
    gr_accl->SetName("gr_accl");
    gr_accl->SetLineColor(kBlack);
    gr_accl->SetLineWidth(2);
    //----------------------------------------------------------------//
    TFile* _file0 = TFile::Open("accelerator.root");

    const int nGraph = 10;

//    const double _cry3_51799_ua9_pos = 5180.8295;   // [m]
//    const double _cry4_51799_ua9_pos = 5181.3245;   // [m]
//    const double _q1_51810_pos = 5185.1704;         // [m]
//    const double _q2_51910_pos = 5217.1681;         // [m]
    const double _xrph_51937_ua9_pos = 5223.92725;  // [m]
//    const double _lsf_52005_pos = 5246.6163;        // [m]
//    const double _q3_52010_pos = 5249.1658;         // [m]
//    const double _mba_52030_pos = 5254.1983;        // [m]
//    const double _mba_52050_pos = 5260.8583;        // [m]
//    const double _mbb_52070_pos = 5267.5083;        // [m]
//    const double _mbb_52090_pos = 5274.1483;        // [m]
//    const double _q4_52110_pos = 5281.1635;         // [m]
//    const double _mbb_52130_pos = 5286.1860;        // [m]
//    const double _mbb_52150_pos = 5292.8260;        // [m]
    const double _xrph_52202_ua9_pos = 5309.700035; // [m]

    TGraph* gr_proton_x[nGraph];
    TGraph* gr_proton_y[nGraph];
    TGraph* gr_kaon_x[nGraph];
    TGraph* gr_kaon_y[nGraph];
    TGraph* gr_pion_x[nGraph];
    TGraph* gr_pion_y[nGraph];

    TMultiGraph* mg_x = new TMultiGraph();
    TMultiGraph* mg_y = new TMultiGraph();
    mg_x->SetName("mg_x");
    mg_y->SetName("mg_y");

    for(Int_t i = 0; i < nGraph; i++)
    {
        TString gr_proton_x_name = "gr_proton_x_";
        TString gr_proton_y_name = "gr_proton_y_";
        gr_proton_x_name += i;
        gr_proton_y_name += i;
        gr_proton_x[i] = (TGraph*)_file0->Get(gr_proton_x_name.Data());
        gr_proton_y[i] = (TGraph*)_file0->Get(gr_proton_y_name.Data());
        gr_proton_x[i]->SetLineColor(kRed);
        gr_proton_y[i]->SetLineColor(kRed);

        TString gr_kaon_x_name = "gr_kaon_x_";
        TString gr_kaon_y_name = "gr_kaon_y_";
        gr_kaon_x_name += i;
        gr_kaon_y_name += i;
        gr_kaon_x[i] = (TGraph*)_file0->Get(gr_kaon_x_name.Data());
        gr_kaon_y[i] = (TGraph*)_file0->Get(gr_kaon_y_name.Data());
        gr_kaon_x[i]->SetLineColor(kBlue);
        gr_kaon_y[i]->SetLineColor(kBlue);

        TString gr_pion_x_name = "gr_pion_x_";
        TString gr_pion_y_name = "gr_pion_y_";
        gr_pion_x_name += i;
        gr_pion_y_name += i;
        gr_pion_x[i] = (TGraph*)_file0->Get(gr_pion_x_name.Data());
        gr_pion_y[i] = (TGraph*)_file0->Get(gr_pion_y_name.Data());
        gr_pion_x[i]->SetLineColor(kGreen);
        gr_pion_y[i]->SetLineColor(kGreen);

        mg_x->Add(gr_proton_x[i]);
        mg_x->Add(gr_kaon_x[i]);
        mg_x->Add(gr_pion_x[i]);

        mg_y->Add(gr_proton_y[i]);
        mg_y->Add(gr_kaon_y[i]);
        mg_y->Add(gr_pion_y[i]);
    }

    mg_x->Add(gr_accl);

    TLine* line_xrph_51937_ua9 = new TLine(_xrph_51937_ua9_pos,-1,_xrph_51937_ua9_pos,1);
    TLine* line_xrph_52202_ua9 = new TLine(_xrph_52202_ua9_pos,-1,_xrph_52202_ua9_pos,1);

    TCanvas* c_1_x = new TCanvas("c_1_x","c_1_x");
    c_1_x->cd();
    mg_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_1_y = new TCanvas("c_1_y","c_1_y");
    c_1_y->cd();
    mg_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TH2D* h_proton_1 = (TH2D*)_file0->Get("h_3");
    TH2D* h_kaon_1 = (TH2D*)_file0->Get("h_4");
    TH2D* h_pion_1 = (TH2D*)_file0->Get("h_5");

    TH2D* h_proton_2 = (TH2D*)_file0->Get("h_6");
    TH2D* h_kaon_2 = (TH2D*)_file0->Get("h_7");
    TH2D* h_pion_2 = (TH2D*)_file0->Get("h_8");

    Int_t first_bin = h_proton_1->GetXaxis()->FindBin(-4);
    Int_t last_bin = h_proton_1->GetXaxis()->FindBin(4);

    h_proton_1->GetXaxis()->SetRange(first_bin,last_bin);
    h_kaon_1->GetXaxis()->SetRange(first_bin,last_bin);
    h_pion_1->GetXaxis()->SetRange(first_bin,last_bin);
    h_proton_2->GetXaxis()->SetRange(first_bin,last_bin);
    h_kaon_2->GetXaxis()->SetRange(first_bin,last_bin);
    h_pion_2->GetXaxis()->SetRange(first_bin,last_bin);

    h_proton_1->SetMarkerStyle(7);
    h_kaon_1->SetMarkerStyle(7);
    h_pion_1->SetMarkerStyle(7);
    h_proton_2->SetMarkerStyle(7);
    h_kaon_2->SetMarkerStyle(7);
    h_pion_2->SetMarkerStyle(7);

    TCanvas* c_2 = new TCanvas("c_2","c_2");
    c_2->Divide(3,2);

    c_2->cd(1);
    h_proton_1->Draw();
    c_2->cd(2);
    h_kaon_1->Draw();
    c_2->cd(3);
    h_pion_1->Draw();

    c_2->cd(4);
    h_proton_2->Draw();
    c_2->cd(5);
    h_kaon_2->Draw();
    c_2->cd(6);
    h_pion_2->Draw();
}
