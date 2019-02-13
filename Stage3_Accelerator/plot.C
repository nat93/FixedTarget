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

void function_1();
void function_2();

int plot()
{
    cout<<endl;
    cout<<"--> function_1() -- threebody decay Lc -> p K- pi+"<<endl;
    cout<<"--> function_2() -- twobody decay Lc -> L0 pi+ -> (p pi-) pi+"<<endl;
    return 0;
}

void function_1()
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

void function_2()
{
    //----------------------------------------------------------------//
    // ACCEL 270 GeV/c proton 15 mrad
    //----------------------------------------------------------------//
    Double_t s_accl[]={5180.829499999999825,5185.170399999999972,5217.168099999999868,5223.927249999999731,5246.616299999999683,
                       5249.165799999999763,5257.328300000000127,5257.728299999999763,5263.988299999999981,5264.378300000000309,
                       5270.638300000000527,5271.018299999999726,5277.278299999999945,5281.163499999999658,5281.163499999999658,
                       5283.055999999999585,5289.315999999999804,5289.695999999999913,5295.956000000000131,5309.700034999999843};

    Double_t x_accl[]={-0.022728400000000,-0.090166508381702,-0.484202038598296,-0.684365651023316,-1.356272948611327,
                       -1.431773141664835,-1.256008376565142,-1.247394977799398,-1.112597690579324,-1.104199626782710,
                       -0.969402339562644,-0.961219610735199,-0.826422323515140,-0.742760381303421,-0.742760381303421,
                       -0.752229171627299,-0.783549420275793,-0.785450683326295,-0.816770931974795,-0.885536789633491};

    TGraph* gr_accl = new TGraph(20,s_accl,x_accl);
    gr_accl->SetName("gr_accl");
    gr_accl->SetLineColor(kCyan);
    gr_accl->SetLineWidth(1);
    //----------------------------------------------------------------//
    // MADX 270 GeV/c proton 15 mrad
    //----------------------------------------------------------------//
    Double_t s_madx[]={5180.829499999999825,                       5181.077000000000226,                       5181.324499999999716,                       5181.940700000000106,
                       5182.556899999999587,                       5182.589399999999841,                       5182.621900000000096,                       5182.864899999999579,
                       5183.107899999999972,                       5183.251650000000154,                       5183.395400000000336,                       5184.282900000000154,
                       5185.170399999999972,                       5186.193400000000111,                       5187.216400000000249,                       5187.628649999999652,
                       5188.040899999999965,                       5188.367400000000089,                       5188.693900000000212,                       5188.693905000000086,
                       5188.693909999999960,                       5189.627029999999650,                       5190.560150000000249,                       5190.560155000000123,
                       5190.560159999999996,                       5190.628655000000435,                       5190.697149999999965,                       5190.697154999999839,
                       5190.697159999999712,                       5191.023654999999962,                       5191.350150000000212,                       5191.879649999999856,
                       5192.409150000000409,                       5193.455775000000358,                       5194.502400000000307,                       5194.907400000000052,
                       5195.312399999999798,                       5195.515400000000227,                       5195.718399999999747,                       5196.045399999999972,
                       5196.372400000000198,                       5196.372405000000072,                       5196.372409999999945,                       5196.777404999999817,
                       5197.182399999999689,                       5197.182404999999562,                       5197.182410000000345,                       5197.250904999999875,
                       5197.319400000000314,                       5197.319405000000188,                       5197.319410000000062,                       5197.869405000000370,
                       5198.419399999999769,                       5198.419404999999642,                       5198.419410000000425,                       5204.147904999999810,
                       5209.876400000000103,                       5210.373150000000351,                       5210.869899999999689,                       5211.409899999999652,
                       5211.949899999999616,                       5213.284749999999804,                       5214.619599999999991,                       5214.867599999999584,
                       5215.115600000000086,                       5215.254350000000159,                       5215.393100000000231,                       5216.280600000000049,
                       5217.168099999999868,                       5217.939349999999649,                       5218.710600000000341,                       5219.659349999999904,
                       5220.608100000000377,                       5220.954424999999901,                       5221.300750000000335,                       5221.550750000000335,
                       5221.800750000000335,                       5222.050750000000335,                       5222.300750000000335,                       5223.113999999999578,
                       5223.927249999999731,                       5231.944000000000415,                       5239.960750000000189,                       5241.006524999999783,
                       5242.052300000000287,                       5242.297300000000178,                       5242.542300000000068,                       5243.028800000000047,
                       5243.515300000000025,                       5243.781275000000278,                       5244.047249999999622,                       5244.197250000000167,
                       5244.347249999999804,                       5244.497250000000349,                       5244.647249999999985,                       5244.921774999999798,
                       5245.196299999999610,                       5245.556800000000294,                       5245.917300000000068,                       5246.266800000000330,
                       5246.616299999999683,                       5246.859800000000178,                       5247.103299999999763,                       5247.247049999999945,
                       5247.390800000000127,                       5248.278299999999945,                       5249.165799999999763,                       5250.117044999999962,
                       5251.068290999999590,                       5252.633294999999634,                       5254.198300000000017,                       5255.763305000000400,
                       5257.328309000000445,                       5257.528299999999945,                       5257.728291000000354,                       5259.293295000000398,
                       5260.858299999999872,                       5262.423305000000255,                       5263.988309000000299,                       5264.183299999999690,
                       5264.378290999999990,                       5265.943295000000035,                       5267.508300000000418,                       5269.073304999999891,
                       5270.638308999999936,                       5270.828300000000127,                       5271.018291000000318,                       5272.583295000000362,
                       5274.148299999999836,                       5275.713305000000219,                       5277.278309000000263,                       5278.194655000000239,
                       5279.110999999999876,                       5279.249749999999949,                       5279.388500000000022,                       5280.275999999999840,
                       5281.163499999999658,                       5282.109744999999748,                       5283.055991000000176,                       5284.620995000000221,
                       5286.185999999999694,                       5287.751005000000077,                       5289.316009000000122,                       5289.399005000000216,
                       5289.481999999999971,                       5289.588995000000068,                       5289.695990999999594,                       5291.260994999999639,
                       5292.826000000000022,                       5294.391005000000405,                       5295.956009000000449,                       5296.316504999999779,
                       5296.676999999999680,                       5300.171999999999571,                       5303.667000000000371,                       5305.705374999999549,
                       5307.743749999999636,                       5308.136249999999563,                       5308.528750000000400,                       5309.114392000000407,
                       5309.700034999999843};

    Double_t x_madx[]={-0.022728427212106,                       -0.026573484188917,                       -0.030418557822177,                       -0.039991467484264,
                       -0.049564453788793,                       -0.050069420361482,                       -0.050574296780747,                       -0.054349426422964,
                       -0.058124532750680,                       -0.060357783278068,                       -0.062591057745294,                       -0.076378761276200,
                       -0.090166527608473,                       -0.102764271110012,                       -0.115362047933835,                       -0.120438636731063,
                       -0.125515285284237,                       -0.129535997672883,                       -0.133556674551950,                       -0.133556753064254,
                       -0.133556831576522,                       -0.145047710575398,                       -0.156538628053597,                       -0.156538612815026,
                       -0.156538685378209,                       -0.157382226327567,                       -0.158225672591739,                       -0.158225744711667,
                       -0.158225816831591,                       -0.162246409450908,                       -0.166267077790823,                       -0.172787593045511,
                       -0.179308087231086,                       -0.192196765429471,                       -0.205085429283513,                       -0.210072792685941,
                       -0.215060113381225,                       -0.217560023973037,                       -0.220059851973005,                       -0.224086673702168,
                       -0.228113533225661,                       -0.228113586560711,                       -0.228113639895726,                       -0.233100965199632,
                       -0.238088236958173,                       -0.238088287544047,                       -0.238088338129887,                       -0.238931819924103,
                       -0.239775330993105,                       -0.239775381108142,                       -0.239775431223172,                       -0.246548340209630,
                       -0.253321296140047,                       -0.253321342492965,                       -0.253321388865060,                       -0.323864931914523,
                       -0.394408435590366,                       -0.400525652682990,                       -0.406642908640927,                       -0.413292727655655,
                       -0.419942524102291,                       -0.436380528503826,                       -0.452818525177476,                       -0.455872526995636,
                       -0.458926522833714,                       -0.460635190466858,                       -0.462343801129857,                       -0.473272919570095,
                       -0.484202044337046,                       -0.507041600357206,                       -0.529881185007865,                       -0.557977206406522,
                       -0.586073233503378,                       -0.596329204112667,                       -0.606585191992421,                       -0.613988632448542,
                       -0.621392052275018,                       -0.628795489121275,                       -0.636198927888863,                       -0.660282259467974,
                       -0.684365617991780,                       -0.921771458140324,                       -1.159177277810282,                       -1.190146556012369,
                       -1.221115863528267,                       -1.228371223750518,                       -1.235626581157799,                       -1.250033670888310,
                       -1.264440745061257,                       -1.272317250560714,                       -1.280193757073902,                       -1.284635813099939,
                       -1.289077866462157,                       -1.293519940551122,                       -1.297962002254059,                       -1.306091692630468,
                       -1.314221398081604,                       -1.324897142426019,                       -1.335572894142961,                       -1.345922895212984,
                       -1.356272897217983,                       -1.363483829034094,                       -1.370694769672466,                       -1.374951749572011,
                       -1.379208719129373,                       -1.405490896572119,                       -1.431773074377488,                       -1.411289447720157,
                       -1.390805798501584,                       -1.357118189180240,                       -1.323430560086443,                       -1.289719342221892,
                       -1.256008116992649,                       -1.251701619683206,                       -1.247395125374978,                       -1.213706221122331,
                       -1.180017316212012,                       -1.146307375690078,                       -1.112597432728055,                       -1.108398603235782,
                       -1.104199770925061,                       -1.070509590167247,                       -1.036819413363518,                       -1.003110747558533,
                       -0.969402081131856,                       -0.965310917229703,                       -0.961219752903968,                       -0.927528301269797,
                       -0.893836845087334,                       -0.860129458165354,                       -0.826422069877422,                       -0.806689946160769,
                       -0.786957822553024,                       -0.783970048746739,                       -0.780982279472034,                       -0.761871296018659,
                       -0.742760329537595,                       -0.747494694145131,                       -0.752229064801607,                       -0.760065987505823,
                       -0.767902909664833,                       -0.775726158233366,                       -0.783549405285889,                       -0.783964657937042,
                       -0.784379884424968,                       -0.784915226307985,                       -0.785450575152985,                       -0.793287788641925,
                       -0.801125010142071,                       -0.808947945186216,                       -0.816770877340590,                       -0.818574539457614,
                       -0.820378246534546,                       -0.837864870259442,                       -0.855351435901682,                       -0.865550129167210,
                       -0.875748716481268,                       -0.877712514880960,                       -0.879676371332692,                       -0.882606491971069,
                       -0.885536693347593};

    TGraph* gr_madx = new TGraph(161,s_madx,x_madx);
    gr_madx->SetName("gr_madx");
    gr_madx->SetLineColor(kBlack);
    gr_madx->SetLineWidth(3);
    //----------------------------------------------------------------//
    TFile* _file0 = TFile::Open("accelerator.root");

    const int nGraph = 100;

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

    TGraph* gr_lambda0_x[nGraph];
    TGraph* gr_lambda0_y[nGraph];
    TGraph* gr_pion_p_x[nGraph];
    TGraph* gr_pion_p_y[nGraph];
    TGraph* gr_proton_x[nGraph];
    TGraph* gr_proton_y[nGraph];
    TGraph* gr_pion_m_x[nGraph];
    TGraph* gr_pion_m_y[nGraph];

    TMultiGraph* mg_x = new TMultiGraph();
    TMultiGraph* mg_lambda0_x = new TMultiGraph();
    TMultiGraph* mg_pion_p_x = new TMultiGraph();
    TMultiGraph* mg_proton_x = new TMultiGraph();
    TMultiGraph* mg_pion_m_x = new TMultiGraph();
    mg_x->SetName("mg_x");
    mg_lambda0_x->SetName("mg_lambda0_x");
    mg_pion_p_x->SetName("mg_pion_p_x");
    mg_proton_x->SetName("mg_proton_x");
    mg_pion_m_x->SetName("mg_pion_m_x");

    TMultiGraph* mg_y = new TMultiGraph();
    TMultiGraph* mg_lambda0_y = new TMultiGraph();
    TMultiGraph* mg_pion_p_y = new TMultiGraph();
    TMultiGraph* mg_proton_y = new TMultiGraph();
    TMultiGraph* mg_pion_m_y = new TMultiGraph();
    mg_y->SetName("mg_y");
    mg_lambda0_y->SetName("mg_lambda0_y");
    mg_pion_p_y->SetName("mg_pion_p_y");
    mg_proton_y->SetName("mg_proton_y");
    mg_pion_m_y->SetName("mg_pion_m_y");

    for(Int_t i = 0; i < nGraph; i++)
    {
        TString gr_lambda0_x_name = "gr_lambda0_x_";
        TString gr_lambda0_y_name = "gr_lambda0_y_";
        gr_lambda0_x_name += i;
        gr_lambda0_y_name += i;
        gr_lambda0_x[i] = (TGraph*)_file0->Get(gr_lambda0_x_name.Data());
        gr_lambda0_y[i] = (TGraph*)_file0->Get(gr_lambda0_y_name.Data());
        gr_lambda0_x[i]->SetLineColor(kMagenta);
        gr_lambda0_y[i]->SetLineColor(kMagenta);
        gr_lambda0_x[i]->SetMarkerStyle(20);
        gr_lambda0_y[i]->SetMarkerStyle(20);
        gr_lambda0_x[i]->SetLineStyle(9);
        gr_lambda0_y[i]->SetLineStyle(9);

        TString gr_pion_p_x_name = "gr_pion_p_x_";
        TString gr_pion_p_y_name = "gr_pion_p_y_";
        gr_pion_p_x_name += i;
        gr_pion_p_y_name += i;
        gr_pion_p_x[i] = (TGraph*)_file0->Get(gr_pion_p_x_name.Data());
        gr_pion_p_y[i] = (TGraph*)_file0->Get(gr_pion_p_y_name.Data());
        gr_pion_p_x[i]->SetLineColor(kGreen);
        gr_pion_p_y[i]->SetLineColor(kGreen);

        TString gr_proton_x_name = "gr_proton_x_";
        TString gr_proton_y_name = "gr_proton_y_";
        gr_proton_x_name += i;
        gr_proton_y_name += i;
        gr_proton_x[i] = (TGraph*)_file0->Get(gr_proton_x_name.Data());
        gr_proton_y[i] = (TGraph*)_file0->Get(gr_proton_y_name.Data());
        gr_proton_x[i]->SetLineColor(kRed);
        gr_proton_y[i]->SetLineColor(kRed);

        TString gr_pion_m_x_name = "gr_pion_m_x_";
        TString gr_pion_m_y_name = "gr_pion_m_y_";
        gr_pion_m_x_name += i;
        gr_pion_m_y_name += i;
        gr_pion_m_x[i] = (TGraph*)_file0->Get(gr_pion_m_x_name.Data());
        gr_pion_m_y[i] = (TGraph*)_file0->Get(gr_pion_m_y_name.Data());
        gr_pion_m_x[i]->SetLineColor(kBlue);
        gr_pion_m_y[i]->SetLineColor(kBlue);

        mg_x->Add(gr_lambda0_x[i]);
        mg_lambda0_x->Add(gr_lambda0_x[i]);

        mg_x->Add(gr_pion_p_x[i]);
        mg_pion_p_x->Add(gr_pion_p_x[i]);

        mg_x->Add(gr_proton_x[i]);
        mg_proton_x->Add(gr_proton_x[i]);

        mg_x->Add(gr_pion_m_x[i]);
        mg_pion_m_x->Add(gr_pion_m_x[i]);

        mg_y->Add(gr_lambda0_y[i]);
        mg_lambda0_y->Add(gr_lambda0_y[i]);

        mg_y->Add(gr_pion_p_y[i]);
        mg_pion_p_y->Add(gr_pion_p_y[i]);

        mg_y->Add(gr_proton_y[i]);
        mg_proton_y->Add(gr_proton_y[i]);

        mg_y->Add(gr_pion_m_y[i]);
        mg_pion_m_y->Add(gr_pion_m_y[i]);
    }

    /*mg_x->Add(gr_madx);
    mg_x->Add(gr_accl);

    mg_lambda0_x->Add(gr_madx);
    mg_lambda0_x->Add(gr_accl);

    mg_pion_p_x->Add(gr_madx);
    mg_pion_p_x->Add(gr_accl);

    mg_proton_x->Add(gr_madx);
    mg_proton_x->Add(gr_accl);

    mg_pion_m_x->Add(gr_madx);
    mg_pion_m_x->Add(gr_accl);*/

    TLine* line_xrph_51937_ua9 = new TLine(_xrph_51937_ua9_pos,-1,_xrph_51937_ua9_pos,1);
    TLine* line_xrph_52202_ua9 = new TLine(_xrph_52202_ua9_pos,-1,_xrph_52202_ua9_pos,1);

    TCanvas* c_1_x = new TCanvas("c_1_x","c_1_x");
    c_1_x->cd();
    mg_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_lambda0_x = new TCanvas("c_lambda0_x","c_lambda0_x");
    c_lambda0_x->cd();
    mg_lambda0_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_pion_p_x = new TCanvas("c_pion_p_x","c_pion_p_x");
    c_pion_p_x->cd();
    mg_pion_p_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_proton_x = new TCanvas("c_proton_x","c_proton_x");
    c_proton_x->cd();
    mg_proton_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_pion_m_x = new TCanvas("c_pion_m_x","c_pion_m_x");
    c_pion_m_x->cd();
    mg_pion_m_x->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_1_y = new TCanvas("c_1_y","c_1_y");
    c_1_y->cd();
    mg_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_lambda0_y = new TCanvas("c_lambda0_y","c_lambda0_y");
    c_lambda0_y->cd();
    mg_lambda0_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_pion_p_y = new TCanvas("c_pion_p_y","c_pion_p_y");
    c_pion_p_y->cd();
    mg_pion_p_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_proton_y = new TCanvas("c_proton_y","c_proton_y");
    c_proton_y->cd();
    mg_proton_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");

    TCanvas* c_pion_m_y = new TCanvas("c_pion_m_y","c_pion_m_y");
    c_pion_m_y->cd();
    mg_pion_m_y->Draw("APL");
    line_xrph_51937_ua9->Draw("same");
    line_xrph_52202_ua9->Draw("same");
}
