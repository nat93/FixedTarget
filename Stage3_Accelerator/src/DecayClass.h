#ifndef DecayClass_h
#define DecayClass_h

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
#include "TLorentzVector.h"

using namespace std;

class DecayClass
{
public:

   DecayClass();
   ~DecayClass();
   Bool_t threeBody(TLorentzVector *pProd,Double_t *mProd);
   Bool_t twoBody(TLorentzVector *pProd, Double_t *mProd);

private:
   void Boost(TLorentzVector& pIn, Double_t mIn, TLorentzVector& pOut);

   TRandom3* _rnd;
};

#endif

/*
 * http://www.helsinki.fi/~www_sefo/phenomenology/Schlippe_relativistic_kinematics.pdf
 * http://www.phys.ufl.edu/~avery/course/4390/f2015/lectures/relativistic_kinematics_2.pdf
 * http://www.phys.spbu.ru/content/File/Library/studentlectures/schlippe/pp05-04.pdf
 * https://github.com/mortenpi/pythia8/blob/master/src/ParticleDecays.cc
 */
