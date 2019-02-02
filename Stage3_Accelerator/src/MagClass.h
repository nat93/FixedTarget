#ifndef MagClass_h
#define MagClass_h

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

class MagClass
{
public :

   MagClass();
   ~MagClass();

   void GetDriftMatrixR(Double_t L, Double_t **M);                                                                          // Drift Matrix First Order
   void GetDriftMatrixT(Double_t L, Double_t ***M);                                                                         // Drift Matrix Second Order
   void GetDipoleRectangularMatrixHorizontalR(Double_t Phi, Double_t Q, Double_t L, Double_t **M);                          // Bending in Horizontal direction for Dipole Rectangular Magnet Matrix First Order
   void GetDipoleRectangularMatrixHorizontalT(Double_t Phi, Double_t Q, Double_t L, Double_t ***M);                         // Bending in Horizontal direction for Dipole Rectangular Magnet Matrix Second Order
   void GetQuadrupoleMatrixR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M);                               // Quadrupole Matrix First Order
   void GetQuadrupoleMatrixThinR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M);                           // Quadrupole Matrix in Thin-Lens Approximation First Order
   void GetQuadrupoleMatrixThinT(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t ***M);                          // Quadrupole Matrix in Thin-Lens Approximation Second Order
   void GetSextupoleMatrixThinR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M);                            // Sextupole Matrix in Thin-Lens Approximation First Order
   void GetSextupoleMatrixThinT(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t ***M);                           // Sextupole Matrix in Thin-Lens Approximation Second Order
   void MatrixMultiplication(Double_t **M_1, Double_t **M_2, Double_t **MULT);                                              // Matrix multiplication function
   void GetNewCoord(Double_t *X0, Double_t **R, Double_t ***T, Int_t order, Double_t *X);                                   // Matrix multiplication for the new coordinates
   bool GetNewCoordDrift(Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV);                                               // Get new coordinates after passing the DRIFT
   bool GetNewCoordDipole(Double_t Phi, Double_t Q, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV);                    // Get new coordinates after passing the DIPOLE
   bool GetNewCoordQuadrupole(Double_t K0, Double_t P, Double_t P0, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV);    // Get new coordinates after passing the QUADRUPOLE
   bool GetNewCoordSextupole(Double_t K0, Double_t P, Double_t P0, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV);     // Get new coordinates after passing the SEXTUPOLE
   bool CheckApertureEllipse(Double_t *X0, Double_t *X, Double_t APH, Double_t APV);                                        // Check aperture at the Entrance and Exit of the element with ELLIPSOID shape

   static const Int_t _mtrx_size = 6;
};

#endif
