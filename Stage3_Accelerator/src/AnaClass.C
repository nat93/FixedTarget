#define AnaClass_cxx
#include "AnaClass.h"
#include "Constants.h"
#include "MagClass.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"

#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

void AnaClass::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L AnaClass.C
//      Root > AnaClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"--> Number of entries: "<<nentries<<endl;

   _h1 = new TH1D("h1","X position on the detector plane [All]",200000,-1000,1000);
   _h2 = new TH1D("h2","Y position on the detector plane [All]",200000,-1000,1000);
   _h3 = new TH2D("h3","YvsX position on the detector plane [All]",2000,-100,100,2000,-100,100);

   _h4 = new TH1D("h4","X position on the detector plane [Primaries]",2500,-1000,1000);
   _h5 = new TH1D("h5","Y position on the detector plane [Primaries]",2500,-1000,1000);
   _h6 = new TH2D("h6","YvsX position on the detector plane [Primaries]",1600,-40,40,1600,-40,40);

   _h7 = new TH1D("h7","X position on the detector plane [Secondaries]",2500,-1000,1000);
   _h8 = new TH1D("h8","Y position on the detector plane [Secondaries]",2500,-1000,1000);
   _h9 = new TH2D("h9","YvsX position on the detector plane [Secondaries]",1600,-40,40,1600,-40,40);

   _h10 = new TH1D("h10","Secondary Particle PDG encoding (not decay)",12000,-6000,6000);
   _h11 = new TH1D("h11","Secondary Particle PDG encoding (decay)",12000,-6000,6000);

   //_h12 = new TH2D("h12","PID vs ToF[ns] (All)",10000,0,500000,800,-4000,4000);
   //_h13 = new TH2D("h13","PID vs ToF[ns] (Primaries)",1000,0,1000,8000,-4000,4000);

   _h12 = new TH2D("h12","Charge*Mass[MeV] vs ToF[ns] (All)",1000,0,1000,4000,-2000,2000);
   _h13 = new TH2D("h13","Charge*Mass[MeV] vs ToF[ns] (All without decay)",1000,0,1000,4000,-2000,2000);

   _h14 = new TH2D("h14","PID vs ToF[ns] (Secondaries, not decay)",10000,0,500000,800,-4000,4000);
   _h15 = new TH2D("h15","PID vs ToF[ns] (Secondaries, decay)",10000,0,10000,800,-4000,4000);

   _h16 = new TH1D("h16","X position on the crystal plane [All]",20000,-1000,1000);
   _h17 = new TH1D("h17","Y position on the crystal plane [All]",20000,-1000,1000);
   _h18 = new TH1D("h18","Z position on the crystal plane [All]",100,0,10);

   _h19 = new TH1D("h19","X position on the crystal plane [Primaries]",20000,-1000,1000);
   _h20 = new TH1D("h20","Y position on the crystal plane [Primaries]",20000,-1000,1000);
   _h21 = new TH1D("h21","Z position on the crystal plane [Primaries]",100,0,10);

   _h22 = new TH1D("h22","X position on the crystal plane [Secondaries]",20000,-1000,1000);
   _h23 = new TH1D("h23","Y position on the crystal plane [Secondaries]",20000,-1000,1000);
   _h24 = new TH1D("h24","Z position on the crystal plane [Secondaries]",100,0,10);

   _h25 = new TH1D("h25","Momentum (GeV) [All]",3100,-10,300);
   _h26 = new TH1D("h26","Momentum (GeV) [Primaries]",3100,-10,300);
   _h27 = new TH1D("h27","Momentum (GeV) [Secondaries]",3100,-10,300);

   _h28 = new TH1D("h28","ToF K+/-[ns]",300000,0,300);
   _h29 = new TH1D("h29","ToF Pi+/-[ns]",300000,0,300);

   _h30 = new TH1D("h30","Px/Pz [rad] (All)",4000000,-2.0,2.0);
   _h31 = new TH1D("h31","Px/Pz [rad] (Primaries)",4000000,-2.0,2.0);
   _h32 = new TH1D("h32","Px/Pz [rad] (Secondaries)",4000000,-2.0,2.0);

   _h33 = new TH2D("h33","Momentum[GeV/c] vs ToF [ns] Pi+/-",3000,0,300,3000,-10,300);
   _h34 = new TH2D("h34","Momentum[GeV/c] vs ToF [ns] K+/-",3000,0,300,3000,-10,300);

   _h35 = new TH1D("h35","X position on the detector plane [Charge > 0]",200000,-1000,1000);
   _h36 = new TH1D("h36","X position on the detector plane [Charge < 0]",200000,-1000,1000);

   _h37 = new TH2D("h37","X vs L",1500,0,150,2000,-1000,1000);
   _h38 = new TH2D("h38","Y vs L",1500,0,150,2000,-1000,1000);

   Long64_t nbytes = 0, nb = 0;
   Long64_t nAll = 0, nPrimaries = 0, nSecondaries = 0;
   Long64_t decay = 0;
   Int_t nParticles = 0;

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

   Double_t DRIFTL0 = Constants::_xrph_52202_ua9_pos - Constants::_cry4_51799_ua9_pos;   // [m]
   Double_t DRIFTL1 = Constants::_q1_51810_pos-Constants::_cry4_51799_ua9_pos;           // [m]
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

   Double_t* coord_x0          = new Double_t[mtrx_size];
   Double_t* coord_x           = new Double_t[mtrx_size];
   Double_t* coord_temp_1_x    = new Double_t[mtrx_size];
   Double_t* coord_temp_2_x    = new Double_t[mtrx_size];

   nentries /= 10;
   Bool_t single = true;
   if(single){nentries = 1;}

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      nParticles = nParticles_C;
      if(single){nParticles = 1;}

      for(Int_t i = 0; i < nParticles; i++)
      {
          //if(Mass_C[i] < 1.322 && Mass_C[i] > 1.321) cout<<"mass: "<<Mass_C[i]<<" pid: "<<Pid_C[i]<<" ctau: "<<c*Tau_C[i]<<" ch: "<<Charge_C[i]<<endl;

          // Real cuts
          //if(Pz_C[i] <= 0) continue;

          // Not Real cuts
          //if(ParentID_C[i] != 0) continue;
          //if(Pz_C[i] <= 0 || Charge_C[i] == 0 || Z_C[i] < 8.49) continue;
          //if(Pz_C[i] <= 0 || Charge_C[i] != 0) continue;

          Double_t distance = (DRIFTL0)*1000.0; // [mm]
          Double_t p = TMath::Sqrt(Px_C[i]*Px_C[i] + Py_C[i]*Py_C[i] + Pz_C[i]*Pz_C[i]);
          Double_t tof = calculateFlightTime(distance, Mass_C[i], p);

          if(single)
          {
              Charge_C[i]   = 1;
              p             = 20.0;                 // [GeV/c]
              /*x*/     coord_x0[0] = 0.0;        // [m]
              /*x'*/    coord_x0[1] = -0.0012;     // [rad]
              /*y*/     coord_x0[2] = 0.0;        // [m]
              /*y'*/    coord_x0[3] = 0.0;        // [rad]
              /*l*/     coord_x0[4] = 0.0;
              /*dp*/    coord_x0[5] = 0.0;
          }
          else
          {
              /*x*/     coord_x0[0] = Constants::_crystal_offset_x + X_C[i]/1000.0; // [m]
              /*x'*/    coord_x0[1] = TMath::Tan(TMath::ATan(Px_C[i]/Pz_C[i]) - Constants::_crystalAngle); // [rad]
              /*y*/     coord_x0[2] = Constants::_crystal_offset_y + Y_C[i]/1000.0; // [m]
              /*y'*/    coord_x0[3] = Py_C[i]/Pz_C[i]; // [rad]
              /*l*/     coord_x0[4] = 0.0;
              /*dp*/    coord_x0[5] = (p - Constants::_nominal_momentum)/Constants::_nominal_momentum;
          }

          //-----------------------------------------//
          // Magnets section
          //-----------------------------------------//

          if(Constants::_switch_magnets && Charge_C[i] != 0)
          {
              // INITIAL
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_x0[0],coord_x0[1],0.0);

              // DRIFT1
              if(!magnet->GetNewCoordDrift(DRIFTL1,order,coord_x0,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // QUAD1: FOC in X, DEFOC in Y
              if(!magnet->GetNewCoordQuadrupole(Charge_C[i]*KQ1,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<Constants::_q1_51810_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q1_51810_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT2
              if(!magnet->GetNewCoordDrift(DRIFTL2,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // QUAD2: FOC in Y, DEFOC in X
              if(!magnet->GetNewCoordQuadrupole(Charge_C[i]*KQ2,p,P0,QL,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<Constants::_q2_51910_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_2_x[0],coord_temp_2_x[1],Constants::_q2_51910_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT3
              if(!magnet->GetNewCoordDrift(DRIFTL3,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_xrph_51937_ua9_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT4
              if(!magnet->GetNewCoordDrift(DRIFTL4,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // SEXT1
              if(!magnet->GetNewCoordSextupole(Charge_C[i]*KS1,p,P0,SL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_lsf_52005_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_lsf_52005_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT5
              if(!magnet->GetNewCoordDrift(DRIFTL5,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // QUAD3: FOC in X, DEFOC in Y
              if(!magnet->GetNewCoordQuadrupole(Charge_C[i]*KQ3,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_q3_52010_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q3_52010_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT6
              if(!magnet->GetNewCoordDrift(DRIFTL6,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE1: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mba_52030_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52030_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT7
              if(!magnet->GetNewCoordDrift(DRIFTL7,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE2: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mba_52050_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mba_52050_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT8
              if(!magnet->GetNewCoordDrift(DRIFTL8,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE3: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mbb_52070_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52070_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT9
              if(!magnet->GetNewCoordDrift(DRIFTL9,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE4: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mbb_52090_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52090_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT10
              if(!magnet->GetNewCoordDrift(DRIFTL10,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // QUAD4: FOC in Y, DEFOC in X
              if(!magnet->GetNewCoordQuadrupole(Charge_C[i]*KQ4,p,P0,QL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_q4_52110_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_q4_52110_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT11
              if(!magnet->GetNewCoordDrift(DRIFTL11,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE5: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mbb_52130_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52130_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT12
              if(!magnet->GetNewCoordDrift(DRIFTL12,order,coord_temp_1_x,coord_temp_2_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}

              // DIPOLE6: BEND & DEFOC in X, FOC in Y
              if(!magnet->GetNewCoordDipole(ADL,Charge_C[i],DL,order,coord_temp_2_x,coord_temp_1_x,APH,APV))
              {cout<<"## At : "<<Constants::_mbb_52150_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_temp_1_x[0],coord_temp_1_x[1],Constants::_mbb_52150_pos-Constants::_cry4_51799_ua9_pos);

              // DRIFT13
              if(!magnet->GetNewCoordDrift(DRIFTL13,order,coord_temp_1_x,coord_x,APH,APV))
              {cout<<"## At : "<<Constants::_xrph_52202_ua9_pos-Constants::_cry4_51799_ua9_pos<<" [m] ##"<<endl; continue;}
              if(single) printf("%10.15f , %10.15f , %10.15f\n",coord_x[0],coord_x[1],Constants::_xrph_52202_ua9_pos-Constants::_cry4_51799_ua9_pos);
          }
          else
          {
              // DRIFT0
              if(!magnet->GetNewCoordDrift(DRIFTL0,order,coord_x0,coord_x,APH,APV))
              {cout<<"## At : "<<"DRIFT"<<" [m] ##"<<endl; continue;}
          }
          //-----------------------------------------//

          // To investigate specified region
          //Double_t x_min = -40.0;
          //Double_t x_max = -20.0;
          //if(coord_x[0]*1000.0 < x_min || coord_x[0]*1000.0 > x_max) continue;

          _h1->Fill(coord_x[0]*1000.0);
          _h2->Fill(coord_x[2]*1000.0);
          _h3->Fill(coord_x[0]*1000.0,coord_x[2]*1000.0);
          _h12->Fill(tof,Charge_C[i]*Mass_C[i]*1000.0);
          _h16->Fill(X_C[i]);
          _h17->Fill(Y_C[i]);
          _h18->Fill(Z_C[i]);
          _h25->Fill(p);

          if(Mass_C[i] < 0.500 && Mass_C[i] > 0.480)
          {
              _h28->Fill(tof);
              _h34->Fill(tof,p);
          }
          if(Mass_C[i] < 0.150 && Mass_C[i] > 0.130)
          {
              _h29->Fill(tof);
              _h33->Fill(tof,p);
          }
          _h30->Fill(Px_C[i]/Pz_C[i] - Constants::_crystalAngle);
          nAll++;

          if(ParentID_C[i] == 0)
          {
              _h4->Fill(coord_x[0]*1000.0);
              _h5->Fill(coord_x[2]*1000.0);
              _h6->Fill(coord_x[0]*1000.0,coord_x[2]*1000.0);
              _h13->Fill(tof,Charge_C[i]*Mass_C[i]*1000.0);
              _h19->Fill(X_C[i]);
              _h20->Fill(Y_C[i]);
              _h21->Fill(Z_C[i]);
              _h26->Fill(p);
              _h31->Fill(Px_C[i]/Pz_C[i] - Constants::_crystalAngle);
              nPrimaries++;
          }
          else
          {
              if(ifNotDecay(distance, Tau_C[i], Mass_C[i], p))
              {
                  _h7->Fill(coord_x[0]*1000.0);
                  _h8->Fill(coord_x[2]*1000.0);
                  _h9->Fill(coord_x[0]*1000.0,coord_x[2]*1000.0);
                  _h10->Fill(Pid_C[i]);
                  _h13->Fill(tof,Charge_C[i]*Mass_C[i]*1000.0);
                  _h14->Fill(tof,Pid_C[i]);
                  _h22->Fill(X_C[i]);
                  _h23->Fill(Y_C[i]);
                  _h24->Fill(Z_C[i]);
                  _h27->Fill(p);
                  _h32->Fill(Px_C[i]/Pz_C[i] - Constants::_crystalAngle);
                  nSecondaries++;
              }
              else
              {
                  _h11->Fill(Pid_C[i]);
                  _h15->Fill(tof,Pid_C[i]);
                  decay++;
              }
          }
          if(Charge_C[i] > 0)
              _h35->Fill(coord_x[0]*1000.0);
          if(Charge_C[i] < 0)
              _h36->Fill(coord_x[0]*1000.0);
      }

      Long64_t div = nentries*0.01;
      if(div < 1) div = 1;
      if(jentry%div == 0)
      {
          printf("\r--> Progress: %3.f %%",100.0*jentry/nentries);
          fflush(stdout);
      }
   }
   cout<<endl;
   cout<<"--> Number of bytes: "<<nbytes<<" [bytes]"<<endl;
   cout<<"--> Secondaries (not decay): "<<100.0*nSecondaries/nAll
      <<"% | Primaries: "<<100.0*nPrimaries/nAll
     <<"% | Secondaries (decay): "<<100.0*decay/nAll
    <<"%"<<endl;

   //Free the array of pointers
   delete[] coord_x0;
   delete[] coord_x;
   delete[] coord_temp_1_x;
   delete[] coord_temp_2_x;
}

void AnaClass::WriteHistos(TString outputRootFile)
{
   TFile* outfile = new TFile(outputRootFile.Data(),"RECREATE");
   _h1->Write();
   _h2->Write();
   _h3->Write();
   _h4->Write();
   _h5->Write();
   _h6->Write();
   _h7->Write();
   _h8->Write();
   _h9->Write();
   _h10->Write();
   _h11->Write();
   _h12->Write();
   _h13->Write();
   _h14->Write();
   _h15->Write();
   _h16->Write();
   _h17->Write();
   _h18->Write();
   _h19->Write();
   _h20->Write();
   _h21->Write();
   _h22->Write();
   _h23->Write();
   _h24->Write();
   _h25->Write();
   _h26->Write();
   _h27->Write();
   _h28->Write();
   _h29->Write();
   _h30->Write();
   _h31->Write();
   _h32->Write();
   _h33->Write();
   _h34->Write();
   _h35->Write();
   _h36->Write();
   _h37->Write();
   _h38->Write();

   outfile->Close();
   cout<<"--> Histograms were written to: "<<outputRootFile<<endl;
}

Bool_t AnaClass::ifNotDecay(Double_t L, Double_t tau, Double_t m, Double_t p)
{
    Double_t N0 = 1.0;
    Double_t probNotDec = notDecayRate(L, N0, tau, m, p);
    if(probNotDec>_rnd->Uniform(0.0,1.0))
        return true;
    if(tau <= 0)
        return true;
    return false;
}

Double_t AnaClass::notDecayRate(Double_t l, Double_t N0, Double_t tau, Double_t m, Double_t p)
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

Double_t AnaClass::calculateGamma(Double_t m, Double_t p){
  if(p>=0.0){
    if(m>0.0){
      return TMath::Sqrt(1.0 + p*p/m/m);
    }
  }
  return -999.0;
}

Double_t AnaClass::calculateFlightTime(Double_t l, Double_t m, Double_t p){
  Double_t betta;
  if(l>0.0){
    if(m>=0.0){
      if(p>0.0){
    betta = calculateBetta(m, p);
    if(betta>0.0){
      return l/c/betta;
    }
      }
    }
  }
  return -999.0;
}

Double_t AnaClass::calculateBetta(Double_t m, Double_t p){
  if(p>=0.0){
    if(m>0.0){
      return 1.0/(TMath::Sqrt(1.0 + m*m/p/p));
    }
    else if( m == 0.0){
      return 1.0;
    }
    else{
      return -999.0;
    }
  }
  return -999.0;
}
