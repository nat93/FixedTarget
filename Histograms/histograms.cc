#include "includes.hh"
#include "constants.hh"

using namespace std;
using namespace constants;

TRandom3* _rnd;

void function_1();

void getThetaAndPhi(Float_t px, Float_t py, Float_t pz, Double_t &theta, Double_t &phi, Double_t &mag);
Double_t getParticleBetaGamma(Double_t pmag, Double_t mass);
Double_t getParticleFlightLength(Double_t ctau_LambdaC, Double_t pmag, Double_t mass);
Double_t getGamma(Double_t mass, Double_t pmag);
Double_t getBeta(Double_t gamma);
Double_t getParticleFlightTime(Double_t l, Double_t beta);
Double_t getTotEnergy( Double_t pmag, Double_t mass);
Double_t getAngleBetweenLatticeAndtrack(Double_t px, Double_t py,Double_t pz, Double_t nx, Double_t ny,Double_t nz);
Bool_t isNotDecayInTarget(Double_t l, Double_t  tau, Double_t beta, Double_t gamma, Double_t &lInTarget);
Double_t notDecayRate(Double_t l, Double_t N0, Double_t tau, Double_t beta, Double_t gamma);

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout<<endl;
        cout<<"--> 1 -- function_1() : histograms"<<endl;
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
      function_1();
      break;
    default:
      cout<<"--> Nothing to do =)"<<endl;
      break;
    }

    return 0;
}

void function_1()
{
    TString input_file_name     = "hardccbar.root";
    TString output_file_name    = "histograms.root";

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
    fChain1->Add(input_file_name.Data());

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

    cout<<"--> Input file: "<<input_file_name<<endl;
    Double_t nEntries = fChain1->GetEntries();
    cout<<"--> nEntries: "<<nEntries<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//
    TH1D* h_1 = new TH1D("h_1","L_{c} E [GeV/c^{2}]",1000,0.0,1000.0);
    TH1D* h_2 = new TH1D("h_2","L_{c} #theta [rad]",8000,-4.0,4.0);
    TH1D* h_3 = new TH1D("h_3","L_{c} #phi [rad]",8000,-4.0,4.0);
    TH1D* h_4 = new TH1D("h_4","L_{c} P [GeV/c]",1000,0.0,1000.0);
    TH1D* h_5 = new TH1D("h_5","L_{c} L [mm]",1000000,0.0,1000.0);
    TH1D* h_6 = new TH1D("h_6","L_{c} #gamma",1000,0.0,1000.0);
    TH1D* h_7 = new TH1D("h_7","L_{c} #beta",2000000,0.0,2.0);
    TH1D* h_8 = new TH1D("h_8","L_{c} Angle between Track and Plane",8000,-4.0,4.0);
    //--------------------------------------------------------------------------//

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
                Double_t theta, phi, pmag;
                getThetaAndPhi(_Px[i],_Py[i],_Pz[i],theta,phi,pmag);
                Double_t length = getParticleFlightLength(ctau_LambdaC,pmag,_mass[i]); // [mm]
                Double_t gamma = getGamma(_mass[i],pmag);
                Double_t beta = getBeta(gamma);
                Double_t angleBetweenLatticeAndtrack = getAngleBetweenLatticeAndtrack( _Px[i], _Py[i], _Pz[i], nx, ny, nz);

                h_1->Fill(_E[i]);
                h_2->Fill(theta);
                h_3->Fill(phi);
                h_4->Fill(pmag);
                h_5->Fill(length);
                h_6->Fill(gamma);
                h_7->Fill(beta);
                h_8->Fill(angleBetweenLatticeAndtrack);
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
    h_8->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

void getThetaAndPhi(Float_t px, Float_t py, Float_t pz, Double_t &theta, Double_t &phi, Double_t &mag)
{
    TVector3 *pv = new TVector3(px, py, pz);
    theta   = pv->Theta();
    phi     = pv->Phi();
    mag     = pv->Mag();
}

Double_t getParticleBetaGamma(Double_t pmag, Double_t mass)
{
    if(mass > 0.0)
    {
        return pmag/mass;
    }
    else
    {
        cout<<" ERROR ---> mass <= 0.0"<<endl
           <<" mass = "<<mass<<endl;
        assert(0);
    }
    return -999.0;
}

Double_t getParticleFlightLength(Double_t ctau_LambdaC, Double_t pmag,Double_t mass)
{
    if(mass > 0.0)
    {
        return ctau_LambdaC*pmag/mass;
    }
    else
    {
        cout<<" ERROR ---> mass <= 0.0"<<endl
           <<" mass = "<<mass<<endl;
        assert(0);
    }
    return -999.0;
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

Bool_t isNotDecayInTarget(Double_t l, Double_t  tau, Double_t beta, Double_t gamma, Double_t &lInTarget)
{
    Double_t N0 = 1.0;
    lInTarget = l - _rnd->Uniform(0.0,l);
    Double_t probNotDec = notDecayRate(lInTarget, N0, tau, beta, gamma);

    if(probNotDec > _rnd->Uniform(0.0,1.0))
    {
        return true;
    }
    return false;
}

Double_t notDecayRate(Double_t l, Double_t N0, Double_t tau, Double_t beta, Double_t gamma)
{
    Double_t t;
    if(l >= 0.0 && N0 >= 0.0 && tau >= 0.0 && beta >= 0.0 && gamma >= 1.0)
    {
        t = getParticleFlightTime(l, beta);
        if(t >= 0.0)
        {
            return N0*TMath::Exp(-t/(tau*gamma));
        }
        else
        {
            cout<<" ERROR ---> t < 0"<<endl
               <<" t = "<<t<<endl;
            assert(0);
        }
    }
    else
    {
        cout<<" ERROR: "<<endl<<" l = "
           <<l<<endl
          <<" beta = "<<beta<<endl
         <<" N0 = "<<N0<<endl
        <<" tau = "<<tau<<endl
        <<" gamma = "<<gamma<<endl;
        assert(0);
    }
    return -999.0;
}

Double_t getParticleFlightTime(Double_t l, Double_t beta)
{
    if(l >= 0.0 && beta > 0.0)
    {
        return l/(c_lightv*beta);
    }
    else
    {
        cout<<" ERROR ---> l or beta < 0"<<endl
           <<" l = "<<l<<endl
          <<" beta = "<<beta<<endl;
        assert(0);
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
