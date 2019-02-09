#define DecayClass_cxx
#include "DecayClass.h"

using namespace std;

DecayClass::DecayClass()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int seed = tp.tv_sec*1000 + tp.tv_usec/1000;
    _rnd = new TRandom3(seed);
}

DecayClass::~DecayClass()
{}

// Do a two-body decay.
Bool_t DecayClass::twoBody(TLorentzVector *pProd, Double_t *mProd)
{
    // Masses.
    Double_t m0   = mProd[0];
    Double_t m1   = mProd[1];
    Double_t m2   = mProd[2];

    const Double_t mSafety = 0.0005;
    /*
     * Minimum mass difference required between the decaying mother mass and the sum of the daughter masses,
     * kept as a safety margin to avoid numerical problems in the decay generation.
     * (default = 0.0005; minimum = 0.; maximum = 0.01)
     */

    // Energies and absolute momentum in the rest frame.
    if (m1 + m2 + mSafety > m0) return kFALSE;

    Double_t e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
    Double_t e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
    Double_t pAbs = 0.5 * TMath::Sqrt( (m0 - m1 - m2) * (m0 + m1 + m2) * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;

    // Isotropic angles give three-momentum.
    Double_t cosTheta = 2. * _rnd->Uniform(0.0,1.0) - 1.;
    Double_t sinTheta = TMath::Sqrt(1. - cosTheta*cosTheta);
    Double_t phi      = 2. * TMath::Pi() * _rnd->Uniform(0.0,1.0);
    Double_t pX       = pAbs * sinTheta * TMath::Cos(phi);
    Double_t pY       = pAbs * sinTheta * TMath::Sin(phi);
    Double_t pZ       = pAbs * cosTheta;

    // Fill four-momenta and boost them away from mother rest frame.
    pProd[1].SetPxPyPzE(  pX,  pY,  pZ, e1);
    pProd[2].SetPxPyPzE( -pX, -pY, -pZ, e2);

    Boost(pProd[0],m0,pProd[1]);
    Boost(pProd[0],m0,pProd[2]);

    // Done.
    return kTRUE;
}

// Do a three-body decay (except Dalitz decays).
Bool_t DecayClass::threeBody(TLorentzVector *pProd, Double_t *mProd)
{
    // Mother and sum daughter masses. Fail if too close.
    Double_t m0      = mProd[0];
    Double_t m1      = mProd[1];
    Double_t m2      = mProd[2];
    Double_t m3      = mProd[3];
    Double_t mSum    = m1 + m2 + m3;
    Double_t mDiff   = m0 - mSum;

    const Double_t mSafety = 0.0005;
    /*
     * Minimum mass difference required between the decaying mother mass and the sum of the daughter masses,
     * kept as a safety margin to avoid numerical problems in the decay generation.
     * (default = 0.0005; minimum = 0.; maximum = 0.01)
     */

    if (mDiff < mSafety) return kFALSE;

    // Kinematical limits for 2+3 mass. Maximum phase-space weight.
    Double_t m23Min  = m2 + m3;
    Double_t m23Max  = m0 - m1;
    Double_t p1Max   = 0.5 * TMath::Sqrt( (m0 - m1 - m23Min) * (m0 + m1 + m23Min) * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0;
    Double_t p23Max  = 0.5 * TMath::Sqrt( (m23Max - m2 - m3) * (m23Max + m2 + m3) * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
    Double_t wtPSmax = 0.5 * p1Max * p23Max;

    // Begin loop over matrix-element corrections.
    Double_t wtPS, m23, p1Abs, p23Abs;

    // Pick an intermediate mass m23 flat in the allowed range.
    while( 1 )
    {
        m23 = m23Min + _rnd->Uniform(0.0,1.0) * mDiff;

        // Translate into relative momenta and find phase-space weight.
        p1Abs  = 0.5 * TMath::Sqrt( (m0 - m1 - m23) * (m0 + m1 + m23) * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0;
        p23Abs = 0.5 * TMath::Sqrt( (m23 - m2 - m3) * (m23 + m2 + m3) * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
        wtPS   = p1Abs * p23Abs;

        // If rejected, try again with new invariant masses.
        if( wtPS >= _rnd->Uniform(0.0,1.0) * wtPSmax ) break;
    }

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    Double_t cosTheta = 2. * _rnd->Uniform(0.0,1.0) - 1.;
    Double_t sinTheta = sqrt(1. - cosTheta*cosTheta);
    Double_t phi      = 2. * TMath::Pi() * _rnd->Uniform(0.0,1.0);
    Double_t pX       = p23Abs * sinTheta * TMath::Cos(phi);
    Double_t pY       = p23Abs * sinTheta * TMath::Sin(phi);
    Double_t pZ       = p23Abs * cosTheta;
    Double_t e2       = sqrt( m2*m2 + p23Abs*p23Abs);
    Double_t e3       = sqrt( m3*m3 + p23Abs*p23Abs);

    pProd[2].SetPxPyPzE( pX,  pY,  pZ, e2 );
    pProd[3].SetPxPyPzE( -pX,  -pY,  -pZ, e3 );

    // Set up m0 -> m1 + m23 isotropic in its rest frame.
    cosTheta        = 2. * _rnd->Uniform(0.0,1.0) - 1.;
    sinTheta        = TMath::Sqrt(1. - cosTheta*cosTheta);
    phi             = 2. * TMath::Pi() * _rnd->Uniform(0.0,1.0);
    pX              = p1Abs * sinTheta * cos(phi);
    pY              = p1Abs * sinTheta * sin(phi);
    pZ              = p1Abs * cosTheta;
    Double_t e1       = TMath::Sqrt( m1*m1 + p1Abs*p1Abs);
    Double_t e23      = TMath::Sqrt( m23*m23 + p1Abs*p1Abs);

    pProd[1].SetPxPyPzE( pX, pY, pZ, e1);

    // Boost 2 + 3 to the 0 rest frame.
    TLorentzVector prod23;
    prod23.SetPxPyPzE( -pX, -pY, -pZ, e23);

    Boost(prod23,m23,pProd[2]);
    Boost(prod23,m23,pProd[3]);

    // Boost 1 + 2 + 3 to the current frame.
    Boost(pProd[0],m0,pProd[1]);
    Boost(pProd[0],m0,pProd[2]);
    Boost(pProd[0],m0,pProd[3]);

    // Done.
    return kTRUE;
}

void DecayClass::Boost(TLorentzVector& pIn, Double_t mIn, TLorentzVector& pOut)
{
    Double_t xx = pOut.Px();
    Double_t yy = pOut.Py();
    Double_t zz = pOut.Pz();
    Double_t tt = pOut.E();

    Double_t betaX = pIn.Px() / pIn.E();
    Double_t betaY = pIn.Py() / pIn.E();
    Double_t betaZ = pIn.Pz() / pIn.E();
    Double_t gamma = pIn.E() / mIn;
    Double_t prod1 = betaX * xx + betaY * yy + betaZ * zz;
    Double_t prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
    xx += prod2 * betaX;
    yy += prod2 * betaY;
    zz += prod2 * betaZ;
    tt  = gamma * (tt + prod1);

    pOut.SetPxPyPzE(xx,yy,zz,tt);
}
