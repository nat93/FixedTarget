#define MagClass_cxx
#include "MagClass.h"

using namespace std;

MagClass::MagClass()
{}

MagClass::~MagClass()
{}

// Drift Matrix First Order
void MagClass::GetDriftMatrixR(Double_t L, Double_t **M)
{
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG DRIFT LENGTH!!! ##"<<endl;
        assert(0);
    }

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            if(i == j)
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;
            }
        }
    }

    M[0][1] = L;
    M[2][3] = L;
}

// Drift Matrix Second Order
void MagClass::GetDriftMatrixT(Double_t L, Double_t ***M)
{
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG DRIFT LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            for(Int_t k = 0; k < _mtrx_size; k++)
            {
                M[i][j][k] = 0.0;
            }
        }
    }

    M[4][1][1] = (-1)*L/2.0; M[0][5][1] = (-1)*L/2.0; M[4][3][3] = (-1)*L/2.0;
    M[2][5][3] = (-1)*L/2.0; M[0][1][5] = (-1)*L/2.0; M[2][3][5] = (-1)*L/2.0;
}

// Bending in Horizontal direction for Dipole Rectangular Magnet Matrix First Order
void MagClass::GetDipoleRectangularMatrixHorizontalR(Double_t Phi, Double_t Q, Double_t L, Double_t **M)
{
    // Check the angle
    if(Phi == 0)
    {
        cout<<"## ERROR: WRONG DIPOLE ANGLE!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    // Check the charge
    if(Q == 0)
    {
        cout<<"## ERROR: WRONG CHARGE!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            if(i == j)
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;
            }
        }
    }

    Double_t** MQ_X = new Double_t*[_mtrx_size/2];
    Double_t** MD_X = new Double_t*[_mtrx_size/2];
    Double_t** MR_X = new Double_t*[_mtrx_size/2];
    Double_t** M1_X = new Double_t*[_mtrx_size/2];
    Double_t** M2_X = new Double_t*[_mtrx_size/2];
    Double_t** MM_X = new Double_t*[_mtrx_size/2];

    Double_t** MQ_Y = new Double_t*[_mtrx_size/2];
    Double_t** MR_Y = new Double_t*[_mtrx_size/2];
    Double_t** M1_Y = new Double_t*[_mtrx_size/2];
    Double_t** MM_Y = new Double_t*[_mtrx_size/2];

    for(Int_t ii = 0; ii < _mtrx_size/2; ii++)
    {
        MQ_X[ii] = new Double_t[_mtrx_size/2];
        MD_X[ii] = new Double_t[_mtrx_size/2];
        MR_X[ii] = new Double_t[_mtrx_size/2];
        M1_X[ii] = new Double_t[_mtrx_size/2];
        M2_X[ii] = new Double_t[_mtrx_size/2];
        MM_X[ii] = new Double_t[_mtrx_size/2];

        MQ_Y[ii] = new Double_t[_mtrx_size/2];
        MR_Y[ii] = new Double_t[_mtrx_size/2];
        M1_Y[ii] = new Double_t[_mtrx_size/2];
        MM_Y[ii] = new Double_t[_mtrx_size/2];
    }

    Double_t Rho = L/Phi; // [m]
    Double_t K = Q/TMath::Power(Rho,2); // [m-2]]
    Double_t Delta = Phi/2.0; // [rad]

    if(K > 0)
    {
        // QUAD DEFOC X
        MQ_X[0][0] = 1.0;                     MQ_X[0][1] = 0.0; MQ_X[0][2] = 0.0;
        MQ_X[1][0] = TMath::Tan(Delta)/Rho;   MQ_X[1][1] = 1.0; MQ_X[1][2] = 0.0;
        MQ_X[2][0] = 0.0;                     MQ_X[2][1] = 0.0; MQ_X[2][2] = 1.0;

        // SECTOR DIPOLE X
        MD_X[0][0] = 1.0;                         MD_X[0][1] = 0.0; MD_X[0][2] = 0.0;
        MD_X[1][0] = (-1)*TMath::Sin(Phi)/Rho;    MD_X[1][1] = 1.0; MD_X[1][2] = TMath::Sin(Phi);
        MD_X[2][0] = 0.0;                         MD_X[2][1] = 0.0; MD_X[2][2] = 1.0;

        // DRIFT X
        MR_X[0][0] = 1.0;  MR_X[0][1] = L/2.0;  MR_X[0][2] = 0.0;
        MR_X[1][0] = 0.0;  MR_X[1][1] = 1.0;    MR_X[1][2] = 0.0;
        MR_X[2][0] = 0.0;  MR_X[2][1] = 0.0;    MR_X[2][2] = 1.0;

        MatrixMultiplication(MR_X,MQ_X,M1_X);
        MatrixMultiplication(MD_X,M1_X,M2_X);
        MatrixMultiplication(MR_X,M2_X,M1_X);
        MatrixMultiplication(MQ_X,M1_X,MM_X);

        // QUAD FOC Y
        MQ_Y[0][0] = 1.0;                      MQ_Y[0][1] = 0.0; MQ_Y[0][2] = 0.0;
        MQ_Y[1][0] = -TMath::Tan(Delta)/Rho;   MQ_Y[1][1] = 1.0; MQ_Y[1][2] = 0.0;
        MQ_Y[2][0] = 0.0;                      MQ_Y[2][1] = 0.0; MQ_Y[2][2] = 1.0;

        // DRIFT Y
        MR_Y[0][0] = 1.0;  MR_Y[0][1] = L;      MR_Y[0][2] = 0.0;
        MR_Y[1][0] = 0.0;  MR_Y[1][1] = 1.0;    MR_Y[1][2] = 0.0;
        MR_Y[2][0] = 0.0;  MR_Y[2][1] = 0.0;    MR_Y[2][2] = 1.0;

        MatrixMultiplication(MR_Y,MQ_Y,M1_Y);
        MatrixMultiplication(MQ_Y,M1_Y,MM_Y);
    }
    else
    {
        // QUAD FOC X
        MQ_X[0][0] = 1.0;                      MQ_X[0][1] = 0.0; MQ_X[0][2] = 0.0;
        MQ_X[1][0] = -TMath::Tan(Delta)/Rho;   MQ_X[1][1] = 1.0; MQ_X[1][2] = 0.0;
        MQ_X[2][0] = 0.0;                      MQ_X[2][1] = 0.0; MQ_X[2][2] = 1.0;

        // SECTOR DIPOLE X
        MD_X[0][0] = 1.0;                     MD_X[0][1] = 0.0; MD_X[0][2] = 0.0;
        MD_X[1][0] = TMath::SinH(Phi)/Rho;    MD_X[1][1] = 1.0; MD_X[1][2] = TMath::SinH(Phi);
        MD_X[2][0] = 0.0;                     MD_X[2][1] = 0.0; MD_X[2][2] = 1.0;

        // DRIFT X
        MR_X[0][0] = 1.0;  MR_X[0][1] = L/2.0;  MR_X[0][2] = 0.0;
        MR_X[1][0] = 0.0;  MR_X[1][1] = 1.0;    MR_X[1][2] = 0.0;
        MR_X[2][0] = 0.0;  MR_X[2][1] = 0.0;    MR_X[2][2] = 1.0;

        MatrixMultiplication(MR_X,MQ_X,M1_X);
        MatrixMultiplication(MD_X,M1_X,M2_X);
        MatrixMultiplication(MR_X,M2_X,M1_X);
        MatrixMultiplication(MQ_X,M1_X,MM_X);

        // QUAD DEFOC Y
        MQ_Y[0][0] = 1.0;                     MQ_Y[0][1] = 0.0; MQ_Y[0][2] = 0.0;
        MQ_Y[1][0] = TMath::Tan(Delta)/Rho;   MQ_Y[1][1] = 1.0; MQ_Y[1][2] = 0.0;
        MQ_Y[2][0] = 0.0;                     MQ_Y[2][1] = 0.0; MQ_Y[2][2] = 1.0;

        // DRIFT Y
        MR_Y[0][0] = 1.0;  MR_Y[0][1] = L;      MR_Y[0][2] = 0.0;
        MR_Y[1][0] = 0.0;  MR_Y[1][1] = 1.0;    MR_Y[1][2] = 0.0;
        MR_Y[2][0] = 0.0;  MR_Y[2][1] = 0.0;    MR_Y[2][2] = 1.0;

        MatrixMultiplication(MR_Y,MQ_Y,M1_Y);
        MatrixMultiplication(MQ_Y,M1_Y,MM_Y);
    }

    M[0][0] = MM_X[0][0]; M[0][1] = MM_X[0][1]; M[0][5] = MM_X[0][2];
    M[1][0] = MM_X[1][0]; M[1][1] = MM_X[1][1]; M[1][5] = MM_X[1][2];
    M[2][2] = MM_Y[0][0]; M[2][3] = MM_Y[0][1];
    M[3][2] = MM_Y[1][0]; M[3][3] = MM_Y[1][1];

    for(Int_t ii = 0; ii < _mtrx_size/2; ii++)
    {
        delete[] MQ_X[ii];
        delete[] MD_X[ii];
        delete[] MR_X[ii];
        delete[] M1_X[ii];
        delete[] M2_X[ii];
        delete[] MM_X[ii];

        delete[] MQ_Y[ii];
        delete[] MR_Y[ii];
        delete[] M1_Y[ii];
        delete[] MM_Y[ii];
    }
    delete[] MQ_X;
    delete[] MD_X;
    delete[] MR_X;
    delete[] M1_X;
    delete[] M2_X;
    delete[] MM_X;

    delete[] MQ_Y;
    delete[] MR_Y;
    delete[] M1_Y;
    delete[] MM_Y;
}

// Bending in Horizontal direction for Dipole Rectangular Magnet Matrix Second Order
void MagClass::GetDipoleRectangularMatrixHorizontalT(Double_t Phi, Double_t Q, Double_t L, Double_t ***M)
{
    // Check the angle
    if(Phi == 0)
    {
        cout<<"## ERROR: WRONG DIPOLE ANGLE!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    // Check the charge
    if(Q == 0)
    {
        cout<<"## ERROR: WRONG CHARGE!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            for(Int_t k = 0; k < _mtrx_size; k++)
            {
                M[i][j][k] = 0.0;
            }
        }
    }
}

// Quadrupole Matrix First Order
void MagClass::GetQuadrupoleMatrixR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M)
{
    // Check the gradient
    if(K0 == 0)
    {
        cout<<"## ERROR: WRONG QUADRUPOLE GRADIENT!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            if(i == j)
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;
            }
        }
    }

    Double_t dK = (-1)*K0*(P-P0)/P0; // [m-2]
    Double_t K = K0 + dK; // [m-2]
    Double_t Phi = L*TMath::Sqrt(TMath::Abs(K)); // [rad]

    if(K > 0)
    {
        M[0][0] = TMath::Cos(Phi);                          M[0][1] = TMath::Sin(Phi)/TMath::Sqrt(K);
        M[1][0] = (-1)*TMath::Sqrt(K)*TMath::Sin(Phi);      M[1][1] = TMath::Cos(Phi);
        M[2][2] = TMath::CosH(Phi);                         M[2][3] = TMath::SinH(Phi)/TMath::Sqrt(K);
        M[3][2] = TMath::Sqrt(K)*TMath::SinH(Phi);          M[3][3] = TMath::CosH(Phi);
    }
    else
    {
        K = TMath::Abs(K);

        M[0][0] = TMath::CosH(Phi);                     M[0][1] = TMath::SinH(Phi)/TMath::Sqrt(K);
        M[1][0] = TMath::Sqrt(K)*TMath::SinH(Phi);      M[1][1] = TMath::CosH(Phi);
        M[2][2] = TMath::Cos(Phi);                      M[2][3] = TMath::Sin(Phi)/TMath::Sqrt(K);
        M[3][2] = (-1)*TMath::Sqrt(K)*TMath::Sin(Phi);  M[3][3] = TMath::Cos(Phi);
    }
}

// Quadrupole Matrix in Thin-Lens Approximation First Order
void MagClass::GetQuadrupoleMatrixThinR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M)
{
    // Check the gradient
    if(K0 == 0)
    {
        cout<<"## ERROR: WRONG QUADRUPOLE GRADIENT!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            if(i == j)
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;
            }
        }
    }

    Double_t dK = (-1)*K0*(P-P0)/P0; // [m-2]
    Double_t K = K0 + dK; // [m-2]

    M[1][0] = (-1)*(K*L);
    M[3][2] = (K*L);
}

// Quadrupole Matrix in Thin-Lens Approximation Second Order
void MagClass::GetQuadrupoleMatrixThinT(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t ***M)
{
    // Check the gradient
    if(K0 == 0)
    {
        cout<<"## ERROR: WRONG QUADRUPOLE GRADIENT!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            for(Int_t k = 0; k < _mtrx_size; k++)
            {
                M[i][j][k] = 0.0;
            }
        }
    }
}

// Sextupole Matrix in Thin-Lens Approximation First Order
void MagClass::GetSextupoleMatrixThinR(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t **M)
{
    // Check the gradient
    if(K0 == 0)
    {
        cout<<"## ERROR: WRONG QUADRUPOLE GRADIENT!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            if(i == j)
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;
            }
        }
    }
}

// Sextupole Matrix in Thin-Lens Approximation Second Order
void MagClass::GetSextupoleMatrixThinT(Double_t K0, Double_t P, Double_t P0, Double_t L, Double_t ***M)
{
    // Check the gradient
    if(K0 == 0)
    {
        cout<<"## ERROR: WRONG QUADRUPOLE GRADIENT!!! ##"<<endl;
        assert(0);
    }
    // Check the length
    if(L < 0)
    {
        cout<<"## ERROR: WRONG MAGNET LENGTH!!! ##"<<endl;
        assert(0);
    }
    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            for(Int_t k = 0; k < _mtrx_size; k++)
            {
                M[i][j][k] = 0.0;
            }
        }
    }

    Double_t dK = (-1)*K0*(P-P0)/P0; // [m-2]
    Double_t K = K0 + dK; // [m-2]

    M[1][0][0] = (-1)*K*L/2.0;
    M[3][2][0] = K*L/2.0;
    M[3][0][2] = K*L/2.0;
    M[1][2][2] = K*L/2.0;
}

// Matrix multiplication function
void MagClass::MatrixMultiplication(Double_t **M_1, Double_t **M_2, Double_t **MULT)
{
    Int_t i, j, k;

    // Initializing elements of matrix MULT to 0.
    for(i = 0; i < _mtrx_size/2; i++)
    {
        for(j = 0; j < _mtrx_size/2; j++)
        {
            MULT[i][j] = 0;
        }
    }

    for(i = 0; i < _mtrx_size/2; i++)
    {
        for(j = 0; j < _mtrx_size/2; ++j)
        {
            for(k = 0; k < _mtrx_size/2; ++k)
            {
                MULT[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

// Matrix multiplication for the new coordinates
void MagClass::GetNewCoord(Double_t *X0, Double_t **R, Double_t ***T, Int_t order, Double_t *X)
{
    // Check the order
    if(order != 1 && order != 2)
    {
        cout<<"## ERROR: WRONG ORDER OF THE TRANSPORT MATRIX!!! ##"<<endl;
        assert(0);
    }

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        X[i] = 0.0;
    }

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; ++j)
        {
            X[i] += R[i][j]*X0[j];
            if(order > 1)
            {
                for(Int_t k = 0; k < _mtrx_size; k++)
                {
                    X[i] += T[i][j][k]*X0[j]*X0[k];
                }
            }
        }
    }
}

// Get new coordinates after passing the DRIFT
bool MagClass::GetNewCoordDrift(Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV)
{
    // Check the order
    if(order != 1 && order != 2)
    {
        cout<<"## ERROR: WRONG ORDER OF THE TRANSPORT MATRIX!!! ##"<<endl;
        assert(0);
    }

    Double_t** R = new Double_t*[_mtrx_size];
    Double_t*** T = new Double_t**[_mtrx_size];

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        R[i] = new Double_t[_mtrx_size];
        T[i] = new Double_t*[_mtrx_size];

        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            T[i][j] = new Double_t[_mtrx_size];
        }
    }

    GetDriftMatrixR(L,R);
    GetDriftMatrixT(L,T);
    GetNewCoord(X0,R,T,order,X);

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            delete[] T[i][j];
        }

        delete[] R[i];
        delete[] T[i];
    }
    delete[] R;
    delete[] T;

    return CheckApertureEllipse(X0,X,APH,APV);
}

// Get new coordinates after passing the DIPOLE
bool MagClass::GetNewCoordDipole(Double_t Phi, Double_t Q, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV)
{
    // Check the order
    if(order != 1 && order != 2)
    {
        cout<<"## ERROR: WRONG ORDER OF THE TRANSPORT MATRIX!!! ##"<<endl;
        assert(0);
    }

    Double_t** R = new Double_t*[_mtrx_size];
    Double_t*** T = new Double_t**[_mtrx_size];

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        R[i] = new Double_t[_mtrx_size];
        T[i] = new Double_t*[_mtrx_size];

        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            T[i][j] = new Double_t[_mtrx_size];
        }
    }

    GetDipoleRectangularMatrixHorizontalR(Phi,Q,L,R);
    GetDipoleRectangularMatrixHorizontalT(Phi,Q,L,T);
    GetNewCoord(X0,R,T,order,X);

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            delete[] T[i][j];
        }

        delete[] R[i];
        delete[] T[i];
    }
    delete[] R;
    delete[] T;

    return CheckApertureEllipse(X0,X,APH,APV);
}

// Get new coordinates after passing the QUADRUPOLE
bool MagClass::GetNewCoordQuadrupole(Double_t K0, Double_t P, Double_t P0, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV)
{
    // Check the order
    if(order != 1 && order != 2)
    {
        cout<<"## ERROR: WRONG ORDER OF THE TRANSPORT MATRIX!!! ##"<<endl;
        assert(0);
    }

    Double_t** R = new Double_t*[_mtrx_size];
    Double_t*** T = new Double_t**[_mtrx_size];

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        R[i] = new Double_t[_mtrx_size];
        T[i] = new Double_t*[_mtrx_size];

        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            T[i][j] = new Double_t[_mtrx_size];
        }
    }

    GetQuadrupoleMatrixThinR(K0,P,P0,L,R);
    GetQuadrupoleMatrixThinT(K0,P,P0,L,T);

    GetNewCoord(X0,R,T,order,X);

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            delete[] T[i][j];
        }

        delete[] R[i];
        delete[] T[i];
    }
    delete[] R;
    delete[] T;

    return CheckApertureEllipse(X0,X,APH,APV);
}

// Get new coordinates after passing the SEXTUPOLE
bool MagClass::GetNewCoordSextupole(Double_t K0, Double_t P, Double_t P0, Double_t L, Int_t order, Double_t *X0, Double_t *X, Double_t APH, Double_t APV)
{
    // Check the order
    if(order != 1 && order != 2)
    {
        cout<<"## ERROR: WRONG ORDER OF THE TRANSPORT MATRIX!!! ##"<<endl;
        assert(0);
    }

    Double_t** R = new Double_t*[_mtrx_size];
    Double_t*** T = new Double_t**[_mtrx_size];

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        R[i] = new Double_t[_mtrx_size];
        T[i] = new Double_t*[_mtrx_size];

        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            T[i][j] = new Double_t[_mtrx_size];
        }
    }

    GetSextupoleMatrixThinR(K0,P,P0,L,R);
    GetSextupoleMatrixThinT(K0,P,P0,L,T);
    GetNewCoord(X0,R,T,order,X);

    for(Int_t i = 0; i < _mtrx_size; i++)
    {
        for(Int_t j = 0; j < _mtrx_size; j++)
        {
            delete[] T[i][j];
        }

        delete[] R[i];
        delete[] T[i];
    }
    delete[] R;
    delete[] T;

    return CheckApertureEllipse(X0,X,APH,APV);
}

// Check aperture at the Entrance and Exit of the element
bool MagClass::CheckApertureEllipse(Double_t *X0, Double_t *X, Double_t APH, Double_t APV)
{
    if(TMath::Power(X0[0]/APH,2) + TMath::Power(X0[2]/APV,2) <= 1.0) // Entrance
    {
        if(TMath::Power(X[0]/APH,2) + TMath::Power(X[2]/APV,2) <= 1.0) // Exit
        {
            return true;
        }
    }
    //cout<<"## WARNING:: The particle is lost !!! ###"<<endl;
    return false;
}
