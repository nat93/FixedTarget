#ifndef constants_hh
#define constants_hh

namespace constants
{
  //--------------------------------------------------------//
  static const double ctau_LambdaC = 59.9/1000.0;           // [mm]
  static const double c_lightv     = 299792458000.0;        // [mm/s]
  static const double tau_LambdaC  = ctau_LambdaC/c_lightv; // [s]
  static const double targetLength = 3.0;                   // [mm]
  static const double targetWidth = 4.0;                    // [mm]
  static const double targetHeight = 50.0;                  // [mm]
  static const double targetCrystalGap = 0.5;               // [mm]
  static const double crystalLength = 20.0;                 // [mm]
  static const double crystalWidth = 4.0;                   // [mm]
  static const double crystalHeight = 50.0;                 // [mm]
  static const double crystalRadius = 4.0;                  // [m]
  //--------------------------------------------------------//
  // Normal vector to lattice                               //
  //--------------------------------------------------------//
  static const double nx = 1.0;                             //
  static const double ny = 0.0;                             //
  static const double nz = 0.0;                             //
  //--------------------------------------------------------//
}

#endif
