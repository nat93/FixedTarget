#ifndef Constants_h
#define Constants_h

namespace Constants
{

// Position
const double _cry3_51799_ua9_pos = 5180.8295;   // [m]
const double _cry4_51799_ua9_pos = 5181.3245;   // [m]
const double _q1_51810_pos = 5185.1704;         // [m]
const double _q2_51910_pos = 5217.1681;         // [m]
const double _xrph_51937_ua9_pos = 5223.92725;  // [m]
const double _lsf_52005_pos = 5246.6163;        // [m]
const double _q3_52010_pos = 5249.1658;         // [m]
const double _mba_52030_pos = 5254.1983;        // [m]
const double _mba_52050_pos = 5260.8583;        // [m]
const double _mbb_52070_pos = 5267.5083;        // [m]
const double _mbb_52090_pos = 5274.1483;        // [m]
const double _q4_52110_pos = 5281.1635;         // [m]
const double _mbb_52130_pos = 5286.1860;        // [m]
const double _mbb_52150_pos = 5292.8260;        // [m]
const double _xrph_52202_ua9_pos = 5309.700035; // [m]

const double _beamAngleInitialAtCryPosition     = -423.573e-6; // [rad]
const double _beamPositionInitialAtCryPosition  = -0.0227284;  // [m]

const double _crystalAngle = -5000.0e-6;        // [rad]
const double _channeledAngleRange = 100.0e-6;   // [rad]

const double _nominal_momentum = 270.0;         // [GeV/c]

const double _dipole_length = 6.26;             // [m]
const double _quadrupole_length = 3.085;        // [m]
const double _sextupole_length = 0.423;         // [m]

//Q26
/*
const double _angle_dipole = 8.445141542e-3;    // [rad] kMBA, kMBB
const double _grad_quad_1 = 0.0144360370;       // [m-2] kQF1
const double _grad_quad_2 = -0.0144394607;      // [m-2] kQD
const double _grad_quad_3 = 0.0144360370;       // [m-2] kQF2
const double _grad_quad_4 = -0.0144394607;      // [m-2] kQD
const double _grad_sext_1 = 0.06325645859;      // [m-3] kLSFA
*/

//Q20

const double _angle_dipole = 8.445141542e-3;    // [rad] kMBA, kMBB
const double _grad_quad_1 = 0.01157957644;      // [m-2] kQF1
const double _grad_quad_2 = -0.01158097147;     // [m-2] kQD
const double _grad_quad_3 = 0.01157957644;      // [m-2] kQF2
const double _grad_quad_4 = -0.01158097147;     // [m-2] kQD
const double _grad_sext_1 = 0.04516855;         // [m-3] kLSFA

const double _crystal_offset_x = 0.0;       // [m]
const double _crystal_offset_y = 0.0;       // [m]
const int _order_transport_matrix = 1;

const double _aph = 0.075;                      // [m]
const double _apv = 0.075;                      // [m]

const bool _switch_magnets = true;
}

#endif
