

#ifndef PHYSCONST_HPP
#define PHYSCONST_HPP

namespace physconst {

  const double m_to_Ang = 1e10;
  const double Ang_to_m = 1e-10;

  const double amu_to_kg = 1.660538782e-27; // kg

  const double eV_to_J = 1.60217653e-19; // J=Nm

  const double kB_J  = 1.412694e-23; // J/K
  const double kB_eV = 8.817343e-5; // eV/K

  const double GPa_to_eVA3 = 6.24150947961e-3; // GPa to eV/Ang^3
  const double eVA3_to_GPa = 160.217653; // eV/Ang^3 => GPa

  const double GPa_to_kbar = 10.0; // GPa to kbar
  const double kbar_to_GPa =  0.1; // 1 kbar = 10^8 Pa = 0.1 GPa

  const double Ekin_atomic_to_eV = 103.642685491; // amu * Ang^2/fs^2 => J => eV

}

#endif


