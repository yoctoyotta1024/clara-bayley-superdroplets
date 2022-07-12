// Author: Clara Bayley
// File: claras_SDconstants.hpp
/* File containing constants used by
SDM and/or CVODE ODE solver */

#ifndef CLARAS_SDCONSTANTS_HPP
#define CLARAS_SDCONSTANTS_HPP

/* Constants used in Superdrops Model NOTE: Any
variation in below quantities are ignored
(e.g. temperature depedences) */
namespace dimmed_constants
/* namespace containing values of constants with dimensions */
{
  const double G = 9.80665;                   // acceleration due to gravity [m/s^2]
  const double RGAS_UNIV = 8.314462618;       // universal molar gas constant [J/Kg/K]
  const double MR_WATER = 0.01801528;         // molecular mass of water [Kg/mol]
  const double MR_DRY = 0.0289647;            // molecular mass of dry air [Kg/mol]
  const double RGAS_DRY = RGAS_UNIV / MR_DRY; // specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
  const double RGAS_V = RGAS_UNIV / MR_WATER; // specific gas constant for water [J/Kg/K]

  const double LATENT_V = 2437300; // specific latent heat of vapourisation of water [J/Kg]  (at 300K)
  const double CP_DRY = 1004.9;    // specific heat capacity (dry) air at const. pressure [J/Kg/K] (taken at 300K)   = 1.400*cv_dry
  const double CP_V = 1864;        // specific heat capacity of water vapour [J/Kg/K] (at 300K)
  const double C_L = 4180;         // specific heat capacity of liquid water[J/Kg/K] (at 300K)

  const double RHO_DRY = 1.177; // density of dry air [Kg/m^3] (at 300K)
  const double RHO_L = 1000;    // density of liquid water condensing [kg/m^3] (at 300K)
  // const double RHO_L        = 996.57;            // density of liquid water condensing [kg/m^3] (at 300K)
  const double DYNVISC = 18.45 * 1e-6; // dynamic viscosity of air [Pa s] (at 300K)

  const double RHO_SOL = 2200;    // density of (dry) areosol [Kg/m^3]
  const double MR_SOL = 0.058443; //  molecular mass of areosol [Kg/mol]
  const int IONIC = 2;            //  degree ionic dissociation (van't Hoff factor) []

  const double LATENT_RGAS_V = LATENT_V / RGAS_V; // for fkl diffusin factors calculation
}

namespace dimless_constants
{
  /* characterstic time, velocity, temperate etc.
  scales used to make ODEs dimensionless constants
  with all CAPITALS constants have dimensions */
  namespace DC = dimmed_constants;

  const double W0 = 0.5;                           // characteristic velocity [m/s]
  const double TIME0 = 4000;                       // timescale [s]
  const double P0 = 100000;                        // pressure [Pa]
  const double TEMP0 = 273.15;                     // temperature [K]
  const double RHO0 = P0 / (DC::RGAS_DRY * TEMP0); // density [Kg/m^3]

  const double CP0 = DC::CP_DRY; // Heat capacity [J/Kg/K]
  const double MR0 = DC::MR_DRY; // molecular molar mass [Kg/mol]
  const double R0 = 1e-6;        // droplet radius lengthscale [m]
  // const double VOL0       = 1e6;                         // parcel volume scale [m^3]
  const double F0 = TIME0 / (RHO0 * R0 * R0); // droplet condensation-diffusion factors []

  /* de-dimensionalise constants for ODEs
  constants with only First letter capitalised
  are dimensionless */
  const double Mr_ratio = DC::MR_WATER / DC::MR_DRY;
  const double Cp_dry = DC::CP_DRY / CP0;
  const double Cp_v = DC::CP_V / CP0;
  const double C_l = DC::C_L / CP0;
  const double Latent_v = DC::LATENT_V / (TEMP0 * CP0);
  const double Rgas_dry = DC::RGAS_DRY / CP0;
  const double Rho_dry = DC::RHO_DRY / RHO0;
  const double Rho_l = DC::RHO_L / RHO0;
  const double Rho_sol = DC::RHO_SOL / RHO0;
  const double Mr_sol = DC::MR_SOL / MR0;
  const int IONIC = DC::IONIC;
}

#endif // CLARAS_SDCONSTANTS_HPP