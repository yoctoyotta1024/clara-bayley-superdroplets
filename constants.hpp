
#ifndef CONSTANTS
#define CONSTANTS


// some constants used in Superdrops Model
/* NOTE: Any variation in below quantities are ignored (e.g. temperature depedences) */



namespace DropletConstants
{
const double RGAS_UNIV    = 8.314462618;       // universal molar gas constant
const double MR_WATER     = 0.01801528;        // molecular mass of water
const double MR_DRY       = 0.0289647;         // molecular mass (per 1 mole) of dry air (M_r)                 
const double MR_RATIO     = MR_WATER/MR_DRY;
const double RGAS_DRY     = RGAS_UNIV/MR_DRY;          // specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
const double RGAS_V       = RGAS_UNIV/MR_WATER;        // specific gas constant for water
 
const double CP_DRY       = 1004.9;            // specific heat capacity (dry) air at const. pressure [J/Kg/K] (taken at 300K)   = 1.400*cv_dry 
const double CP_V         = 1864;              // specific heat capacity of water vapour [J/Kg/K] (at 300K)
const double C_L          = 4180;              // specific heat capacity of liquid water[J/Kg/K] (at 300K)
const double LATENT_V     = 2437300;           // specific latent heat of vapourisation of water [J/Kg]  (at 300K)

const double RHO_DRY      = 1.177;             // density of dry air [Kg/m^3] (at 300K)
const double RHO_L        = 996.57;            // density of liquid water condensing [kg/m^3] (at 300K)
const double DYNVISC      = 18.45*1e-6;        // dynamic viscosity of air [Pa s] (at 300K)
}





#endif //CONSTANTS