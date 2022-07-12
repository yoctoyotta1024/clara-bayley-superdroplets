// Author: Clara Bayley
// File: superdrops.cpp
/* Functionality file for superdrops class */

#include "superdrop.hpp"

Superdrop::Superdrop(double aEps, double aR, double aM_sol, double aRho_l, double aRho_sol,
										 double aMr_sol, double aIonic) : Common2AllSuperdrops(aRho_l, aRho_sol,
																																					 aMr_sol, aIonic)
/* constructor function with parent and instance initialised */
{
	eps = aEps;
	r = aR;
	m_sol = aM_sol;
	setSuperdropInitials(aEps, aR, aM_sol);
}

void Superdrop::setSuperdropInitials(double aEps, double aR, double aM_sol)
/* Combined Setter functions to set private attributes */
{
	eps0 = aEps;
	r0 = aR;
	m_sol0 = aM_sol;
}

double Superdrop::dry_r()
/* calculate radius of dry droplet, ie.
 radius if it was entirely made of solute */
{
	return pow(3 * m_sol / (4 * M_PI * getRho_sol()), 1.0 / 3);
}

double Superdrop::mass()
/* calculate total mass of droplet
	mass = (water + dry areosol)  */
{
	double m;
	m = m_sol * (1.0 - getRho_l() / getRho_sol()); // mass contribution of solute
	m = 4.0 / 3 * M_PI * getRho_l() * pow(r, 3.0) + m;

	return m;
}

double Superdrop::vol()
/* volume of droplet */
{
	return 4.0 / 3 * M_PI * pow(r, 3.0);
}

// double Superdrop::rhoeff()
// 	/* calculates effective density of droplet
// 	so mass_droplet = m = 4/3*pi*r^3 * rhoeff */
// {
// 	double effsol; //effect of solute on density
// 	effsol = 1.0 - getRho_l()/getRho_sol();
// 	effsol = 3*m_sol/(4.0*M_PI*pow(r, 3.0))*effsol;

// 	return getRho_l() + effsol;
// }

// double Superdrop::m_w()
// 	/* mass of only water in droplet */
// {
// 	double v_w;
// 	v_w = 4.0/3 * M_PI *  (pow(r,3.0) - pow(dry_r, 3.0));

// 	return getRho_l() * v_w;
// }

double Superdrop::akohler_factor(const double temp)
/* calculate a in raoult factor (exp^(a/r)) to
account for effect of dissolved solute
on radial growth of droplet. Using eq.6.24
and eq.6.22 of lohmann, luond
and mahrt intro 2 clouds textbook */
{
	static const double akoh = 3.3e-7 / (dlc::TEMP0 * dlc::R0); // used in dimensionless kohler factor a calculation

	return akoh / temp; // 2*surface_tens (=0.0756N/m)/(rho_l*rgas_v*temp) [eq.6.24]
}

double Superdrop::bkohler_factor()
/* calculate b in kelvin factor (1-b/r^3)
to account for curvature on radial growth
of droplet. Using eq.6.24 and eq.6.22 of lohmann,
luond and mahrt intro 2 clouds textbook */
{
	static const double bkoh = 43e-6 * dlc::RHO0 / dlc::MR0; // used in dimensionless kohler factor b calculation

	return bkoh * m_sol * getIonic() / getMr_sol();
}

void Superdrop::diffusion_factors(double &fkl, double &fdl,
																	const double &p, const double &temp, const double &psat)
/* Calculate dimensionless Fkl and Fdl
* heat and vapour diffusion factors in
equation for radial growth of droplets
according to eq.7.27 lohmann, luond and
mahrt intro 2 clouds textbook
*/
{
	double Thermk, Diffuse_v;
	const double Temp = temp * dlc::TEMP0;
	const double P = p * dlc::P0;
	const double Psat = psat * dlc::P0;

	Thermk = 7.11756e-5 * pow(Temp, 2.0) + Temp * 4.38127686e-3; // [eq.7.24]
	Diffuse_v = (4.012182971e-5 / P * pow(Temp, 1.94)) / DC::RGAS_V;

	fkl = (DC::LATENT_RGAS_V / Temp - 1) * DC::LATENT_V / (Thermk * dlc::F0);
	fdl = Temp / (Diffuse_v * Psat) / dlc::F0;
}

double Superdrop::condensation_growth(const double &p, const double &temp,
																			const double &psat, const double &s_ratio, const double &delt)
/* radial growth/shrink of each superdroplet due to
condensation and diffusion of water vapour.
[eq.7.27 lohmann intro 2 clouds] */
{

	double fkl, fdl, delr;
	const double a = akohler_factor(temp);
	const double b = bkohler_factor();

	fkl = 0;
	fdl = 0;
	diffusion_factors(fkl, fdl, p, temp, psat);

	delr = (s_ratio - 1 - a / r + b / pow(r, 3.0)) / (dlc::Rho_l * (fkl + fdl) * r) * delt;
	// delr = (s_ratio-1) / (dlc::Rho_l * (fkl+fdl) * r) * delt;

	/*  if droplets are dry, do not shrink further */
	if (r <= dry_r() && delr <= 0.0)
	{
		delr = 0.0;
	}

	r += delr;

	return eps * delr;
}
