// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "CaSeMaConfig.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <algorithm>
#include <fstream>
#include <tclap/CmdLine.h>

#include "MPReal.hpp"

#include "ModelData.hpp"
#include "AnalyticalMoments.hpp"


casema::ModelData<mpfr::mpreal> qamarModel(const mpfr::mpreal& velocity)
{
	casema::ModelData<mpfr::mpreal> data;

	data.nComponents = 1;
	
	data.initialLiquidConcentration.push_back(mpfr::mpreal(0));
	data.initialSolidConcentration.push_back(mpfr::mpreal(0));

	data.colDispersion = mpfr::mpreal("0.002") / mpfr::mpreal(600000);
	data.colLength = mpfr::mpreal("0.017");
	data.colPorosity = mpfr::mpreal("0.4");

	data.parRadius = mpfr::mpreal("0.00004");
	data.parPorosity = mpfr::mpreal("0.333");

	data.filmDiffusion.push_back(mpfr::mpreal("0.01") / mpfr::mpreal(6000));
	data.particleDiffusion.push_back(mpfr::mpreal("1e-6") / (mpfr::mpreal(600000)*data.parPorosity));
	data.surfaceDiffusion.push_back(mpfr::mpreal(0));

	data.velocity = velocity;

	data.kineticBinding = false;

	data.bindingModel = casema::LINEAR;
	data.linearKA.push_back(mpfr::mpreal(25));
	data.linearKD.push_back(mpfr::mpreal(10));

	data.nInletSections = 2;
	data.sectionTimes.push_back(mpfr::mpreal(0));
	data.sectionTimes.push_back(mpfr::mpreal(1200));
	data.sectionTimes.push_back(mpfr::mpreal(10000));
	data.constCoeff.push_back(mpfr::mpreal(1000)); // 1 g / l in mol / m^3
	data.constCoeff.push_back(mpfr::mpreal(0));

	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));
	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));

	return data;
}

casema::ModelData<mpfr::mpreal> qamarModelCmMin(const mpfr::mpreal& velocity)
{
	casema::ModelData<mpfr::mpreal> data;

	data.nComponents = 1;
	
	data.initialLiquidConcentration.push_back(mpfr::mpreal(0));
	data.initialSolidConcentration.push_back(mpfr::mpreal(0));

	data.colDispersion = mpfr::mpreal("0.002");
	data.colLength = mpfr::mpreal("1.7");
	data.colPorosity = mpfr::mpreal("0.4");

	data.parRadius = mpfr::mpreal("0.004");
	data.parPorosity = mpfr::mpreal("0.333");

	data.filmDiffusion.push_back(mpfr::mpreal("0.01"));
	data.particleDiffusion.push_back(mpfr::mpreal("1e-6") / data.parPorosity);
	data.surfaceDiffusion.push_back(mpfr::mpreal(0));

	data.velocity = velocity;

	data.kineticBinding = false;

	data.bindingModel = casema::LINEAR;
	data.linearKA.push_back(mpfr::mpreal(25));
	data.linearKD.push_back(mpfr::mpreal(10));

	data.nInletSections = 2;
	data.sectionTimes.push_back(mpfr::mpreal(0));
	data.sectionTimes.push_back(mpfr::mpreal(20));
	data.sectionTimes.push_back(mpfr::mpreal(10000));
	data.constCoeff.push_back(mpfr::mpreal("0.001")); // 1 g / l in mol / cm^3
	data.constCoeff.push_back(mpfr::mpreal(0));

	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));
	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));

	return data;
}

void writeMoments(const casema::ModelData<mpfr::mpreal>& model)
{
	std::cout << "=========== MOMENTS (uncentered) ============" << std::endl; 

	const std::vector<mpfr::mpreal> moments = casema::analyticalMoments(model);
	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		std::cout << "MU_" << i << " = " << moments[i] << std::endl;
	}

	std::cout << "\n";


	std::cout << "=========== MOMENTS (central) ============" << std::endl; 

	const std::vector<mpfr::mpreal> momCent = casema::analyticalCentralMoments(model);
	for (std::size_t i = 0; i < momCent.size(); ++i)
	{
		std::cout << "mu_" << i << " = " << momCent[i] << std::endl;
	}

	std::cout << "\n";

	std::cout << "===== Check (central from uncentral) =====" << std::endl; 

	std::cout << "Mom " << 0 << ": " << abs(momCent[0] - moments[0]) << std::endl;
	std::cout << "Mom " << 1 << ": " << abs(momCent[1] - moments[1]) << std::endl;
	std::cout << "Mom " << 2 << ": " << abs(momCent[2] - moments[2] + moments[1] * moments[1]) << std::endl;
	std::cout << "Mom " << 3 << ": " << abs(momCent[3] - moments[3] + moments[1] * ( mpfr::mpreal(3) * moments[2] - mpfr::mpreal(2) * moments[1] * moments[1])) << std::endl;
	std::cout << "Mom " << 4 << ": " << abs(momCent[4] - moments[4] + moments[1] * ( mpfr::mpreal(4) * moments[3] - moments[1] * ( mpfr::mpreal(6) * moments[2] - mpfr::mpreal(3) * moments[1] * moments[1]))) << std::endl;

	std::cout << "\n";
}

int main(int argc, char** argv)
{

	std::size_t precision;
	double velocity;

	try
	{
		TCLAP::CmdLine cmd("Uses formulas from Qamar et al. (2014) to compute moments of models from the same publication", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<double>("v", "velocity", "Interstitial velocity in cm/min (default: 0.4 cm/min)", false, 0.4, "Double"))->storeIn(&velocity);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	std::cout.flags(std::ios::scientific);
	std::cout.precision(precision);

	std::cout << "#################### SI units #######################" << std::endl;
	writeMoments(qamarModel(velocity / mpfr::mpreal(6000)));

	std::cout << "\n################# cm, min units #####################" << std::endl;
	writeMoments(qamarModelCmMin(velocity));

	::mpfr_free_cache();

	return 0;
}

