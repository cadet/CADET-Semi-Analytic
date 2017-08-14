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

#include "CliUtil.hpp"
#include "ModelData.hpp"
#include "ErrorEstimate.hpp"


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
	data.constCoeff.push_back(mpfr::mpreal(1)); // 1 g / l in mol / m^3
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
	data.constCoeff.push_back(mpfr::mpreal("0.000001")); // 1 g / l in mol / cm^3
	data.constCoeff.push_back(mpfr::mpreal(0));

	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));
	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));

	return data;
}


void computeParamsFromError(const casema::ModelData<mpfr::mpreal>& model, const mpfr::mpreal& error, const mpfr::mpreal& weight)
{
	std::cout << "Desired error: " << error << std::endl;
	std::cout << "       Weight: " << weight << std::endl;

	std::cout << "Max inlet: " << casema::detail::maxInletSpline(model) << std::endl;;

	// Compute parameters    
	{
		unsigned int nSummands = 0;
		mpfr::mpreal abscissa;
		casema::paramsFromErrorStep(mpfr::mpreal(0), model, weight, error, abscissa, nSummands);
		std::cout << "== Step ==============" << std::endl;
		std::cout << "Abscissa: " << abscissa << std::endl;
		std::cout << "Summands: " << nSummands << std::endl;

		std::cout << "== Step ==============" << std::endl;

		const mpfr::mpreal ce = casema::consistencyError(abscissa, model);
		std::cout << "Consistency: " << ce << std::endl;

		const mpfr::mpreal te = casema::truncationErrorStep(abscissa, nSummands, model);
		std::cout << "Truncation : " << te << std::endl;

		std::cout << "Sum: " << te + ce << std::endl;    
	}
	{
		unsigned int nSummands = 0;
		mpfr::mpreal abscissa;
		casema::paramsFromErrorBounded(mpfr::mpreal(0), model, weight, error, abscissa, nSummands);
		std::cout << "== Bounded ===========" << std::endl;
		std::cout << "Abscissa: " << abscissa << std::endl;
		std::cout << "Summands: " << nSummands << std::endl;

		std::cout << "== Bounded ===========" << std::endl;

		const mpfr::mpreal ce = casema::consistencyError(abscissa, model);
		std::cout << "Consistency: " << ce << std::endl;

		const mpfr::mpreal te = casema::truncationErrorBounded(abscissa, nSummands, model);
		std::cout << "Truncation : " << te << std::endl;

		std::cout << "Sum: " << te + ce << std::endl;    
	}
}


void computeError(const casema::ModelData<mpfr::mpreal>& model, const mpfr::mpreal& abscissa, unsigned int nSummands)
{
	std::cout << "Abscissa: " << abscissa << std::endl;
	std::cout << "Summands: " << nSummands << std::endl;
	std::cout << "Max inlet: " << casema::detail::maxInletSpline(model) << std::endl;;

	std::cout << "== Step ==============" << std::endl;

	{
		const mpfr::mpreal ce = casema::consistencyError(abscissa, model);
		std::cout << "Consistency: " << ce << std::endl;

		const mpfr::mpreal te = casema::truncationErrorStep(abscissa, nSummands, model);
		std::cout << "Truncation : " << te << std::endl;

		std::cout << "Sum: " << te + ce << std::endl;    
	}

	std::cout << "== Bounded ===========" << std::endl;

	{
		const mpfr::mpreal ce = casema::consistencyError(abscissa, model);
		std::cout << "Consistency: " << ce << std::endl;

		const mpfr::mpreal te = casema::truncationErrorBounded(abscissa, nSummands, model);
		std::cout << "Truncation : " << te << std::endl;

		std::cout << "Sum: " << te + ce << std::endl;    
	}
}


int main(int argc, char** argv)
{
	std::size_t precision;
	double velocity;
	mpfr::mpreal error;
	mpfr::mpreal weight;
	mpfr::mpreal abscissa;
	unsigned int nSummands;

	try
	{
		TCLAP::CmdLine cmd("Tests error analysis code", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<double>("v", "velocity", "Interstitial velocity in cm/min (default: 0.5 cm/min)", false, 0.5, "Double"))->storeIn(&velocity);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Desired error", false, 0.0, "Double"))->storeIn(&error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Distribution of error onto consistency and truncation (default: 0.5)", false, mpfr::mpreal(0.5), "Double"))->storeIn(&weight);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa", false, 0.0, "Double"))->storeIn(&abscissa);
		cmd >> (new TCLAP::ValueArg<unsigned int>("n", "summands", "Number of summands", false, 0, "Int"))->storeIn(&nSummands);

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

	const casema::ModelData<mpfr::mpreal> model = qamarModel(velocity / mpfr::mpreal(6000));

	if (error > 0)
		computeParamsFromError(model, error, weight);
	else if ((nSummands > 0) && (abscissa > 0))
		computeError(model, abscissa, nSummands);
	else
		std::cout << "ERROR: Either error (-e) or abscissa (-a) and number of summands (-n) required" << std::endl;

	::mpfr_free_cache();

	return 0;
}

