// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
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
#include "CenteredMoments.hpp"


mpfr::mpreal testCenteredFromNonCentered(std::size_t rounds)
{
	mpfr::mpreal maxDeviation(0);
	for (std::size_t i = 0; i < rounds; ++i)
	{
		// Construct random vector
		std::vector<mpfr::mpreal> data(5, mpfr::mpreal(0));

		for (int k = 0; k < 5; ++k)
		{
			// Random numbers between 0 and 100
			data[k] = mpfr::random() * mpfr::mpreal(100);
		}

		const std::vector<mpfr::mpreal> cent = casema::centralMomentsFromNonCentral(data);

		maxDeviation = max(abs(cent[0] - data[0]), maxDeviation);
		maxDeviation = max(abs(cent[1] - data[1]), maxDeviation);
		maxDeviation = max(abs(cent[2] - data[2] + data[1] * data[1]), maxDeviation);
		maxDeviation = max(abs(cent[3] - data[3] + data[1] * ( mpfr::mpreal(3) * data[2] - mpfr::mpreal(2) * data[1] * data[1])), maxDeviation);
		maxDeviation = max(abs(cent[4] - data[4] + data[1] * ( mpfr::mpreal(4) * data[3] - data[1] * ( mpfr::mpreal(6) * data[2] - mpfr::mpreal(3) * data[1] * data[1]))), maxDeviation);
	}
	return maxDeviation;
}

int main(int argc, char** argv)
{

	std::size_t precision;
	std::size_t rounds;

	try
	{
		TCLAP::CmdLine cmd("Tests the computation of centered moments from non-centered ones", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "rounds", "Number of test rounds (default: 10)", false, 10, "Int"))->storeIn(&rounds);
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

	mpfr::mpreal maxDev = testCenteredFromNonCentered(rounds);
	std::cout << "Maximum deviation: " << maxDev << std::endl;

	::mpfr_free_cache();

	return 0;
}

