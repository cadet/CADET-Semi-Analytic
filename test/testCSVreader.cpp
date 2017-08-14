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
#include "CSVReader.hpp"


int main(int argc, char** argv)
{

	std::string inFile;
	size_t precision;

	try
	{
		TCLAP::CmdLine cmd("Reads the given CSV file", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("file", "CSV file", true, "", "File"))->storeIn(&inFile);
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

	casema::CSVReader reader;
	if (reader.readFile(inFile))
	{
		const std::size_t len = reader.numColumns();
		std::cout << "Columns: " << len << std::endl;

		std::cout << "Header: ";
		for (std::size_t i = 0; i < len - 1; ++i)
		{
			std::cout << reader.header(i) << ", ";
		}
		std::cout << reader.header(len-1) << std::endl;

		std::cout << "Length: ";
		for (std::size_t i = 0; i < len - 1; ++i)
		{
			std::cout << reader.column(i).size() << ", ";
		}
		std::cout << reader.column(len-1).size() << std::endl;

		for (std::size_t j = 0; j < std::min(std::size_t(10), reader.column(0).size()); ++j)
		{
			for (std::size_t i = 0; i < len - 1; ++i)
			{
				std::cout << reader.column(i)[j] << ", ";
			}
			std::cout << reader.column(len-1)[j] << std::endl;
		}
	}


	::mpfr_free_cache();

	return 0;
}

