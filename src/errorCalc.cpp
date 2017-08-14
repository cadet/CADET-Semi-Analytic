// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2017: Samuel Leweke¹
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

#include "MPReal.hpp"

#include <tclap/CmdLine.h>

#include "StringUtil.hpp"
#include "TclapCustomOutput.hpp"
#include "VersionInfo.hpp"

#include "CadetEnumeration.hpp"
#include "CSVReader.hpp"


void writeFile(const std::vector<std::string>& colHeader, const std::vector<std::vector<mpfr::mpreal> const*>& colPtrs, std::ostream& fs, std::size_t precision)
{
	// Set output flags
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs.flags(std::ios::scientific);
	fs.precision(precision);

	// Write CSV header
	for (std::size_t i = 0; i < colHeader.size() - 1; ++i)
	{
		fs << colHeader[i] << ",";
	}
	fs << colHeader[colHeader.size() - 1] << "\n";

	// Get maximum number of lines
	std::size_t lines = 0;
	for (std::vector<std::vector<mpfr::mpreal> const*>::const_iterator it = colPtrs.begin(); it != colPtrs.end(); ++it)
	{
		lines = std::max(lines, (*it)->size());
	}

	// Write data
	const std::vector<mpfr::mpreal>& ptrA = *colPtrs[0];
	const std::vector<mpfr::mpreal>& ptrB = *colPtrs[1];
	for (std::size_t i = 0; i < lines; ++i)
	{
		// Calculate error
		mpfr::mpreal error(0);
		mpfr::mpreal relerror(0);

		if ((ptrA.size() > i) && (ptrB.size() > i))
		{
			error = abs(ptrA[i] - ptrB[i]);
			relerror = error / ptrA[i];
		}
		else
		{
			if (ptrA.size() > i)
				error = abs(ptrA[i]);
			else if (ptrB.size() > i)
				error = abs(ptrB[i]);

			relerror = error;
		}

		fs << error << "," << log10(error) << "," << relerror << "," << log10(relerror);

		if (colPtrs.size() > 2)
		{
			for (std::size_t j = 2; j < colPtrs.size(); ++j)
			{
				const std::vector<mpfr::mpreal>& ptrCur = *colPtrs[j];
				if (ptrCur.size() > i)
					fs << "," << ptrCur[i];
				else
					fs << ",";
			}
		}

		fs << "\n";
	}

	fs.precision(curPrec);
	fs.flags(curFlags);
}


int main(int argc, char** argv)
{
	std::size_t precision;
	std::size_t outPrecision;
	std::vector<std::string> files;
	std::string outFile;
	std::string keepColumnStr;
	std::string errorColumnStr;

	try
	{
		TCLAP::CustomOutput customOut("errorCalc");
		TCLAP::CmdLine cmd("Calculates point-wise difference of arbitrary CSV file columns", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write output to file (default: disabled)", false, std::string(), "File"))->storeIn(&outFile);
		cmd >> (new TCLAP::UnlabeledMultiArg<std::string>("data", "Files", true, "", "File"))->storeIn(&files);

		cmd >> (new TCLAP::ValueArg<std::string>("k", "keep", "Keep columns (default: none)", false, std::string(), "Column1,Column2"))->storeIn(&keepColumnStr);
		cmd >> (new TCLAP::ValueArg<std::string>("e", "error", "Calculate error between these two columns", true, std::string(), "Column1,Column2"))->storeIn(&errorColumnStr);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	if (files.size() != 2)
	{
		std::cerr << "ERROR: Need exactly two files, " << files.size() << " given" << std::endl;
		return 1;
	}

	if (outPrecision == 0)
		outPrecision = precision;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	// Check error column names
	const std::vector<std::string> errorColumns = casema::util::split(errorColumnStr, ",");
	if (errorColumns.size() != 2)
	{
		std::cerr << "ERROR: Need exactly two columns in --error / -e argument, got " << errorColumns.size() << std::endl;
		return 1;
	}

	// Read files
	casema::CSVReader readerA;
	if (!readerA.readFile(files[0]))
	{
		std::cerr << "ERROR reading " << files[0] << std::endl;
		return 1;
	}

	casema::CSVReader readerB;
	if (!readerB.readFile(files[1]))
	{
		std::cerr << "ERROR reading " << files[1] << std::endl;
		return 1;
	}

	const std::vector<std::string> keepColumns = casema::util::split(keepColumnStr, ",");

	// Get pointers to error columns and kept columns
	std::vector<std::string> colHeader;
	std::vector<std::vector<mpfr::mpreal> const*> colPtrs;

	// First four items are error columns
	colHeader.push_back("error");
	colHeader.push_back("logerror");
	colHeader.push_back("relerror");
	colHeader.push_back("logrelerror");

	bool colNotFound = true;
	if (readerA.contains(errorColumns[0]))
	{
		if (readerB.contains(errorColumns[1]))
		{
			std::cout << "Column " << errorColumns[0] << " in file " << files[0] << std::endl;
			std::cout << "Column " << errorColumns[1] << " in file " << files[1] << std::endl;
			colNotFound = false;
			colPtrs.push_back(&readerA.column(errorColumns[0]));
			colPtrs.push_back(&readerB.column(errorColumns[1]));
		}
	}
	if (colNotFound && readerB.contains(errorColumns[0]))
	{
		if (readerA.contains(errorColumns[1]))
		{
			std::cout << "Column " << errorColumns[0] << " in file " << files[1] << std::endl;
			std::cout << "Column " << errorColumns[1] << " in file " << files[2] << std::endl;
			colNotFound = false;
			colPtrs.push_back(&readerA.column(errorColumns[1]));
			colPtrs.push_back(&readerB.column(errorColumns[0]));
		}
	}
	if (colNotFound)
	{
		std::cerr << "ERROR: Could not find columns " << errorColumns[0] << " and " << errorColumns[1] << " in files" << std::endl;
		return 1;
	}

	// Try to find kept columns
	for (std::vector<std::string>::const_iterator it = keepColumns.begin(); it != keepColumns.end(); ++it)
	{
		if (readerA.contains(*it))
		{
			colHeader.push_back(*it);
			colPtrs.push_back(&readerA.column(*it));
			std::cout << "Column " << *it << " in file " << files[0] << std::endl;
		}
		else if (readerB.contains(*it))
		{
			colHeader.push_back(*it);
			colPtrs.push_back(&readerB.column(*it));
			std::cout << "Column " << *it << " in file " << files[1] << std::endl;
		}
		else
			std::cout << "WARNING: Could not find column " << *it << std::endl;
	}

	std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);

	// Write file header
	fs << "# Sources: " << files[0] << " and " << files[1] << "\n";
	fs << "# Columns: " << errorColumns[0] << " (ref) and " << errorColumns[1] << std::endl;

	// Write data
	writeFile(colHeader, colPtrs, fs, outPrecision);

	::mpfr_free_cache();

	return 0;
}

