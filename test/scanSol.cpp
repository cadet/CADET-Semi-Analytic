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

#include "MPReal.hpp"

#include <tclap/CmdLine.h>
#include "CliUtil.hpp"

#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

typedef casema::laplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


inline void progressBar(unsigned int cur, unsigned int total, unsigned int w = 50)
{
	if ( (cur != total) && (cur % (total/100+1) != 0) ) return;
 
	double ratio   =  cur / static_cast<double>(total);
	unsigned int c =  ratio * w;
 
	std::cout << std::setw(3) << static_cast<int>(ratio * 100.0) << "% [";
	for (unsigned int i = 0; i < c; ++i) 
		std::cout << "=";
	for (unsigned int i = c; i < w; ++i) 
		std::cout << " ";
	std::cout << "]\r" << std::flush;
}

template <class lapFunc_t, typename real_t, typename complex_t>
void writeResult(std::ostream& fs, std::size_t precision, const lapFunc_t& f, std::size_t nRe, std::size_t nIm, const real_t& maxRe, const real_t& minRe, const real_t& maxIm, const real_t& minIm, bool useLogLinear, std::size_t outPrecision, bool progBar)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "x,y,real,imag,abs,logabs,logabsreal,logabsimag\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	const std::size_t total = nRe * nIm;
	std::size_t cur = 0;

	real_t dx = (nRe > 1) ? (maxRe - minRe) / (nRe - 1) : real_t(0);
	real_t dy = (nIm > 1) ? (maxIm - minIm) / (nIm - 1) : real_t(0);
	real_t startRe = minRe;
	real_t startIm = minIm;

	if (useLogLinear)
	{
		dx = (nRe > 1) ? (log10(maxRe) - log10(minRe)) / (nRe - 1) : real_t(0);
		dy = (nIm > 1) ? (log10(maxIm) - log10(minIm)) / (nIm - 1) : real_t(0);
		if ((nRe > 1) || ((nRe == 1) && (minRe != 0)))
			startRe = log10(minRe);
		if ((nIm > 1) || ((nIm == 1) && (minIm != 0)))
			startIm = log10(minIm);
	}

	for (std::size_t i = 0; i < nRe; ++i)
	{
		for (std::size_t j = 0; j < nIm; ++j)
		{
			const real_t xv = startRe + i * dx;
			const real_t yv = startIm + j * dy;
			
			complex_t fv;
			if (useLogLinear)
			{
				const real_t curX = ((nRe > 1) || ((nRe == 1) && (minRe != 0))) ? pow(10, xv) : xv;
				const real_t curY = ((nIm > 1) || ((nIm == 1) && (minIm != 0))) ? pow(10, yv) : yv;
				fv = f(complex_t(curX, curY));
				fs << curX << "," << curY << ",";
			}
			else
			{
				fv = f(complex_t(xv, yv));
				fs << xv << "," << yv << ",";
			}

			fs << fv.real() << "," << fv.imag() << "," << abs(fv) << "," << log10(abs(fv)) << "," << log10(abs(fv.real())) << "," << log10(abs(fv.imag())) << "\n"; 
			
			++cur;
			if (progBar)
				progressBar(cur, total);
		}
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);

	std::cout << std::endl;
}


void run(casema::ModelData<mpfr::mpreal>& model, std::size_t nRe, std::size_t nIm, const mpfr::mpreal& maxRe, const mpfr::mpreal& minRe, const mpfr::mpreal& maxIm, const mpfr::mpreal& minIm, bool useLogLinear, std::size_t precision, std::size_t outPrecision, const std::string& outFile)
{
	std::cout << "Computing " << nRe * nIm << " function evaluations" << std::endl;

	Inlet_t inlet(model);
	if (model.kineticBinding)
	{
		typedef casema::laplaceSolution::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		if (outFile.empty())
		{
			writeResult<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(std::cout, outPrecision, solution, nRe, nIm, maxRe, minRe, maxIm, minIm, useLogLinear, precision, false);
		}
		else
		{
			std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);
			writeResult<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(fs, outPrecision, solution, nRe, nIm, maxRe, minRe, maxIm, minIm, useLogLinear, precision, true);
		}
	}
	else
	{
		typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		if (outFile.empty())
		{
			writeResult<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(std::cout, outPrecision, solution, nRe, nIm, maxRe, minRe, maxIm, minIm, useLogLinear, precision, false);
		}
		else
		{
			std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);
			writeResult<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(fs, outPrecision, solution, nRe, nIm, maxRe, minRe, maxIm, minIm, useLogLinear, precision, true);
		}
	}
}


int main(int argc, char** argv)
{
	std::size_t precision;
	std::size_t outPrecision;
	mpfr::mpreal maxRe;
	mpfr::mpreal minRe;
	mpfr::mpreal maxIm;
	mpfr::mpreal minIm;
	std::size_t nRe;
	std::size_t nIm;
	std::string inFile;
	std::string outFile;
	bool useLogLinear;

	try
	{
		TCLAP::CmdLine cmd("Scans the Laplace transform of GRM models", ' ', "");

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: 16)", false, 16, "Int"))->storeIn(&outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "maxre", "Maximum real value (default: 5.0)", false, 5.0, "Float"))->storeIn(&maxRe);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "maxim", "Maximum imaginary value (default: 1.0e3)", false, 1000.0, "Float"))->storeIn(&maxIm);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "minre", "Minimum real value (default: 0.0)", false, 0.0, "Float"))->storeIn(&minRe);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "minim", "Minimum imaginary value (default: -1.0e-3)", false, -1000.0, "Float"))->storeIn(&minIm);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "nre", "Number of real grid points (default: 100)", false, 100, "Int"))->storeIn(&nRe);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "nim", "Number of imaginary grid points (default: 20000)", false, 20001, "Int"))->storeIn(&nIm);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&outFile);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file", true, "", "File"))->storeIn(&inFile);
		cmd >> (new TCLAP::SwitchArg("l", "loglinear", "Use log-linear spacing"))->storeIn(&useLogLinear);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	std::cout << "Precision: " << precision << " digits (base 10) = " << mpfr::digits2bits(precision) << " bit" << std::endl;

	casema::ModelData<mpfr::mpreal> model;
	bool loaded = casema::readModel<mpfr::mpreal>(inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		if (casema::checkModelForCompatibility(model))        
			run(model, nRe, nIm, maxRe, minRe, maxIm, minIm, useLogLinear, precision, outPrecision, outFile);
	}

	::mpfr_free_cache();

	return 0;
}
