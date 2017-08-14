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

#include <omp.h>

#include "MPReal.hpp"

#include <tclap/CmdLine.h>
#include "CliUtil.hpp"

#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "Constants.hpp"
#include "LaplaceSolution.hpp"
#include "LaplaceInlet.hpp"
#include "CliErrorParser.hpp"

typedef casema::LaplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


template <typename real_t>
struct ProgramOptions
{
	std::size_t numThreads;
	std::size_t summands;
	real_t sigma;
	real_t errorTrunc;
	real_t errorCons;
	real_t error;
	real_t errorWeight;
};

template <class inletFunc_t, typename real_t, typename complex_t>
std::vector<complex_t> fourierCoeff(std::size_t precision, std::size_t summands, const real_t& tMax, const real_t& sigma, const inletFunc_t& f)
{
	std::vector<complex_t> funcEvals(summands + 1);

	const real_t T = tMax / real_t(2) * real_t("1.01");
	funcEvals[0] = f(complex_t(sigma)).real() / real_t(2);

	#pragma omp parallel
	{
		casema::Constants<real_t>::init(precision);
		const real_t zero(0);

		#pragma omp for schedule(static)
		for (std::size_t k = 1; k <= summands; ++k)
		{
			const real_t kpit = k * casema::Constants<real_t>::pi() / T;
			funcEvals[k] = f(sigma + complex_t(zero, kpit));
		}

		casema::Constants<real_t>::clear();
	}
	return funcEvals;
}


template <class inletFunc_t, typename real_t>
std::vector<real_t> timeDomainInlet(std::size_t precision, const real_t& tMax, const inletFunc_t& f)
{
	std::vector<real_t> funcEvals(static_cast<unsigned int>(tMax) + 1);

	#pragma omp parallel
	{
		casema::Constants<real_t>::init(precision);
		const real_t zero(0);

		#pragma omp for schedule(static)
		for (std::size_t k = 0; k < funcEvals.size(); ++k)
		{
			funcEvals[k] = f.timeDomain(real_t(k));
		}

		casema::Constants<real_t>::clear();
	}
	return funcEvals;
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpcomplex>& con, const std::vector<mpfr::mpreal>& td, std::size_t precision)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "k,real,imag,abs,logabs,logabsreal,logabsimag\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < con.size(); ++i)
	{
		fs << i << "," << con[i].real() << "," << con[i].imag() << "," << abs(con[i]) << "," << log10(abs(con[i])) << "," << log10(abs(con[i].real())) << "," << log10(abs(con[i].imag())) << "\n"; 
	}

	fs << std::endl;
	fs << "##### Time Domain #####\n";
	fs << "t,inlet\n";

	for (std::size_t i = 0; i < td.size(); ++i)
	{
		fs << i << "," << td[i] << "\n"; 
	}

	fs.precision(curPrec);
	fs.flags(curFlags);
}


void run(const casema::ModelData<mpfr::mpreal>& model, std::size_t summands, const mpfr::mpreal& sigma, std::size_t precision, std::size_t outPrecision, const std::string& outFile)
{
	const mpfr::mpreal tMax = casema::maxSimulationTime(model);

	Inlet_t inlet(model);
	const std::vector<mpfr::mpcomplex> output = fourierCoeff<Inlet_t, mpfr::mpreal, mpfr::mpcomplex>(precision, summands, tMax, sigma, inlet);
	const std::vector<mpfr::mpreal> timeDomain = timeDomainInlet<Inlet_t, mpfr::mpreal>(precision, tMax, inlet);

	// Write infos
	const mpfr::mpreal T = tMax / mpfr::mpreal(2) * mpfr::mpreal("1.01");
	
	std::ios_base::fmtflags curFlags = std::cout.flags();
	std::streamsize curPrec = std::cout.precision();
	std::cout.flags(std::ios::scientific);
	std::cout.precision(precision);

	std::cout << "T = " << T << std::endl;
	std::cout << "pi / T = " << casema::Constants<mpfr::mpreal>::pi() / T << std::endl;
	std::cout << "sigma = " << sigma << std::endl;
	std::cout << "Summands = " << summands << std::endl;
	std::cout << "Max arg = " << casema::Constants<mpfr::mpreal>::pi() / T * summands << std::endl;
	std::cout << "Error factor = " << sqrt(mpfr::mpreal(2)) * exp(model.colLength * model.velocity / (2 * model.colDispersion)) << std::endl;
	std::cout << "Log10(Error factor) = " << log10(sqrt(mpfr::mpreal(2)) * exp(model.colLength * model.velocity / (2 * model.colDispersion))) << std::endl;
	std::cout << "Error exponent factor = " << model.colLength * sqrt(casema::Constants<mpfr::mpreal>::pi() / (2 * T * model.colDispersion)) << std::endl;

	std::cout.precision(curPrec);
	std::cout.flags(curFlags);
	
	if (outFile.empty())
	{
		writeResult(std::cout, output, timeDomain, outPrecision);
	}
	else
	{
		std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);
		writeResult(fs, output, timeDomain, outPrecision);
	}
}


int main(int argc, char** argv)
{
	std::size_t precision;
	std::size_t outPrecision;
	ProgramOptions<mpfr::mpreal> opts;
	std::string inFile;
	std::string outFile;

	try
	{
		TCLAP::CmdLine cmd("Evaluates Laplace transform and time domain of piecewise cubic polynomial inlet profile", ' ', "");

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: 16)", false, 16, "Int"))->storeIn(&outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method, used as safety margin if error (-e) is given", false, 0, "Float"))->storeIn(&opts.sigma);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Error threshold (default: 1e-10)", false, 1e-10, "Float"))->storeIn(&opts.error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Weight used to distribute error onto consistency and truncation (default: 0.5)", false, 0.5, "Float"))->storeIn(&opts.errorWeight);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "sum", "Maximum number of summands in Durbin's method", false, 0, "Int"))->storeIn(&opts.summands);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&outFile);
		cmd >> (new TCLAP::ValueArg<std::size_t>("t", "threads", "Number of threads (default: 4)", false, 4, "Int"))->storeIn(&opts.numThreads);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file", true, "", "File"))->storeIn(&inFile);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
	omp_set_num_threads(opts.numThreads);

	std::cout << "Precision: " << precision << " digits (base 10) = " << mpfr::digits2bits(precision) << " bit" << std::endl;
	
	casema::ModelData<mpfr::mpreal> model;
	bool loaded = casema::readModel<mpfr::mpreal>(inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		if (casema::checkModelForCompatibility(model))
		{
			if (!casema::processErrorOptions(opts, model))
			{
				std::cout << "ERROR: Valid parameter combinations are error (-e) and weight (-w) (abscissa (-a) is allowed), or abscissa (-a) and summands (-n)" << std::endl;
				::mpfr_free_cache();
				return 1;
			}

			std::cout << "Detected inlet: " << (casema::inletIsStep(model) ? "Step like" : "Bounded") << std::endl;
			if (model.kineticBinding)
			{
				std::cout << "WARNING: Detected kinetic binding; error estimates may not hold" << std::endl;
			}

			run(model, opts.summands, opts.sigma, precision, outPrecision, outFile);
		}
	}

	::mpfr_free_cache();

	return 0;
}
