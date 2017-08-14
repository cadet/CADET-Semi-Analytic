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
#include "CliUtil.hpp"

#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "Constants.hpp"
#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

typedef casema::laplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


template <typename real_t>
struct ProgramOptions
{
	std::size_t precision;
	std::size_t outPrecision;
	std::size_t summands;
	std::size_t skipSummands;
	mpfr::mpreal sigma;
	mpfr::mpreal tMax;
	std::string inFile;
	std::string outFile;
	bool withoutInlet;
};


template <class lapFunc_t, typename real_t, typename complex_t>
void fourierCoeff(std::size_t precision, std::size_t summands, const real_t& tMax, const real_t& sigma, const lapFunc_t& f, const real_t& errorExponentFactor, std::vector<complex_t>& funcEvals, std::vector<real_t>& errorEval)
{
	const real_t T = tMax / real_t(2) * real_t("1.01");

	#pragma omp parallel
	{
		casema::Constants<real_t>::init(precision);
		const real_t zero(0);

		#pragma omp for schedule(static)
		for (std::size_t k = 1; k <= summands; ++k)
		{
			const real_t kpit = k * casema::Constants<real_t>::pi() / T;
			funcEvals[k-1] = f.withoutInlet(sigma + complex_t(zero, kpit));
			errorEval[k-1] -= sqrt(real_t(k)) * errorExponentFactor;
		}

		casema::Constants<real_t>::clear();
	}
}


template <class lapFunc_t, typename real_t, typename complex_t>
void fourierCoeff(std::size_t precision, std::size_t summands, const real_t& tMax, const real_t& sigma, const lapFunc_t& f, const Inlet_t& in, const real_t& errorExponentFactor, std::vector<complex_t>& funcEvals, std::vector<real_t>& errorEval)
{
	const real_t T = tMax / real_t(2) * real_t("1.01");

	#pragma omp parallel
	{
		casema::Constants<real_t>::init(precision);
		const real_t zero(0);

		#pragma omp for schedule(static)
		for (std::size_t k = 1; k <= summands; ++k)
		{
			const real_t kpit = k * casema::Constants<real_t>::pi() / T;
			funcEvals[k-1] = f(sigma + complex_t(zero, kpit));
			errorEval[k-1] += log10(abs(in(sigma + complex_t(zero, kpit)))) - sqrt(real_t(k)) * errorExponentFactor;
		}

		casema::Constants<real_t>::clear();
	}
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpcomplex>& con, const std::vector<mpfr::mpreal>& error, std::size_t precision, std::size_t skipSummands)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "k,abs,logabs,error,logerror\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < con.size(); i += skipSummands)
	{
		fs << (i+1) << "," << abs(con[i]) << "," << log10(abs(con[i])) << "," << exp10(error[i]) << "," << error[i] << "\n"; 
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


void writeMeta(std::ostream& os, bool kinetic, ProgramOptions<mpfr::mpreal>& opts)
{
	os << "# Model: " << opts.inFile << std::endl;
	os << "# Binding: " << (kinetic ? "Dynamic" : "Rapid-Equlibrium") << "\n";
	os << "# Precision " << opts.precision << " digits = " << mpfr::digits2bits(opts.precision) << " bits\n";
	os << "# Inlet included: " << (opts.withoutInlet ? "No" : "Yes") << "\n";

	os << "# Summands: " << opts.summands << "\n";
	os << "# Abscissa (a): " << opts.sigma << "\n";
	os << "# Max time T: " << opts.tMax << std::endl;
}


void run(const casema::ModelData<mpfr::mpreal>& model, ProgramOptions<mpfr::mpreal>& opts)
{
	// Log10(errorPrefactor)
	const mpfr::mpreal errorPrefactor = log10(sqrt(mpfr::mpreal(2))) + model.colLength * model.velocity / (2 * model.colDispersion * log(mpfr::mpreal(10)));
	const mpfr::mpreal T = opts.tMax / mpfr::mpreal(2) * mpfr::mpreal("1.01");
	const mpfr::mpreal errorExponentFactor = model.colLength * sqrt(casema::Constants<mpfr::mpreal>::pi() / (2 * T * model.colDispersion)) / log(mpfr::mpreal(10));

	std::vector<mpfr::mpcomplex> output(opts.summands);
	std::vector<mpfr::mpreal> error(opts.summands, errorPrefactor);
	Inlet_t inlet(model);

	if (model.kineticBinding)
	{
		typedef casema::laplaceSolution::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		if (opts.withoutInlet)
			fourierCoeff<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(opts.precision, opts.summands, opts.tMax, opts.sigma, solution, errorExponentFactor, output, error);
		else
			fourierCoeff<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(opts.precision, opts.summands, opts.tMax, opts.sigma, solution, inlet, errorExponentFactor, output, error);
	}
	else
	{
		typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		if (opts.withoutInlet)
			fourierCoeff<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(opts.precision, opts.summands, opts.tMax, opts.sigma, solution, errorExponentFactor, output, error);
		else
			fourierCoeff<Solution_t, mpfr::mpreal, mpfr::mpcomplex>(opts.precision, opts.summands, opts.tMax, opts.sigma, solution, inlet, errorExponentFactor, output, error);
	}

	if (opts.outFile.empty())
	{
		writeMeta(std::cout, model.kineticBinding, opts);
		writeResult(std::cout, output, error, opts.outPrecision, opts.skipSummands);
	}
	else
	{
		std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);
		writeMeta(fs, model.kineticBinding, opts);
		writeResult(fs, output, error, opts.outPrecision, opts.skipSummands);
	}
}


int main(int argc, char** argv)
{
	ProgramOptions<mpfr::mpreal> opts;

	try
	{
		TCLAP::CmdLine cmd("Compares Fourier coefficients to their error estimate", ' ', "");

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: 16)", false, 16, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method (default: 1e-2)", false, 1e-2, "Float"))->storeIn(&opts.sigma);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "sum", "Number of summands in Durbin's method (default: 10000)", false, 10000, "Int"))->storeIn(&opts.summands);
		cmd >> (new TCLAP::ValueArg<std::size_t>("s", "skip", "Number of summands to skip in output (default: 10)", false, 10, "Int"))->storeIn(&opts.skipSummands);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("T", "tmax", "End time (default: taken from model)", false, -1, "Float"))->storeIn(&opts.tMax);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&opts.outFile);
		cmd >> (new TCLAP::SwitchArg("i", "noinlet", "Exclude inlet from Laplace domain solution"))->storeIn(&opts.withoutInlet);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file", true, "", "File"))->storeIn(&opts.inFile);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(opts.precision));

	std::cout << "Precision: " << opts.precision << " digits (base 10) = " << mpfr::digits2bits(opts.precision) << " bit" << std::endl;
	
	casema::ModelData<mpfr::mpreal> model;
	bool loaded = casema::readModel<mpfr::mpreal>(opts.inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		if (casema::checkModelForCompatibility(model))
		{
			std::cout << "Detected inlet: " << (casema::inletIsStep(model) ? "Step like" : "Bounded") << std::endl;
			if (model.kineticBinding)
			{
				std::cout << "WARNING: Detected kinetic binding; error estimates may not hold" << std::endl;
			}

			if (opts.tMax <= 0)
				opts.tMax = casema::maxSimulationTime(model);

			run(model, opts);
		}
	}

	::mpfr_free_cache();

	return 0;
}
