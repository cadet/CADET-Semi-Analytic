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
#include "CliExtrapolationParser.hpp"

#include "TclapCustomOutput.hpp"
#include "VersionInfo.hpp"
#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "DurbinsMethod.hpp"
#include "LaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"

#include "CliErrorParser.hpp"

typedef casema::LaplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;

template <typename real_t>
struct ProgramOptions
{
	std::size_t summands;
	std::size_t precision;
	std::size_t outPrecision;
	real_t sigma;
	real_t errorTrunc;
	real_t errorCons;
	real_t error;
	real_t errorWeight;
	std::string inFile;
	std::string outFile;
	std::size_t numThreads;

	std::string extraMethods;
	std::size_t convTimes;
	real_t convThresholdAbs;
	real_t convThresholdRel;
	std::size_t consAgree;
};


template <typename real_t>
void writeMeta(std::ostream& os, const casema::ModelData<real_t>& model, casema::ConsensusEstimator<real_t> const* const estimator, const ProgramOptions<real_t>& opts)
{
	os << "# Version: " << casema::getVersion() << "\n";
	os << "# Commit: " << casema::getCommitHash() << "\n";
#ifdef CASEMA_KAHAN_SUMMATION
	os << "# Kahan summation: On\n";
#else
	os << "# Kahan summation: Off\n";
#endif
	os << "# Model: " << opts.inFile << "\n";
	os << "# Binding: " << (model.kineticBinding ? "Dynamic" : "Rapid-Equlibrium") << "\n";
	os << "# Precision " << opts.precision << " digits = " << mpfr::digits2bits(opts.precision) << " bits\n";

	if (estimator)
		os << "# Extrapolation: " << estimator->name() << "\n";
	else
		os << "# Extrapolation: Disabled\n";

	os << "# Convergence: " << opts.convTimes << " times below " << opts.convThresholdAbs << " (abs) or " << opts.convThresholdRel << " (rel) \n";
	os << "# Consensus: " << opts.consAgree << " extrapolators\n";
	os << "# MaxSummands: " << opts.summands << "\n";
	os << "# Abscissa (a): " << opts.sigma << "\n";
	os << "# Truncation error: " << opts.errorTrunc << "\n";
	os << "# Consistency error: " << opts.errorCons << "\n";
	os << "# Error: " << opts.error << "\n";
	os << "# Error weight: " << opts.errorWeight << std::endl;
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpreal>& time, const std::vector<mpfr::mpreal>& con, std::size_t precision, std::size_t timeOffset)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "time,outlet\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < con.size(); ++i)
	{
		fs << time[i+timeOffset] << "," << con[i] << "\n"; 
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpreal>& time, const std::vector<mpfr::mpreal>& con, const std::vector<mpfr::mpreal>& radius, std::size_t precision, std::size_t timeOffset)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "time,outlet,radius\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < con.size(); ++i)
	{
		fs << time[i+timeOffset] << "," << con[i] << "," << radius[i] << "\n";
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


template <class Model_t, class Solution_t, class Inlet_t>
std::vector<mpfr::mpreal> invert(const Model_t& model, const Inlet_t& inlet, const std::size_t timeOffset, const ProgramOptions<mpfr::mpreal>& opts)
{
	Solution_t solution(model, inlet);
	std::cout << solution.name() << std::endl;
	casema::DurbinsMethod<mpfr::mpreal, mpfr::mpcomplex> durbin;
	return durbin.invert(opts.precision, opts.summands, casema::maxSimulationTime(model), opts.sigma, model.outletTimes, timeOffset, solution);
}


void run(casema::ModelData<mpfr::mpreal>& model, const ProgramOptions<mpfr::mpreal>& opts, casema::ConsensusEstimator<mpfr::mpreal>& estimator, bool useExtrapolation)
{
	if (!model.writeUserTimes)
	{
		std::cout << "Generating output time points" << std::endl;
		model.outletTimes.reserve(model.sectionTimes[model.sectionTimes.size()-1].toLong() + 1);

		mpfr::mpreal cur = mpfr::mpreal(1);
		while (cur <= model.sectionTimes[model.sectionTimes.size()-1])
		{
			model.outletTimes.push_back(cur);
			cur += mpfr::mpreal(1);
		}
	}

	std::size_t timeOffset = 0;
	if (model.outletTimes[0] <= mpfr::mpreal(0))
	{
		// Skip first time point
		std::cout << "Removed time point t = 0.0 since Laplace solution has a pole there" << std::endl;
		timeOffset = 1;
	}

	std::cout << "Starting Laplace inversion for Model" << std::endl;

	std::vector<mpfr::mpreal> output;
	Inlet_t inlet(model);
	switch (model.modelType)
	{
		case casema::GeneralRateModel:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::GeneralRateModel::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			else
			{
				typedef casema::LaplaceSolution::GeneralRateModel::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			break;
		case casema::LumpedRateModelWithPores:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithPores::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			else
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithPores::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			break;
		case casema::LumpedRateModelWithoutPores:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithoutPores::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			else
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithoutPores::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				output = invert<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, timeOffset, opts);
			}
			break;
	}

	std::cout << "Writing output" << std::endl;
	if (opts.outFile.empty())
	{
		writeResult(std::cout, model.outletTimes, output, opts.outPrecision, timeOffset);
	}
	else
	{
		std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);
		
		if (useExtrapolation)
			writeMeta<mpfr::mpreal>(fs, model, nullptr, opts);
		else
			writeMeta<mpfr::mpreal>(fs, model, nullptr, opts);

		writeResult(fs, model.outletTimes, output, opts.outPrecision, timeOffset);
	}
}


int main(int argc, char** argv)
{
	ProgramOptions<mpfr::mpreal> opts;

	try
	{
		TCLAP::CustomOutput customOut("chrom");
		TCLAP::CmdLine cmd("Uses a numerical inverse Laplace transform to solve GRM models", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<std::size_t>("t", "threads", "Number of threads (default: 4)", false, 4, "Int"))->storeIn(&opts.numThreads);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method, used as safety margin if error (-e) is given", false, 0, "Float"))->storeIn(&opts.sigma);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Error threshold (default: 1e-10)", false, 1e-10, "Float"))->storeIn(&opts.error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Weight used to distribute error onto consistency and truncation (default: 0.5)", false, 0.5, "Float"))->storeIn(&opts.errorWeight);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "sum", "Maximum number of summands in Durbin's method", false, 0, "Int"))->storeIn(&opts.summands);

		cmd >> (new TCLAP::ValueArg<std::string>("x", "extrapolation", "Use extrapolation methods (default: none)", false, std::string(), "Method1,Method2"))->storeIn(&opts.extraMethods);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "convtimes", "Number of successful convergence tests (default: 3)", false, 3, "Int"))->storeIn(&opts.convTimes);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "convabs", "Absolute threshold used in convergence test (default: 100 * MachineEps)", false, mpfr::mpreal(-1), "Double"))->storeIn(&opts.convThresholdAbs);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "convrel", "Relative threshold used in convergence test (default: 0)", false, mpfr::mpreal(0), "Double"))->storeIn(&opts.convThresholdRel);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "agree", "Number of extrapolators used for accuracy estimate (default: 1)", false, 1, "Int"))->storeIn(&opts.consAgree);

		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&opts.outFile);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file (HDF5 or XML)", true, "", "File"))->storeIn(&opts.inFile);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(opts.precision));
	omp_set_num_threads(opts.numThreads);

	// Set default convergence threshold
	if (opts.convThresholdAbs < 0)
		opts.convThresholdAbs = mpfr::mpreal(100) * std::numeric_limits<mpfr::mpreal>::epsilon();

	// Set default output precision
	if (opts.outPrecision == 0)
		opts.outPrecision = opts.precision;

	// Create extrapolation method
	casema::ConsensusEstimator<mpfr::mpreal> estimator;
	const bool extrapActive = casema::addExtrapolationMethods(opts.extraMethods, estimator, opts);

	casema::ModelData<mpfr::mpreal> model;
	const bool loaded = casema::readModel<mpfr::mpreal>(opts.inFile, model);

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

			if (!opts.outFile.empty())
			{
				if (extrapActive)
					writeMeta<mpfr::mpreal>(std::cout, model, &estimator, opts);
				else
					writeMeta<mpfr::mpreal>(std::cout, model, nullptr, opts);
			}

			run(model, opts, estimator, extrapActive);
		}
	}

	::mpfr_free_cache();

	return 0;
}

