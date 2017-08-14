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
#include <quadpack/workspace.hpp>

#include <tclap/CmdLine.h>
#include "CliUtil.hpp"

#include "TclapCustomOutput.hpp"
#include "VersionInfo.hpp"
#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "DurbinsMethod.hpp"
#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "AnalyticalMoments.hpp"

#include "CliErrorParser.hpp"

typedef casema::laplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


class PrecisionGuard
{
public:
	PrecisionGuard(std::ostream& os, std::size_t precision) : _os(os)
	{
		_curFlags = os.flags();
		_curPrec = os.precision();

		os.flags(std::ios::scientific);
		os.precision(precision);
	}

	~PrecisionGuard()
	{
		_os.precision(_curPrec);
		_os.flags(_curFlags);
	}

private:
	std::ostream& _os;
	std::ios_base::fmtflags _curFlags;
	std::streamsize _curPrec;
};


template <typename real_t>
struct ProgramOptions
{
	std::size_t summands;
	std::size_t quadOrder;
	std::size_t quadIter;
	std::size_t moment;
	std::size_t precision;
	std::size_t outPrecision;
	real_t sigma;
	real_t errorTrunc;
	real_t errorCons;
	real_t error;
	real_t errorWeight;
	real_t mean;
	real_t quadRel;
	real_t quadAbs;
	real_t quadMin;
	real_t quadMax;
	std::string inFile;
	std::string outFile;
	std::size_t numThreads;
	bool central;
	bool nonCentral;

	bool noMeta;
};


template <typename real_t>
void writeMeta(std::ostream& os, const casema::ModelData<real_t>& model, const ProgramOptions<real_t>& opts)
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
	os << "# Precision " << opts.precision << " digits = " << mpfr::digits2bits(opts.precision) << " bits; Machine eps = " << std::numeric_limits<mpfr::mpreal>::epsilon() << "\n";
	os << "# Summands: " << opts.summands << "\n";
	os << "# Durbin abscissa (a): " << opts.sigma << "\n";
	os << "# Durbin truncation error: " << opts.errorTrunc << "\n";
	os << "# Durbin consistency error: " << opts.errorCons << "\n";
	os << "# Durbin error: " << opts.error << "\n";
	os << "# Durbin error weight: " << opts.errorWeight << std::endl;
	os << "# Quadrature: Order " << opts.quadOrder << " iterations " << opts.quadIter << "\n";
	os << "# Quad abs error: " << opts.quadAbs << "\n";
	os << "# Quad rel error: " << opts.quadRel << "\n";
	os << "# Quad range: " << opts.quadMin << " to " << opts.quadMax << "\n";

	if (!isnan(opts.mean))
		os << "# First moment: " << opts.mean << "\n";
}


template <typename real_t, typename Durbin_t>
struct FirstMomentFunctional
{
	FirstMomentFunctional(Durbin_t* const d) : durbin(d) { }

	real_t operator()(const real_t& t) const
	{
		return t * std::max(durbin->invert(t), real_t(0));
	}

	Durbin_t* durbin;
};


template <typename real_t, typename Durbin_t>
struct NonCentralMomentFunctional
{
	NonCentralMomentFunctional(Durbin_t* const d, std::size_t ord) : durbin(d), order(ord) { }

	real_t operator()(const real_t& t) const
	{
		return pow(t, order) * std::max(durbin->invert(t), real_t(0));
	}

	Durbin_t* durbin;
	std::size_t order;
};

template <typename real_t, typename Durbin_t>
struct CentralMomentFunctional
{
	CentralMomentFunctional(Durbin_t* const d, std::size_t ord, real_t* m) : durbin(d), order(ord), mean(m) { }

	real_t operator()(const real_t& t) const
	{
		return pow(t - (*mean), order) * std::max(durbin->invert(t), real_t(0));
	}

	Durbin_t* durbin;
	std::size_t order;
	real_t* mean;
};


template <typename real_t, class Solution_t>
void calcMoment(Solution_t& solution, const casema::ModelData<real_t>& model, const ProgramOptions<real_t>& opts, std::ostream& os)
{
	const real_t tMax = casema::maxSimulationTime(model);
	const real_t mass = casema::detail::injectedMass(model);
	const char* const header = "order,central,noncentral,errorcentral,errornoncentral\n";

	// Early out if 0th moment is wanted
	if (opts.moment == 0)
	{
		os << header << 0 << "," << mass << "," << mass << ",0,0" << std::endl;
		return;
	}

	// Precompute function evaluations in Durbin's method
	std::cout << "Precompute function evaluations for Durbin's method" << std::endl;
	casema::PrecomputedDurbinsMethod<real_t, mpfr::mpcomplex> durbin(opts.precision, opts.summands, tMax, opts.sigma, solution);

	// Construct quadrature rule
	std::cout << "Construct quadrature rule" << std::endl;
	QuadPack::Workspace<real_t> ws(opts.quadIter, opts.quadOrder);

	real_t mean = opts.mean;
	real_t meanError(0);
	const bool compMean = (opts.moment > 1) && opts.central && isnan(opts.mean);
	if (compMean || (opts.moment == 1))
	{
		FirstMomentFunctional<real_t, casema::PrecomputedDurbinsMethod<real_t, mpfr::mpcomplex>> func(&durbin);

		// Compute first moment
		std::cout << "Calculate first moment" << std::endl;
		const int errorCode = ws.qag(func, opts.quadMin, opts.quadMax, opts.quadAbs, opts.quadRel, mean, meanError);

		if ((errorCode != QuadPack::GSL_SUCCESS) && (errorCode != QuadPack::GSL_EMAXITER) && (errorCode != QuadPack::GSL_EROUND))
		{
			os << "ERROR: Integration (mean) failed with error code " << errorCode << ", " << QuadPack::returnCodeToString(errorCode) << std::endl;
			return;
		}

		mean /= mass;
	}

	// Early out if first moment (mean) is wanted
	if (opts.moment == 1)
	{
		os << header << 1 << "," << mean << "," << mean << "," << meanError << "," << meanError << std::endl;
		return;
	}

	// Compute n-th moment
	std::cout << "Compute " << opts.moment << "th moment" << std::endl;

	// Write mean if enabled
	if (!opts.noMeta && compMean)
	{
		os << "# Mean: " << mean << " error " << meanError << "\n";
	}
	os << header << opts.moment << ",";

	// Moment order > 1
	real_t momCent(0);
	real_t errorCent(0);
	if (opts.central)
	{
		CentralMomentFunctional<real_t, casema::PrecomputedDurbinsMethod<real_t, mpfr::mpcomplex>> func(&durbin, opts.moment, &mean);
		const int errorCode = ws.qag(func, opts.quadMin, opts.quadMax, opts.quadAbs, opts.quadRel, momCent, errorCent);
		momCent /= mass;

		if ((errorCode != QuadPack::GSL_SUCCESS) && (errorCode != QuadPack::GSL_EMAXITER) && (errorCode != QuadPack::GSL_EROUND))
		{
			os << "ERROR: Integration (central) failed with error code " << errorCode << ", " << QuadPack::returnCodeToString(errorCode) << std::endl;
			return;
		}
	}

	real_t momNonCent(0);
	real_t errorNonCent(0);
	if (opts.nonCentral)
	{
		NonCentralMomentFunctional<real_t, casema::PrecomputedDurbinsMethod<real_t, mpfr::mpcomplex>> func(&durbin, opts.moment);
		const int errorCode = ws.qag(func, opts.quadMin, opts.quadMax, opts.quadAbs, opts.quadRel, momNonCent, errorNonCent);
		momNonCent /= mass;

		if ((errorCode != QuadPack::GSL_SUCCESS) && (errorCode != QuadPack::GSL_EMAXITER) && (errorCode != QuadPack::GSL_EROUND))
		{
			os << "ERROR: Integration (noncentral) failed with error code " << errorCode << ", " << QuadPack::returnCodeToString(errorCode) << std::endl;
			return;
		}
	}

	os << momCent << "," << momNonCent << "," << errorCent << "," << errorNonCent << std::endl;
}


void run(const casema::ModelData<mpfr::mpreal>& model, const ProgramOptions<mpfr::mpreal>& opts)
{
	std::ostream* os = &std::cout;
	if (!opts.outFile.empty())
	{
		std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);
		if (!opts.noMeta)
			writeMeta(fs, model, opts);
	}

	PrecisionGuard guard(*os, opts.outPrecision);
	writeMeta(*os, model, opts);

	// Write meta to stdout too in case of file output
	if (!opts.outFile.empty())
	{
		PrecisionGuard guard2(std::cout, opts.outPrecision);
		writeMeta(std::cout, model, opts);
	}

	Inlet_t inlet(model);
	if (model.kineticBinding)
	{
		typedef casema::laplaceSolution::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		calcMoment<mpfr::mpreal, Solution_t>(solution, model, opts, *os);
	}
	else
	{
		typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		calcMoment<mpfr::mpreal, Solution_t>(solution, model, opts, *os);
	}

	if (!opts.outFile.empty())
	{
		delete os;
		os = &std::cout;
	}
}


int main(int argc, char** argv)
{
	ProgramOptions<mpfr::mpreal> opts;

	try
	{
		TCLAP::CustomOutput customOut("momentQuadrature");
		TCLAP::CmdLine cmd("Uses a numerical inverse Laplace transform and numerical quadrature to compute moments of GRM models", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<std::size_t>("t", "threads", "Number of threads (default: 4)", false, 4, "Int"))->storeIn(&opts.numThreads);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::SwitchArg("", "nometa", "Do not write meta info to file"))->storeIn(&opts.noMeta);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method, used as safety margin if error (-e) is given", false, 0, "Float"))->storeIn(&opts.sigma);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Error threshold (default: 1e-10)", false, 1e-10, "Float"))->storeIn(&opts.error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Weight used to distribute error onto consistency and truncation (default: 0.5)", false, 0.5, "Float"))->storeIn(&opts.errorWeight);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "sum", "Maximum number of summands in Durbin's method", false, 0, "Int"))->storeIn(&opts.summands);

		cmd >> (new TCLAP::ValueArg<std::size_t>("m", "order", "Order of Gauss-Kronrod rule (default: 10)", false, 10, "Int"))->storeIn(&opts.quadOrder);
		cmd >> (new TCLAP::ValueArg<std::size_t>("M", "iter", "Maximum iterations of adaptive quadrature (default: 25)", false, 25, "Int"))->storeIn(&opts.quadIter);
		cmd >> (new TCLAP::ValueArg<std::size_t>("N", "moment", "Order of moment (default: 1)", false, 1, "Int"))->storeIn(&opts.moment);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("r", "mean", "First moment (retention time)", false, mpfr::mpreal().setNan(), "Float"))->storeIn(&opts.mean);
		cmd >> (new TCLAP::SwitchArg("c", "central", "Compute central moment"))->storeIn(&opts.central);
		cmd >> (new TCLAP::SwitchArg("C", "noncentral", "Compute non-central moment"))->storeIn(&opts.nonCentral);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "max", "Upper integration bound (default: inferred from model)", false, mpfr::mpreal().setNan(), "Float"))->storeIn(&opts.quadMax);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "min", "Lower integration bound (default: 0)", false, mpfr::mpreal(0), "Float"))->storeIn(&opts.quadMin);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "quadabs", "Absolute quadrature error", false, mpfr::mpreal().setNan(), "Float"))->storeIn(&opts.quadAbs);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "quadrel", "Relative quadrature error", false, mpfr::mpreal().setNan(), "Float"))->storeIn(&opts.quadRel);

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

	// Set default output precision
	if (opts.outPrecision == 0)
		opts.outPrecision = opts.precision;

	// Compute both moments by default
	if (!(opts.nonCentral || opts.central))
	{
		opts.nonCentral = true;
		opts.central = true;
	}

	// Quadrature order should not be less than 3
	opts.quadOrder = std::max(opts.quadOrder, std::size_t(3));

	// Quadrature iterations should not be less than 5
	opts.quadIter = std::max(opts.quadIter, std::size_t(5));

	casema::ModelData<mpfr::mpreal> model;
	const bool loaded = casema::readModel<mpfr::mpreal>(opts.inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		// Set default quadrature range
		if (isnan(opts.quadMax))
			opts.quadMax = model.sectionTimes[model.sectionTimes.size()-1];

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

			run(model, opts);
		}
	}

	::mpfr_free_cache();

	return 0;
}

