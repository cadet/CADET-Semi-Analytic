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

#include "MPReal.hpp"

#include <tclap/CmdLine.h>
#include "CliUtil.hpp"
#include "CliExtrapolationParser.hpp"

#include "VersionInfo.hpp"
#include "TclapCustomOutput.hpp"

#include "CenteredMoments.hpp"

#include "MPAD.hpp"

#ifdef CASEMA_USE_FADBAD
	#include <tadiff.h>
#endif
#include <cppad/cppad.hpp>

#include "ADMoments.hpp"

#include "Sequence.hpp"
#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"


template <typename real_t>
struct ProgramOptions
{
	std::size_t precision;
	std::size_t outPrecision;
	std::size_t order;
	std::string outFile;
	bool useCppAD;
	bool single;
	std::size_t nIter;
	bool verbose;

	real_t seqStart;
	real_t seqParam;

	std::size_t convTimes;
	real_t convThresholdAbs;
	real_t convThresholdRel;
	std::size_t consAgree;
};


mpfr::mpreal subfactorial(std::size_t n)
{
	return floor(mpfr::fac_ui(n) / exp(mpfr::mpreal(0)) + mpfr::mpreal(0.5));
}


template <typename real_t>
struct UniformDist
{
	UniformDist() : b(2) { }

	const std::string name() const { return "UniformDist(" + std::to_string(b) + ")"; }
	
	template <typename eval_t>
	const eval_t operator()(const eval_t& x) const
	{
		// Real function: Uniform distribution on [0, b]
		return (real_t(1) - exp(-b * x)) / (b * x);
	}

	const std::vector<real_t> nonCentralMoments(std::size_t order) const
	{
		std::vector<real_t> moms(order+1, 1);
		real_t bPower(b);
		for (std::size_t i = 1; i <= order; ++i)
		{
			moms[i] = bPower / real_t(i+1);
			bPower *= b;
		}
		return moms;
	}

	const std::vector<real_t> centralMoments(std::size_t order) const
	{
		real_t bPower(1);

		std::vector<real_t> moms(order+1, 0);
		moms[0] = 1;
		if (order >= 1)
			moms[1] = b / real_t(2);
	   
		const real_t one(1);
		const real_t two(2);
		for (std::size_t i = 2; i <= order; i += 2)
		{
			bPower *= sqr(b / two);
			moms[i] = bPower / (one + i);
		}
		return moms;
	}

	real_t b;
};


template <typename real_t>
struct ExponentialDist
{
	ExponentialDist() : lambda(2) { }

	const std::string name() const { return "ExponentialDist(" + std::to_string(lambda) + ")"; }
	
	template <typename eval_t>
	const eval_t operator()(const eval_t& x) const
	{
		// Real function lambda * exp(-lambda * x)
		return lambda / (lambda + x);
	}

	const std::vector<real_t> nonCentralMoments(std::size_t order) const
	{
		std::vector<real_t> moms(order+1, 1);
		for (std::size_t i = 1; i <= order; ++i)
		{
			moms[i] = moms[i-1] * real_t(i) / lambda;
		}
		return moms;
	}

	const std::vector<real_t> centralMoments(std::size_t order) const
	{
		real_t cur(1);
		real_t prev(0);
		real_t last(1);
		real_t lambdaPower(lambda);

		std::vector<real_t> moms(order+1, 1);
		if (order >= 1)
			moms[1] = real_t(1) / lambda;
	   
		for (std::size_t i = 2; i <= order; ++i)
		{
			// Use recurrence formula to update subfactorial
			cur = (i-1) * (prev + last);
			last = prev;
			prev = cur;

			lambdaPower *= lambda;
			moms[i] = cur / lambdaPower;
		}
		return moms;
	}

	real_t lambda;
};


template <typename real_t>
casema::Sequence<real_t>* createSequence(const std::string& name, const real_t& start, const real_t& param)
{
	std::string lowName = name;
	std::transform(lowName.begin(), lowName.end(), lowName.begin(), ::tolower);

	casema::Sequence<real_t>* seq = nullptr;

	if (lowName == "linear")
		seq = new casema::LinearSequence<real_t>(start, param);
	else if (lowName == "loglin")
		seq = new casema::LogLinearSequence<real_t>(start, param);
	else if (lowName == "geom")
		seq = new casema::GeometricSequence<real_t>(start, param);

	return seq;
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const real_t& moment, const real_t& radius, std::size_t order, const std::vector<real_t>& ref)
{
	const real_t error = abs(moment - ref[order]);
	os << "order,moment,radius,error,logerror\n";
	os << order << "," << moment << "," << radius << "," << error << "," << log10(error) << "\n";
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const real_t& moment, std::size_t order, const std::vector<real_t>& ref)
{
	writeMomentsTable(os, moment, real_t(0), order, ref);
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, const std::vector<real_t>& ref)
{
	os << "order,moment,radius,error,logerror\n";
	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		const real_t error = abs(moments[i] - ref[i]);
		os << i << "," << moments[i] << ",0," << error << "," << log10(error) << "\n";
	}
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, const std::vector<real_t>& radius, const std::vector<real_t>& ref)
{
	os << "order,moment,radius,error,logerror\n";

	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		const real_t error = abs(moments[i] - ref[i]);
		os << i << "," << moments[i] << "," << radius[i] << "," << error << "," << log10(error) << "\n";
	}
}


template <typename real_t>
void writeVerboseInsightHeader(std::ostream& os, const casema::ConsensusEstimator<real_t>& estimator)
{
	os << "## x,";
	for (typename casema::ConsensusEstimator<real_t>::const_iterator it = estimator.begin(); it != estimator.end(); ++it)
		os << (*it)->name() << ",";
	os << "consensus,radius,error,logerror" << std::endl;
}


template <typename real_t>
void writeVerboseInsight(std::ostream& os, const real_t& pos, const casema::ConsensusEstimator<real_t>& estimator, const std::vector<real_t>& ref, std::size_t order)
{
	os << "## " << pos << ",";
	const std::vector<real_t>& vals = estimator.estimates();
	for (typename std::vector<real_t>::const_iterator it = vals.begin(); it != vals.end(); ++it)
		os << (*it) << ",";

	const real_t error = abs(estimator.center() - ref[order]);
	os << estimator.center() << "," << estimator.radius() << "," << error << "," << log10(error) << std::endl;
}


template <typename real_t>
void writeVerboseInsightFooter(std::ostream& os, const casema::ConsensusEstimator<real_t>& estimator)
{
	os << "## Converged,";
	for (typename casema::ConsensusEstimator<real_t>::const_iterator it = estimator.begin(); it != estimator.end(); ++it)
		os << ((static_cast<casema::ConvergenceMonitoredExtrapolationMethod<real_t>*>(*it))->converged() ? "YES" : "NO") << ",";
	os << "," << std::endl;
}


template <typename real_t>
void calcSingleMoments(std::ostream& os, casema::MomentGenerator<real_t>& momGen, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, 
	const real_t& mass, const ProgramOptions<real_t>& opts, const std::vector<real_t>& refNonCentral)
{
	os << "############ Non-central\n";

	if (opts.verbose)
	{
		os << "\n# Order: " << opts.order << "\n";
		writeVerboseInsightHeader(os, estimator);
	}

	for (std::size_t i = 0; i < opts.nIter; ++i)
	{
		const real_t pos = seq.next();
		const real_t val = momGen.moment(pos, mass, opts.order);

		estimator.nextPoint(val, pos);

		if (opts.verbose)
			writeVerboseInsight(os, pos, estimator, refNonCentral, opts.order);
	}

	if (opts.verbose)
	{
		os << std::endl;
		writeVerboseInsightFooter(os, estimator);
	}

	writeMomentsTable(os, estimator.center(), estimator.radius(), opts.order, refNonCentral);
}


template <typename real_t>
void calcAllMoments(std::ostream& os, casema::MomentGenerator<real_t>& momGen, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, 
	const real_t& mass, const ProgramOptions<real_t>& opts, const std::vector<real_t>& refNonCentral, const std::vector<real_t>& refCentral)
{
	// Cache moments in matrix, each moment is contiguous in memory, i.e. nth point of mth moment is momMatrix[m * opts.nIter + n]
	real_t* const momMatrix = new real_t[(opts.order+1) * opts.nIter];
	real_t* const points = new real_t[opts.nIter];
	for (std::size_t i = 0; i < opts.nIter; ++i)
	{
		points[i] = seq.next();
		const std::vector<real_t> output = momGen.moments(points[i], mass, opts.order);
		for (std::size_t j = 0; j < output.size(); ++j)
			momMatrix[j * opts.nIter + i] = output[j];
	}

	// Operate on every moment separately
	real_t const* curSegment = momMatrix;
	
	std::vector<real_t> finalMomEst(opts.order+1);
	std::vector<real_t> finalRadius(opts.order+1);

	os << "############ Non-central\n";

	for (std::size_t m = 0; m <= opts.order; ++m)
	{
		estimator.reset();

		if (opts.verbose)
		{
			os << "\n# Order: " << m << "\n";
			writeVerboseInsightHeader(os, estimator);
		}

		for (std::size_t i = 0; i < opts.nIter; ++i)
		{
			estimator.nextPoint(curSegment[i], points[i]);

			if (opts.verbose)
				writeVerboseInsight(os, points[i], estimator, refNonCentral, m);
		}
		finalMomEst[m] = estimator.center();
		finalRadius[m] = estimator.radius();

		if (opts.verbose)
			writeVerboseInsightFooter(os, estimator);

		// Advance segment pointer
		curSegment += opts.nIter;
	}

	delete[] points;
	delete[] momMatrix;

	if (opts.verbose)
		os << std::endl;
	writeMomentsTable(os, finalMomEst, finalRadius, refNonCentral);

	os << "############ Central\n";
	const std::vector<real_t> centMom = casema::centralMomentsFromNonCentral(finalMomEst);
	writeMomentsTable(os, centMom, refCentral);
}


template <typename real_t, template <class T> class Functor>
void run(std::ostream& os, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, const ProgramOptions<real_t>& opts)
{
	casema::MomentGenerator<real_t>* momGen = nullptr;
	Functor<real_t> f;

	if (opts.useCppAD)
		momGen = new casema::CppADMomentGenerator<real_t, Functor<real_t>>(f);
	else
	{
#ifdef CASEMA_USE_FADBAD
		momGen = new casema::FadBadMomentGenerator<real_t, Functor<real_t>>(f);
#endif            
	}

	const std::vector<real_t> refNonCentral = f.nonCentralMoments(opts.order);

	// Get 0th order moment
	const real_t mass = refNonCentral[0];

	if (opts.single)
		calcSingleMoments(os, *momGen, estimator, seq, mass, opts, refNonCentral);
	else
	{
		const std::vector<real_t> refCentral = f.centralMoments(opts.order);
		calcAllMoments(os, *momGen, estimator, seq, mass, opts, refNonCentral, refCentral);
	}

	delete momGen;
}


int main(int argc, char** argv)
{
	typedef mpfr::mpreal real_t;

	ProgramOptions<real_t> opts;
	std::string extraMethods;
	std::string sequenceType;
	std::string model;

	try
	{
		TCLAP::CmdLine cmd("Uses an inverse Laplace transform, algorithmic differentiation, and extrapolation to compute moments", ' ', "");

		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model", true, "", "Model"))->storeIn(&model);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write output to file (default: disabled)", false, std::string(), "File"))->storeIn(&opts.outFile);
		cmd >> (new TCLAP::SwitchArg("v", "verbose", "Write intermediate results (default: disabled)"))->storeIn(&opts.verbose);

#ifdef CASEMA_USE_FADBAD
		cmd >> (new TCLAP::SwitchArg("", "cppAD", "Use CppAD (default: FADBAD++)"))->storeIn(&opts.useCppAD);
#else
		opts.useCppAD = true;
#endif

		cmd >> (new TCLAP::ValueArg<std::size_t>("m", "order", "Maximum moment order (default: 4)", false, 4, "Int"))->storeIn(&opts.order);
		cmd >> (new TCLAP::SwitchArg("S", "single", "Calculate moment of maximum order only"))->storeIn(&opts.single);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::ValueArg<std::string>("e", "extrapolation", "Use extrapolation methods (default: none)", false, std::string(), "Method1,Method2"))->storeIn(&extraMethods);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "convtimes", "Number of successful convergence tests (default: 3)", false, 3, "Int"))->storeIn(&opts.convTimes);
		cmd >> (new TCLAP::ValueArg<real_t>("", "convabs", "Absolute threshold used in convergence test (default: 100 * MachineEps)", false, real_t(-1), "Double"))->storeIn(&opts.convThresholdAbs);
		cmd >> (new TCLAP::ValueArg<real_t>("", "convrel", "Relative threshold used in convergence test (default: 0)", false, real_t(0), "Double"))->storeIn(&opts.convThresholdRel);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "agree", "Number of extrapolators used for accuracy estimate (default: 1)", false, 1, "Int"))->storeIn(&opts.consAgree);

		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "iterations", "Maximum number of iterations (default: 50)", false, 50, "Int"))->storeIn(&opts.nIter);

		cmd >> (new TCLAP::ValueArg<std::string>("s", "seq", "Sequence type (default: geom)", false, "geom", "Type"))->storeIn(&sequenceType);
		cmd >> (new TCLAP::ValueArg<real_t>("", "seqstart", "Sequence start point (default: 1e-15)", false, real_t("1.0e-15"), "Double"))->storeIn(&opts.seqStart);
		cmd >> (new TCLAP::ValueArg<real_t>("", "seqparam", "Sequence parameter (default: 0.1)", false, real_t("0.1"), "Double"))->storeIn(&opts.seqParam);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(opts.precision));

	// Create extrapolation method
	casema::ConsensusEstimator<real_t> estimator;
	casema::addExtrapolationMethods(extraMethods, estimator, opts);

	// Create sequence generator
	casema::Sequence<real_t>* const seq = createSequence(sequenceType, opts.seqStart, opts.seqParam);
	if (!seq)
	{
		std::cout << "ERROR: Unknown sequence type \"" << sequenceType << "\"" << std::endl;
		return 0;
	}

	// Set default convergence threshold
	if (opts.convThresholdAbs < 0)
		opts.convThresholdAbs = real_t(100) * std::numeric_limits<real_t>::epsilon();

	// Set default output precision
	if (opts.outPrecision == 0)
		opts.outPrecision = opts.precision;

	std::ostream* os = &std::cout;
	if (!opts.outFile.empty())
	{
		os = new std::ofstream(opts.outFile, std::ofstream::out | std::ofstream::trunc);
	}
	std::ios_base::fmtflags curFlags = os->flags();
	std::streamsize curPrec = os->precision();

	os->flags(std::ios::scientific);
	os->precision(opts.outPrecision);

	// Select test function
	casema::util::toLower(model);
	if (model == "exp2")
		run<real_t, ExponentialDist>(*os, estimator, *seq, opts);
	else if (model == "uni2")
		run<real_t, UniformDist>(*os, estimator, *seq, opts);
	else
		std::cout << "ERROR: Unkown model " << model << std::endl;

	os->precision(curPrec);
	os->flags(curFlags);
	if (!opts.outFile.empty())
	{
		delete os;
		os = &std::cout;
	}

	if (seq)
		delete seq;

	::mpfr_free_cache();

	return 0;
}

