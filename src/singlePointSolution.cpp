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

#include <omp.h>

#include "MPReal.hpp"

#include <tclap/CmdLine.h>
#include "CliUtil.hpp"

#include "TclapCustomOutput.hpp"
#include "VersionInfo.hpp"
#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "MPComplex.hpp"
#include "ExtrapolatedDurbinsMethod.hpp"
#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"

#include "CliErrorParser.hpp"

typedef casema::laplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


template <typename real_t>
struct ProgramOptions
{
	std::size_t summands;
	std::size_t minSummands;
	std::size_t precision;
	std::size_t outPrecision;
	real_t timePoint;
	real_t sigma;
	real_t errorTrunc;
	real_t errorCons;
	real_t error;
	real_t errorWeight;
	real_t refSol;
	std::string inFile;
	std::string outFile;
	std::size_t numThreads;

	bool noMeta;
	bool single;
	bool finalOnly;
	std::string extraMethod;
	std::size_t convTimes;
	real_t convThresholdAbs;
	real_t convThresholdRel;
};


template <typename real_t>
class OutputWriter
{
public:

	OutputWriter(std::ostream& outStream, const ProgramOptions<real_t>& opts) : _outStream(outStream), _opts(opts) { }

	void loopHook(const std::size_t summand, const real_t& seriesVal, const real_t& curEst, const real_t& f0, const real_t& factor)
	{
		const real_t endSeriesVal = (f0 + seriesVal) * factor;
		const real_t endCurEst = (f0 + curEst) * factor;

		if (_opts.finalOnly)
			_outStream << summand << "," << endSeriesVal << "," << endCurEst;
		else
			_outStream << summand << "," << seriesVal << "," << curEst << "," << endSeriesVal << "," << endCurEst;

		if (!isnan(_opts.refSol))
		{
			const real_t error = abs(endCurEst - _opts.refSol);
			_outStream << "," << error << "," << log10(error);
		}

		_outStream << "\n";
	}

	void loopHook(const std::size_t summand, const real_t& cosVal, const real_t& cosEst, const real_t& sinVal, const real_t& sinEst, const real_t& f0, const real_t& factor, bool convCos, bool convSin)
	{
		const real_t endSeriesVal = (f0 + cosVal - sinVal) * factor;
		const real_t endCurEst = (f0 + cosEst - sinEst) * factor;

		char conv = '-';
		if (convCos && !convSin)
			conv = 'C';
		else if (convSin && !convCos)
			conv = 'S';
		else if (convSin && convCos)
			conv = 'B';

		if (_opts.finalOnly)
			_outStream << summand << "," << conv << "," << endSeriesVal << "," << endCurEst;
		else
			_outStream << summand << "," << conv << "," << cosVal << "," << cosEst << "," << sinVal << "," << sinEst << "," << endSeriesVal << "," << endCurEst;

		if (!isnan(_opts.refSol))
		{
			const real_t error = abs(endCurEst - _opts.refSol);
			_outStream << "," << error << "," << log10(error);
		}

		_outStream << "\n";
	}

	void loopStartHook(bool combined) const
	{
		if (combined)
		{
			if (_opts.finalOnly)
				_outStream << "summand,finalSeries,finalLimit";
			else
				_outStream << "summand,series,limit,finalSeries,finalLimit";
		}
		else
		{
			if (_opts.finalOnly)
				_outStream << "summand,conv,finalSeries,finalLimit";
			else
				_outStream << "summand,conv,cosSeries,cosLimit,sinSeries,sinLimit,finalSeries,finalLimit";
		}

		if (!isnan(_opts.refSol))
			_outStream << ",error,logerror";

		_outStream << "\n";
	}

protected:
	std::ostream& _outStream;
	const ProgramOptions<real_t>& _opts;
};


template <typename real_t>
class RecordingHookPolicy
{
public:
	RecordingHookPolicy() : _writer(nullptr) { }

	const OutputWriter<real_t>* const writer() const { return _writer; }
	void writer(OutputWriter<real_t>* const writer) { _writer = writer; }

protected:

	void loopHook(const std::size_t summand, const real_t& seriesVal, const real_t& curEst, const real_t& f0, const real_t& factor) const
	{
		_writer->loopHook(summand, seriesVal, curEst, f0, factor);
	}

	void loopHook(const std::size_t summand, const real_t& cosVal, const real_t& cosEst, const real_t& sinVal, const real_t& sinEst, const real_t& f0, const real_t& factor, bool convCos, bool convSin) const
	{
		_writer->loopHook(summand, cosVal, cosEst, sinVal, sinEst, f0, factor, convCos, convSin);
	}

	void loopStartHook(bool combined) const
	{
		_writer->loopStartHook(combined);
	}

	mutable OutputWriter<real_t>* _writer;
};


template <typename real_t, typename complex_t, template <class T> class Extrapolator_t> using InversionMethod = casema::ExtrapolatedDurbinsMethod<real_t, complex_t, casema::SingleExtrapolatorPolicy<real_t, Extrapolator_t>, RecordingHookPolicy>;

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
	os << "# Convergence: " << opts.convTimes << " times below " << opts.convThresholdAbs << " (abs) or " << opts.convThresholdRel << " (rel) \n";
	os << "# Summands: Min " << opts.minSummands << " max " << opts.summands << "\n";
	os << "# Abscissa (a): " << opts.sigma << "\n";
	os << "# Truncation error: " << opts.errorTrunc << "\n";
	os << "# Consistency error: " << opts.errorCons << "\n";
	os << "# Error: " << opts.error << "\n";
	os << "# Error weight: " << opts.errorWeight << "\n";
	os << "# Time point: " << opts.timePoint << "\n";
	os << "# Sin() and cos() separately: " << (opts.single ? "Yes" : "No") << "\n";

	if (!isnan(opts.refSol))
		os << "# Reference solution: " << opts.refSol << "\n";

	os << "# Extrapolation: " << opts.extraMethod << std::endl;
}


template <class Solution_t, template <class T> class Extrapolator_t>
void performInversion(Solution_t& solution, const casema::ModelData<mpfr::mpreal>& model, const ProgramOptions<mpfr::mpreal>& opts, mpfr::mpreal& radius, std::size_t& numIter, OutputWriter<mpfr::mpreal>* const writer)
{
	InversionMethod<mpfr::mpreal, mpfr::mpcomplex, Extrapolator_t> durbin(opts.precision, opts.summands, opts.minSummands, casema::maxSimulationTime(model), opts.sigma);
	durbin.writer(writer);
	durbin.options(opts);

	if (opts.single)
		durbin.invertSingle(model.outletTimes, 0, solution, &radius, &numIter);
	else
		durbin.invertCombined(model.outletTimes, 0, solution, &radius, &numIter);
}


template <class Solution_t>
void invert(Solution_t& solution, const casema::ModelData<mpfr::mpreal>& model, const ProgramOptions<mpfr::mpreal>& opts, mpfr::mpreal& radius, std::size_t& numIter, OutputWriter<mpfr::mpreal>* const writer)
{
	std::cout << solution.name() << std::endl;
	std::string lowName = opts.extraMethod;
	std::transform(lowName.begin(), lowName.end(), lowName.begin(), ::tolower);
	if (lowName == "ads")
		performInversion<Solution_t, casema::AitkenDeltaSquaredMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "wem")
		performInversion<Solution_t, casema::WynnEpsilonMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "wrm")
		performInversion<Solution_t, casema::WynnRhoMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "iad")
		performInversion<Solution_t, casema::IteratedAitkenDeltaSquaredMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "lum")
		performInversion<Solution_t, casema::LevinUMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "ltm")
		performInversion<Solution_t, casema::LevinTMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "ibt")
		performInversion<Solution_t, casema::IteratedBrezinskiThetaMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "btm")
		performInversion<Solution_t, casema::BrezinskiThetaMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "nam")
		performInversion<Solution_t, casema::NevilleAitkenMethod>(solution, model, opts, radius, numIter, writer);
	else if (lowName == "rem")
		performInversion<Solution_t, casema::RichardsonMethod>(solution, model, opts, radius, numIter, writer);
	else
		performInversion<Solution_t, casema::IdentityExtrapolation>(solution, model, opts, radius, numIter, writer);
}


void run(casema::ModelData<mpfr::mpreal>& model, const ProgramOptions<mpfr::mpreal>& opts)
{
	model.outletTimes.clear();
	model.outletTimes.push_back(opts.timePoint);

	std::cout << "Starting Laplace inversion for Model" << std::endl;

	mpfr::mpreal radius = 0;
	std::size_t numIter = 0;

	std::ostream* outStream = &std::cout;
	std::ios_base::fmtflags curFlags = outStream->flags();
	std::streamsize curPrec = outStream->precision();

	outStream->flags(std::ios::scientific);
	outStream->precision(opts.outPrecision);

	writeMeta(std::cout, model, opts);
	if (!opts.outFile.empty())
	{
		outStream = new std::ofstream(opts.outFile, std::ofstream::out | std::ofstream::trunc);
		if (!opts.noMeta)
			writeMeta(*outStream, model, opts);

		outStream->flags(std::ios::scientific);
		outStream->precision(opts.outPrecision);
	}

	OutputWriter<mpfr::mpreal> writer(*outStream, opts);

	Inlet_t inlet(model);
	if (model.kineticBinding)
	{
		typedef casema::laplaceSolution::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		invert<Solution_t>(solution, model, opts, radius, numIter, &writer);
	}
	else
	{
		typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
		Solution_t solution(model, inlet);
		invert<Solution_t>(solution, model, opts, radius, numIter, &writer);
	}

	if (!opts.noMeta)
	{
		(*outStream) << "# Radius: " << radius << std::endl;
		(*outStream) << "# Iterations: " << numIter << std::endl;
	}

	if (!opts.outFile.empty())
	{
		delete outStream;
		outStream = &std::cout;
	}

	outStream->precision(curPrec);
	outStream->flags(curFlags);
}


int main(int argc, char** argv)
{
	ProgramOptions<mpfr::mpreal> opts;

	try
	{
		TCLAP::CustomOutput customOut("singlePointSolution");
		TCLAP::CmdLine cmd("Uses a numerical inverse Laplace transform to solve GRM models for a single time point using a single extrapolation method", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<std::size_t>("t", "threads", "Number of threads (default: 4)", false, 4, "Int"))->storeIn(&opts.numThreads);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::SwitchArg("S", "single", "Extrapolate sin() and cos() terms separately"))->storeIn(&opts.single);
		cmd >> (new TCLAP::SwitchArg("f", "final", "Write extrapolation with factor only (default: both)"))->storeIn(&opts.finalOnly);
		cmd >> (new TCLAP::SwitchArg("", "nometa", "Do not write meta info to file"))->storeIn(&opts.noMeta);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("r", "ref", "Reference solution for error calculation (default: none)", false, mpfr::mpreal().setNan(), "Float"))->storeIn(&opts.refSol);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("T", "time", "Time point (default: 10.0)", false, 10, "Float"))->storeIn(&opts.timePoint);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method, used as safety margin if error (-e) is given", false, 0, "Float"))->storeIn(&opts.sigma);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Error threshold (default: 1e-10)", false, 1e-10, "Float"))->storeIn(&opts.error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Weight used to distribute error onto consistency and truncation (default: 0.5)", false, 0.5, "Float"))->storeIn(&opts.errorWeight);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "sum", "Maximum number of summands in Durbin's method", false, 0, "Int"))->storeIn(&opts.summands);
		cmd >> (new TCLAP::ValueArg<std::size_t>("m", "minsum", "Minimum number of summands in Durbin's method (default: 0)", false, 0, "Int"))->storeIn(&opts.minSummands);

		cmd >> (new TCLAP::ValueArg<std::string>("x", "extrapolation", "Use extrapolation methods (default: none)", false, std::string(), "Method"))->storeIn(&opts.extraMethod);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "convtimes", "Number of successful convergence tests (default: 3)", false, 3, "Int"))->storeIn(&opts.convTimes);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "convabs", "Absolute threshold used in convergence test (default: 100 * MachineEps)", false, mpfr::mpreal(-1), "Double"))->storeIn(&opts.convThresholdAbs);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("", "convrel", "Relative threshold used in convergence test (default: 0)", false, mpfr::mpreal(0), "Double"))->storeIn(&opts.convThresholdRel);

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

			run(model, opts);
		}
	}

	::mpfr_free_cache();

	return 0;
}

