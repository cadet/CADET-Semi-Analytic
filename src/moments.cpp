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
#include "CliExtrapolationParser.hpp"

#include "VersionInfo.hpp"
#include "TclapCustomOutput.hpp"

#include "ModelReaderHelper.hpp"
#include "ModelDataChecker.hpp"

#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "AnalyticalMoments.hpp"
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

#include "CSVReader.hpp"


template <typename real_t>
struct ProgramOptions
{
	std::size_t precision;
	std::size_t outPrecision;
	std::size_t order;
	std::string inFile;
	std::string outFile;
	bool useCppAD;
	bool single;
	std::size_t nIter;
	bool verbose;
	bool analytic;

	real_t seqStart;
	real_t seqParam;

	std::size_t convTimes;
	real_t convThresholdAbs;
	real_t convThresholdRel;
	std::size_t consAgree;
	bool stopSingleConv;

	std::vector<mpfr::mpreal>* refCentral;
	std::vector<mpfr::mpreal>* refNonCentral;
};


std::vector<mpfr::mpreal>* loadReference(std::string& input)
{
	casema::util::trim(input);
	if (input.empty())
		return nullptr;

	const std::vector<std::string> tokens = casema::util::split(input, ",");
	if (tokens.empty())
		return nullptr;

	casema::CSVReader reader;
	if (!reader.readFile(tokens[0]))
		return nullptr;

	if (tokens.size() == 2)
	{
		if (reader.contains(tokens[1]))
		{
			const std::size_t idx = reader.indexOf(tokens[1]);
			return new std::vector<mpfr::mpreal>(reader.column(idx).begin(), reader.column(idx).end());
		}
		
		std::cout << "WARNING: Could not find column " << tokens[1] << " in file " << tokens[0] << ", ignoring reference" << std::endl;
		return nullptr;
	}

	if (reader.contains("moment"))
	{
		const std::size_t idx = reader.indexOf("moment");
		return new std::vector<mpfr::mpreal>(reader.column(idx).begin(), reader.column(idx).end());
	}

	std::cout << "WARNING: Could not identify moments column in file " << tokens[0] << ", ignoring reference" << std::endl;
	return nullptr;
}


template <typename real_t>
void writeMeta(std::ostream& os, const casema::ModelData<real_t>& model, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, const ProgramOptions<real_t>& opts)
{
	const real_t mass = casema::detail::injectedMass(model);

	os << "# Version: " << casema::getVersion() << "\n";
	os << "# Commit: " << casema::getCommitHash() << "\n";
	os << "# Model: " << opts.inFile << "\n";
	os << "# Binding: " << (model.kineticBinding ? "Dynamic" : "Rapid-Equlibrium") << "\n";
	os << "# Injected mass: " << mass << "\n";
	os << "# Precision " << opts.precision << " digits = " << mpfr::digits2bits(opts.precision) << " bits\n";
	os << "# Extrapolation: " << estimator.name() << "\n";
	os << "# Convergence: " << opts.convTimes << " times below " << opts.convThresholdAbs << " (abs) or " << opts.convThresholdRel << " (rel) \n";
	os << "# Consensus: " << opts.consAgree << " extrapolators\n";
	os << "# MaxIter: " << opts.nIter << "\n";
	os << "# Sequence: " << seq.name() << "\n";
	os << "# SeqStart: " << opts.seqStart << " SeqParam " << opts.seqParam << "\n";
	os << "# Moments: " << (opts.analytic ? "Analytical" : "AD-Exp") << std::endl;
	if (!opts.analytic)
		os << "# AD: " << (opts.useCppAD ? "CppAD" : "FADBAD++") << std::endl;
}


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
void writeMomentsTable(std::ostream& os, const real_t& moment, const real_t& radius, std::size_t order, std::vector<real_t> const* const ref)
{
	if (ref)
		os << "order,moment,radius,error,logerror,relerror,logrelerror\n";
	else
		os << "order,moment,radius\n";
	
	os << order << "," << moment << "," << radius << "\n";
	if (ref)
	{
		if (ref->size() > order)
		{
			const real_t error = abs(moment - (*ref)[order]);
			const real_t relError = error / (*ref)[order];
			os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
		}
		else
			os << ",,,,";
	}

	os << "\n";
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const real_t& moment, std::size_t order, std::vector<real_t> const* const ref)
{
	writeMomentsTable(os, moment, real_t(0), order, ref);
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, std::vector<real_t> const* const ref)
{
	if (ref)
		os << "order,moment,radius,error,logerror,relerror,logrelerror\n";
	else
		os << "order,moment,radius\n";
	
	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		os << i << "," << moments[i] << ",0";
		if (ref)
		{
			if (ref->size() > i)
			{
				const real_t error = abs(moments[i] - (*ref)[i]);
				const real_t relError = error / (*ref)[i];
				os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
			}
			else
				os << ",,,,";
		}
		os << "\n";
	}
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, const std::vector<real_t>& radius, std::vector<real_t> const* const ref)
{
	if (ref)
		os << "order,moment,radius,error,logerror,relerror,logrelerror\n";
	else
		os << "order,moment,radius\n";

	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		os << i << "," << moments[i] << "," << radius[i];
		if (ref)
		{
			if (ref->size() > i)
			{
				const real_t error = abs(moments[i] - (*ref)[i]);
				const real_t relError = error / (*ref)[i];
				os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
			}
			else
				os << ",,,,";
		}
		os << "\n";
	}
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, const std::vector<std::size_t>& iters, std::vector<real_t> const* const ref)
{
	if (ref)
		os << "order,moment,radius,iter,error,logerror,relerror,logrelerror\n";
	else
		os << "order,moment,radius,iter\n";
	
	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		os << i << "," << moments[i] << ",0," << iters[i];
		if (ref)
		{
			if (ref->size() > i)
			{
				const real_t error = abs(moments[i] - (*ref)[i]);
				const real_t relError = error / (*ref)[i];
				os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
			}
			else
				os << ",,,,";
		}
		os << "\n";
	}
}


template <typename real_t>
void writeMomentsTable(std::ostream& os, const std::vector<real_t>& moments, const std::vector<real_t>& radius, const std::vector<std::size_t>& iters, std::vector<real_t> const* const ref)
{
	if (ref)
		os << "order,moment,radius,iter,error,logerror,relerror,logrelerror\n";
	else
		os << "order,moment,radius,iter\n";

	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		os << i << "," << moments[i] << "," << radius[i] << "," << iters[i];
		if (ref)
		{
			if (ref->size() > i)
			{
				const real_t error = abs(moments[i] - (*ref)[i]);
				const real_t relError = error / (*ref)[i];
				os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
			}
			else
				os << ",,,,";
		}
		os << "\n";
	}
}


template <typename real_t>
void calcAnalyticalMoments(std::ostream& os, const casema::ModelData<real_t>& model, const ProgramOptions<real_t>& opts)
{
	os << "############ Non-central\n";
	const std::vector<real_t> moments = casema::analyticalMoments(model);
	if (opts.single)
		writeMomentsTable(os, moments[opts.order], opts.order, opts.refNonCentral);
	else
		writeMomentsTable(os, moments, opts.refNonCentral);

	os << "############ Central\n";
	const std::vector<mpfr::mpreal> momCent = casema::analyticalCentralMoments(model);
	if (opts.single)
		writeMomentsTable(os, momCent[opts.order], opts.order, opts.refCentral);
	else
		writeMomentsTable(os, momCent, opts.refCentral);
}


template <typename real_t>
void writeVerboseInsightHeader(std::ostream& os, const casema::ConsensusEstimator<real_t>& estimator, bool ref)
{
	os << "## x,";
	for (typename casema::ConsensusEstimator<real_t>::const_iterator it = estimator.begin(); it != estimator.end(); ++it)
		os << (*it)->name() << ",";
	os << "consensus,radius"; 

	if (ref)
		os << ",error,logerror,relerror,logrelerror";

	os << std::endl;
}


template <typename real_t>
void writeVerboseInsight(std::ostream& os, const real_t& pos, const casema::ConsensusEstimator<real_t>& estimator, std::vector<real_t> const* const ref, std::size_t order)
{
	os << "## " << pos << ",";
	const std::vector<real_t>& vals = estimator.estimates();
	for (typename std::vector<real_t>::const_iterator it = vals.begin(); it != vals.end(); ++it)
		os << (*it) << ",";
	os << estimator.center() << "," << estimator.radius();

	if (ref)
	{
		if (ref->size() > order)
		{
			const real_t error = abs(estimator.center() - (*ref)[order]);
			const real_t relError = error / (*ref)[order];
			os << "," << error << "," << log10(error) << "," << relError << "," << log10(relError);
		}
		else
			os << ",,,,";
	}

	os << std::endl;
}


template <typename real_t>
void writeVerboseInsightFooter(std::ostream& os, const casema::ConsensusEstimator<real_t>& estimator)
{
	os << "## Converged,";
	for (typename casema::ConsensusEstimator<real_t>::const_iterator it = estimator.begin(); it != estimator.end(); ++it)
		os << ((static_cast<casema::ConvergenceMonitoredExtrapolationMethod<real_t>*>(*it))->converged() ? "YES" : "NO") << ",";
	os << std::endl;
}


template <typename real_t>
void calcSingleMoments(std::ostream& os, casema::MomentGenerator<real_t>& momGen, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, const real_t& mass, const ProgramOptions<real_t>& opts)
{
	os << "############ Non-central\n";

	if (opts.verbose)
	{
		os << "\n# Order: " << opts.order << "\n";
		writeVerboseInsightHeader(os, estimator, opts.refNonCentral);
	}

	for (std::size_t i = 0; i < opts.nIter; ++i)
	{
		const real_t pos = seq.next();
		const real_t val = momGen.moment(pos, mass, opts.order);

		estimator.nextPoint(val, pos);

		if (opts.verbose)
			writeVerboseInsight(os, pos, estimator, opts.refNonCentral, opts.order);

		if (opts.stopSingleConv && (estimator.numEstimators() == 1))
		{
			casema::ConvergenceMonitoredExtrapolationMethod<real_t>* convMon = static_cast<casema::ConvergenceMonitoredExtrapolationMethod<real_t>*>(*estimator.begin());
			if (convMon->converged())
			{
				break;
			}
		}
	}

	if (opts.verbose)
	{
		os << std::endl;
		writeVerboseInsightFooter(os, estimator);
	}

	writeMomentsTable(os, estimator.center(), estimator.radius(), opts.order, opts.refNonCentral);
}


template <typename real_t>
void calcAllMoments(std::ostream& os, casema::MomentGenerator<real_t>& momGen, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, const real_t& mass, const ProgramOptions<real_t>& opts)
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
	std::vector<std::size_t> iterations(opts.order+1);

	os << "############ Non-central\n";

	for (std::size_t m = 0; m <= opts.order; ++m)
	{
		estimator.reset();

		if (opts.verbose)
		{
			os << "\n# Order: " << m << "\n";
			writeVerboseInsightHeader(os, estimator, opts.refNonCentral);
		}

		for (std::size_t i = 0; i < opts.nIter; ++i)
		{
			estimator.nextPoint(curSegment[i], points[i]);

			if (opts.verbose)
				writeVerboseInsight(os, points[i], estimator, opts.refNonCentral, m);

			if (opts.stopSingleConv && (estimator.numEstimators() == 1))
			{
				casema::ConvergenceMonitoredExtrapolationMethod<real_t>* convMon = static_cast<casema::ConvergenceMonitoredExtrapolationMethod<real_t>*>(*estimator.begin());
				if (convMon->converged())
				{
					iterations[m] = i+1;
					break;
				}
			}
		}
		finalMomEst[m] = estimator.center();
		finalRadius[m] = estimator.radius();

		if (iterations[m] == 0)
			iterations[m] = opts.nIter;

		if (opts.verbose)
			writeVerboseInsightFooter(os, estimator);

		// Advance segment pointer
		curSegment += opts.nIter;
	}

	delete[] points;
	delete[] momMatrix;

	if (opts.verbose)
		os << std::endl;
	writeMomentsTable(os, finalMomEst, finalRadius, iterations, opts.refNonCentral);

	os << "############ Central\n";
	const std::vector<real_t> centMom = casema::centralMomentsFromNonCentral(finalMomEst);
	writeMomentsTable(os, centMom, opts.refCentral);
}


template <typename real_t>
void run(std::ostream& os, const casema::ModelData<real_t>& model, casema::ConsensusEstimator<real_t>& estimator, casema::Sequence<real_t>& seq, const ProgramOptions<real_t>& opts)
{
	typedef casema::laplaceSolution::Inlet<real_t, real_t> Inlet_t;
	typedef casema::laplaceSolution::SingleComponentLinearDynamic<real_t, real_t, Inlet_t> SolutionDyn_t;
	typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<real_t, real_t, Inlet_t> SolutionREq_t;

	casema::MomentGenerator<real_t>* momGen = nullptr;
	Inlet_t inlet(model);

	if (opts.verbose)
		writeMeta(os, model, estimator, seq, opts);

	if (opts.analytic)
	{
		if (!model.kineticBinding && (opts.order <= 4))
			calcAnalyticalMoments(os, model, opts);
		else
			std::cout << "ERROR: Analytical moments up to 4th order are only supported for rapid-equlibrium binding (see Qamar et al. 2014)" << std::endl;

		return;
	}

	if (model.kineticBinding)
	{
		if (opts.useCppAD)
			momGen = new casema::CppADMomentGenerator<real_t, SolutionDyn_t>(SolutionDyn_t(model, inlet));
		else
		{
#ifdef CASEMA_USE_FADBAD
			momGen = new casema::FadBadMomentGenerator<real_t, SolutionDyn_t>(SolutionDyn_t(model, inlet));
#endif            
		}
	}
	else
	{
		if (opts.useCppAD)
			momGen = new casema::CppADMomentGenerator<real_t, SolutionREq_t>(SolutionREq_t(model,inlet));
		else
		{
#ifdef CASEMA_USE_FADBAD
			momGen = new casema::FadBadMomentGenerator<real_t, SolutionREq_t>(SolutionREq_t(model, inlet));
#endif            
		}
	}

	// Get 0th order moment
	const real_t mass = casema::detail::injectedMass(model);

	if (opts.single)
		calcSingleMoments(os, *momGen, estimator, seq, mass, opts);
	else
		calcAllMoments(os, *momGen, estimator, seq, mass, opts);

	delete momGen;
}


int main(int argc, char** argv)
{
	typedef mpfr::mpreal real_t;

	ProgramOptions<real_t> opts;
	std::string extraMethods;
	std::string sequenceType;
	std::string refCentral;
	std::string refNonCentral;

	try
	{
		TCLAP::CustomOutput customOut("moments");
		TCLAP::CmdLine cmd("Uses Laplace transform, algorithmic differentiation, and extrapolation to compute moments of GRM models", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file", true, "", "File"))->storeIn(&opts.inFile);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write output to file (default: disabled)", false, std::string(), "File"))->storeIn(&opts.outFile);
		cmd >> (new TCLAP::SwitchArg("v", "verbose", "Write intermediate results (default: disabled)"))->storeIn(&opts.verbose);

#ifdef CASEMA_USE_FADBAD
		cmd >> (new TCLAP::SwitchArg("", "cppAD", "Use CppAD (default: FADBAD++)"))->storeIn(&opts.useCppAD);
#else
		opts.useCppAD = true;
#endif

		cmd >> (new TCLAP::ValueArg<std::size_t>("m", "order", "Maximum moment order (default: 4)", false, 4, "Int"))->storeIn(&opts.order);
		cmd >> (new TCLAP::SwitchArg("S", "single", "Calculate moment of maximum order only"))->storeIn(&opts.single);
		cmd >> (new TCLAP::SwitchArg("a", "analytic", "Analytical moments up to 4th order by Qamar et al. (2014) for rapid-equlibrium"))->storeIn(&opts.analytic);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::ValueArg<std::string>("e", "extrapolation", "Use extrapolation methods (default: none)", false, std::string(), "Method1,Method2"))->storeIn(&extraMethods);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "convtimes", "Number of successful convergence tests (default: 3)", false, 3, "Int"))->storeIn(&opts.convTimes);
		cmd >> (new TCLAP::ValueArg<real_t>("", "convabs", "Absolute threshold used in convergence test (default: 100 * MachineEps)", false, real_t(-1), "Double"))->storeIn(&opts.convThresholdAbs);
		cmd >> (new TCLAP::ValueArg<real_t>("", "convrel", "Relative threshold used in convergence test (default: 0)", false, mpfr::mpreal(0), "Double"))->storeIn(&opts.convThresholdRel);
		cmd >> (new TCLAP::ValueArg<std::size_t>("", "agree", "Number of extrapolators used for accuracy estimate (default: 1)", false, 1, "Int"))->storeIn(&opts.consAgree);
		cmd >> (new TCLAP::SwitchArg("", "stopsingleconv", "Stop extrapolation on convergence if single extrapolator is used"))->storeIn(&opts.stopSingleConv);

		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "iterations", "Maximum number of iterations (default: 50)", false, 50, "Int"))->storeIn(&opts.nIter);

		cmd >> (new TCLAP::ValueArg<std::string>("s", "seq", "Sequence type (default: geom)", false, "geom", "Type"))->storeIn(&sequenceType);
		cmd >> (new TCLAP::ValueArg<real_t>("", "seqstart", "Sequence start point (default: 1e-15)", false, real_t("1.0e-15"), "Double"))->storeIn(&opts.seqStart);
		cmd >> (new TCLAP::ValueArg<real_t>("", "seqparam", "Sequence parameter (default: 0.1)", false, real_t("0.1"), "Double"))->storeIn(&opts.seqParam);

		cmd >> (new TCLAP::ValueArg<std::string>("", "refcentral", "Reference central moments for error calculation (default: none; format: CSVfile,columnName)", false, "", "File,Column"))->storeIn(&refCentral);
		cmd >> (new TCLAP::ValueArg<std::string>("", "refnoncentral", "Reference non-central moments for error calculation (default: none; format: CSVfile,columnName)", false, "", "File,Column"))->storeIn(&refNonCentral);

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

	// Load model
	casema::ModelData<mpfr::mpreal> model;
	const bool loaded = casema::readModel<mpfr::mpreal>(opts.inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		if (casema::checkModelForCompatibility(model))        
		{
			if (!opts.verbose || (opts.verbose && !opts.outFile.empty()))
				writeMeta(std::cout, model, estimator, *seq, opts);

			// Load reference moments
			opts.refCentral = nullptr;
			opts.refNonCentral = nullptr;
			if (!refCentral.empty())
				opts.refCentral = loadReference(refCentral);
			if (!refNonCentral.empty())
				opts.refNonCentral = loadReference(refNonCentral);

			if (opts.outFile.empty())
			{
				std::ios_base::fmtflags curFlags = std::cout.flags();
				std::streamsize curPrec = std::cout.precision();

				std::cout.flags(std::ios::scientific);
				std::cout.precision(opts.outPrecision);

				run(std::cout, model, estimator, *seq, opts);

				std::cout.precision(curPrec);
				std::cout.flags(curFlags);
			}
			else
			{
				std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);

				std::ios_base::fmtflags curFlags = fs.flags();
				std::streamsize curPrec = fs.precision();

				fs.flags(std::ios::scientific);
				fs.precision(opts.outPrecision);

				run(fs, model, estimator, *seq, opts);

				fs.precision(curPrec);
				fs.flags(curFlags);
			}

			if (opts.refCentral)
				delete opts.refCentral;
			if (opts.refNonCentral)
				delete opts.refNonCentral;
		}
	}

	if (seq)
		delete seq;

	::mpfr_free_cache();

	return 0;
}

