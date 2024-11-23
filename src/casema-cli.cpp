// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-present: Samuel Leweke¹² and the CADET-Semi-Analytic
//  authors, see the AUTHORS file
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "CompileTimeConfig.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <algorithm>
#include <fstream>
#include <memory>

#include <Eigen/Core>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "MPReal.hpp"

#include <tclap/CmdLine.h>
#include "TclapUtils.hpp"
#include "VersionInfo.hpp"

#include "MPComplex.hpp"
#include "DurbinsMethod.hpp"
#include "LaplaceError.hpp"
#include "BesselZeros.hpp"
#include "SlicedVector.hpp"
#include "TimeReader.hpp"
#include "ProgressBar.hpp"

#include "io/ParameterProviderImpl.hpp"
#include "ModelBuilder.hpp"
#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"


struct ProgramOptions
{
	std::size_t laplaceSummands;
	std::size_t hankelSummands;
	std::size_t precision;
	std::size_t outPrecision;
	mpfr::mpreal abscissa;
	mpfr::mpreal errorTrunc;
	mpfr::mpreal errorCons;
	mpfr::mpreal error;
	mpfr::mpreal errorWeight;
	std::string inFile;
	std::string outFile;
	int numThreads;
	bool kahan;
	bool ignoreCSTR;
};

class ProgressBarUpdater
{
public:
	ProgressBarUpdater() : _prevPhase(0) { _pb.begin(); }
	~ProgressBarUpdater() { _pb.finish(); }

	void update(double progress, int phase)
	{
		if ((phase == 1) && (_prevPhase == 0))
		{
			_prevPhase = 1;
			_pb.finish();
			_pb.begin();
		}

		if (phase == 0)
			_pb.update(progress, "Eval Laplace");
		else
			_pb.update(progress, "Eval Solution");
	}

private:
	casema::ProgressBar _pb;
	int _prevPhase;
};

#ifdef ENABLE_HDF5
void writeMetaAndResultToH5(const ProgramOptions& opts, const casema::model::ModelSystem* sys, const std::vector<mpfr::mpreal>& solutionTimes, const casema::util::SlicedVector<mpfr::mpreal>& output, const double time_sim, const int timeOffset)
{
	casema::io::IFileWriter* writer = casema::io::createWriter("h5");

	bool split_ports = true; // from CADET file format
	bool single_as_multi_port = false; // from CADET file format

	if (opts.outFile.empty()) // write to input file
	{
		casema::io::IFileReader* reader = casema::io::createReader("h5");
		reader->openFile(opts.inFile, "r");
		reader->pushGroup("input");
		reader->pushGroup("return");
		if (reader->exists("SINGLE_AS_MULTI_PORT"))
			single_as_multi_port = reader->getBool("SINGLE_AS_MULTI_PORT");

		reader->closeFile();
		writer->openFile(opts.inFile, "rw");
	}
	else if (opts.outFile == opts.inFile)
	{
		casema::io::IFileReader* reader = casema::io::createReader("h5");
		reader->openFile(opts.outFile, "r");
		reader->pushGroup("input");
		reader->pushGroup("return");
		if (reader->exists("SINGLE_AS_MULTI_PORT"))
			single_as_multi_port = reader->getBool("SINGLE_AS_MULTI_PORT");

		reader->closeFile();
		writer->openFile(opts.outFile, "rw");
	}
	else // create output file
		writer->openFile(opts.outFile, "co");

	// Write meta data
	writer->deleteGroup("meta");
	writer->pushGroup("meta");
	writer->writeString("CASEMA_VERSION", casema::getVersion());
	writer->writeString("CASEMA_BRANCH", casema::getBranchRefspec());
	writer->writeString("CASEMA_COMMIT", casema::getCommitHash());
	writer->writeInt("FILE_FORMAT", casema::getFileFormat());
	writer->writeDouble("TIME_SIM", time_sim);
	{
		std::ostringstream oss;
		oss << opts.errorTrunc + opts.errorCons;
		writer->writeString("ERROR", oss.str());
	}
	writer->popGroup();

	// Write program options
	writer->deleteGroup("program_options");
	writer->pushGroup("program_options");
	writer->writeBool("KAHAN_SUMMATION", opts.kahan);
	writer->writeBool("IGNORE_CSTR", opts.ignoreCSTR);
#ifdef ENABLE_BESSEL_TRUNCATION
	writer->writeBool("BESSEL_FUNCTION_TRUNCATION", true);
#else
	writer->writeBool("BESSEL_FUNCTION_TRUNCATION", false);
#endif
	{
		std::ostringstream oss;
		oss << mpfr::digits2bits(opts.precision);
		writer->writeString("PRECISION", oss.str());
		oss.str("");
		//writer->writeDouble("OUTPUT_PRECISION", std::numeric_limits<double>::digits10 + 1);

		writer->writeBool("EXTRAPOLATION", false);
		writer->writeInt("MAX_LAPLACE_SUMMANDS", opts.laplaceSummands);
		writer->writeInt("MAX_HANKEL_SUMMANDS", opts.hankelSummands);
		oss << opts.abscissa;
		writer->writeString("ABSCISSA", oss.str());
		oss.str("");
		oss << opts.errorCons;
		writer->writeString("CONSISTENCY_ERROR", oss.str());
		oss.str("");
		oss << opts.errorTrunc;
		writer->writeString("TRUNCATION_ERROR", oss.str());
		oss.str("");
		oss << opts.errorTrunc + opts.errorCons;
		writer->writeString("ERROR", oss.str());
		oss.str("");
		oss << opts.error;
		writer->writeString("REQUESTED_ERROR", oss.str());
		oss.str("");
		oss << opts.errorWeight;
		writer->writeString("ERROR_WEIGHT", oss.str());
	}
	writer->popGroup();

	writer->deleteGroup("output");
	writer->pushGroup("output");

	writer->pushGroup("solution");

	const int timePoints = solutionTimes.size() - timeOffset;
	std::vector<double> solutionBuffer(timePoints);

	const int numUnits = sys->numModels();
	const int numOutlets = sys->numOutlets();
	int numWrittenOutlets = 0;

	for (int i = 0; i < timePoints; i++)
		solutionBuffer[i] = static_cast<double>(solutionTimes[i + timeOffset]);

	writer->writeVectorDouble("SOLUTION_TIMES", timePoints, &solutionBuffer[0], 1, timePoints);

	std::ostringstream oss;

	for (int j = 0; j < numUnits; ++j)
	{
		oss << std::setw(3) << std::setfill('0') << j;
		writer->pushGroup("unit_" + oss.str());
		oss.str("");

		const int nPorts = std::max(1, sys->unitOperation(j)->numOutletPorts());

		// always split ports
		for (int k = 0; k < nPorts; ++k)
		{
			for (int i = 0; i < timePoints; i++)
				solutionBuffer[i] = static_cast<double>(output(numWrittenOutlets, i));

			oss << std::setw(3) << std::setfill('0') << k;
			const std::string outName = (single_as_multi_port || nPorts > 1) ? "SOLUTION_OUTLET_PORT_" + oss.str() : "SOLUTION_OUTLET";

			writer->writeVectorDouble(outName, timePoints, &solutionBuffer[0], 1, timePoints);
			oss.str("");

			numWrittenOutlets++;
		}
		writer->popGroup();
	}

	writer->closeFile();
}
#endif

void writeMeta(std::ostream& os, const casema::model::ModelSystem& model, const ProgramOptions& opts, const double simTime)
{
	os << "# CASEMA Version: " << casema::getVersion() << "\n";
	os << "# CASEMA Branch: " << casema::getBranchRefspec() << "\n";
	os << "# CASEMA Commit: " << casema::getCommitHash() << "\n";
	os << "# CASEMA File Format: " << casema::getFileFormat() << "\n";
	if(simTime > 0.0)
		os << "# Compute time: " << simTime << "\n";

	if (opts.kahan)
		os << "# Kahan summation: On\n";
	else
		os << "# Kahan summation: Off\n";
	if (opts.ignoreCSTR)
		os << "# Ignore CSTR: On\n";
	else
		os << "# Ignore CSTR: Off\n";

#ifdef ENABLE_BESSEL_TRUNCATION
	os << "# Bessel function truncation: On\n";
#else
	os << "# Bessel function truncation: Off\n";
#endif

	os << "# Model: " << opts.inFile << "\n";
	for (int i = 0; i < model.numModels(); ++i)
	{
		casema::model::UnitOperation const* const m = model.unitOperation(i);
		os << "#   Unit " << i << ": " << m->unitOperationName() << "\n";
	}
	os << "# Precision " << opts.precision << " digits = " << mpfr::digits2bits(opts.precision) << " bits\n";
	os << "# Extrapolation: Disabled\n";
	os << "# MaxSummandsLaplace: " << opts.laplaceSummands << "\n";
	os << "# MaxSummandsHankel: " << opts.hankelSummands << "\n";
	os << "# Abscissa (a): " << opts.abscissa << "\n";
	os << "# Consistency error: " << opts.errorCons << "\n";
	os << "# Truncation error: " << opts.errorTrunc << "\n";
	os << "# Error: " << opts.errorTrunc + opts.errorCons << " (requested " << opts.error << ")" << "\n";
	os << "# Error weight: " << opts.errorWeight << std::endl;
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpreal>& time, const casema::util::SlicedVector<mpfr::mpreal>& res, std::size_t precision, std::size_t timeOffset, const casema::model::ModelSystem& model)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	const std::size_t numUnits = model.numModels();
	const std::size_t numOutlets = model.numOutlets();

	fs << "time";
	for (std::size_t j = 0; j < numUnits; ++j)
	{
		casema::model::UnitOperation const* const m = model.unitOperation(j);
		for (int k = 0; k < std::max(1, m->numOutletPorts()); ++k)
			fs << ",UNIT" << j << "PORT" << k;
	}
	fs << "\n";

	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < time.size() - timeOffset; ++i)
	{
		fs << time[i+timeOffset];
		for (std::size_t j = 0; j < numOutlets; ++j)
			fs << "," << res(j, i);
		fs << "\n"; 
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


casema::util::SlicedVector<mpfr::mpreal> invert(const casema::model::ModelSystem& model, const mpfr::mpreal& maxTime, const std::vector<mpfr::mpreal>& outletTimes, const std::size_t timeOffset, const ProgramOptions& opts)
{
	ProgressBarUpdater pbu;

#ifdef _OPENMP
	std::vector<casema::model::ModelSystem::Workspace> ws(0);
	ws.reserve(opts.numThreads);
	for (int i = 0; i < opts.numThreads; ++i)
		ws.push_back(model.makeWorkspace());
#else
	casema::model::ModelSystem::Workspace ws = model.makeWorkspace();
#endif

	if (opts.kahan)
		return casema::invertLaplaceKahan(
			[&model, &ws](const mpfr::mpcomplex& s, mpfr::mpcomplex* res)
			{
#ifdef _OPENMP
				model.evaluate(s, res, ws[omp_get_thread_num()]);
#else
				model.evaluate(s, res, ws);
#endif
			},
			model.numOutlets(),
			opts.laplaceSummands,
			opts.precision,
			maxTime,
			opts.abscissa,
			outletTimes.data() + timeOffset,
			outletTimes.size() - timeOffset,
			[&](double progress, int phase)
			{
				pbu.update(progress, phase);
			}
		);
	else
		return casema::invertLaplace(
			[&model, &ws](const mpfr::mpcomplex& s, mpfr::mpcomplex* res)
			{
#ifdef _OPENMP
				model.evaluate(s, res, ws[omp_get_thread_num()]);
#else
				model.evaluate(s, res, ws);
#endif
			},
			model.numOutlets(),
			opts.laplaceSummands,
			opts.precision,
			maxTime,
			opts.abscissa,
			outletTimes.data() + timeOffset,
			outletTimes.size() - timeOffset,
			[&](double progress, int phase)
			{
				pbu.update(progress, phase);
			}
		);
}


int main(int argc, char** argv)
{
#ifdef _OPENMP
	Eigen::initParallel();
#endif

	ProgramOptions opts;

	try
	{
		TCLAP::CustomOutput customOut("chrom");
		TCLAP::CmdLine cmd("Uses a numerical inverse Laplace transform to solve GRM models", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<int>("t", "threads", "Number of threads (default: all available)", false, -1, "Int"))->storeIn(&opts.numThreads);

		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: same as working precision)", false, 0, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Working precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", 
													false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);

		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("a", "abscissa", "Abscissa in Durbin's method, used as safety margin if error (-e) is given", false, 0, "Float"))->storeIn(&opts.abscissa);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "error", "Error threshold (default: 1e-10)", false, 1e-10, "Float"))->storeIn(&opts.error);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("w", "weight", "Weight used to distribute error onto consistency and truncation (default: 0.5)", false, 0.5, "Float"))->storeIn(&opts.errorWeight);
		cmd >> (new TCLAP::ValueArg<std::size_t>("N", "lapsum", "Maximum number of (Laplace) summands in Durbin's method (Laplace inversion)", false, 0, "Int"))->storeIn(&opts.laplaceSummands);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "hankelsum", "Number of (Hankel) summands in Dini's expansion (Hankel inversion)", false, 0, "Int"))->storeIn(&opts.hankelSummands);
		cmd >> (new TCLAP::SwitchArg("", "kahan", "Use Kahan summation"))->storeIn(&opts.kahan);
		cmd >> (new TCLAP::SwitchArg("", "ignorecstr", "Ignore CSTRs in error estimates"))->storeIn(&opts.ignoreCSTR);

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
#ifdef _OPENMP
	if (opts.numThreads <= 0)
		opts.numThreads = omp_get_max_threads();
	else
		opts.numThreads = std::min(opts.numThreads, omp_get_max_threads());

	omp_set_num_threads(opts.numThreads);
	std::cout << "Using " << opts.numThreads << " threads\n";
#endif

	// Set default output precision
	if (opts.outPrecision == 0)
		opts.outPrecision = opts.precision;

	// Extract suffix
	const std::size_t fileEnding = opts.inFile.find_last_of('.') + 1;
	std::string ext = opts.inFile.substr(fileEnding);
	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	const std::unique_ptr<casema::io::IFileReader> rd(casema::io::createReader(ext));
	if (!rd)
	{
		std::cout << "ERROR: Input file format ." << ext << " not recognized" << std::endl;
		::mpfr_free_cache();
		return 1;
	}

	std::cout << "Reading model... " << std::flush;

	try
	{
		rd->openFile(opts.inFile, "r");
	}
	catch (const std::exception& e)
	{
		std::cout << "\nERROR: " << e.what() << std::endl;
		::mpfr_free_cache();
		return 1;
	}

	casema::io::FileParameterProvider fpp(*rd);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* sys;
	casema::SimulationTime st;

	try
	{
		fpp.pushScope("model");
		sys = mb.createSystem(fpp);
		fpp.popScope();

		st = casema::readTimes(fpp);
		sys->setSectionTimes(st.sectionTimes.data(), st.sectionTimes.size());
	}
	catch (const std::exception& e)
	{
		std::cout << "\nERROR: " << e.what() << std::endl;
		rd->closeFile();
		::mpfr_free_cache();
		return 1;
	}

	rd->closeFile();

	std::cout << "[DONE]" << std::endl;

	if (st.solutionTimes.empty())
	{
		std::cout << "Generating output time points... " << std::flush;
		st.solutionTimes.reserve(st.sectionTimes.back().toLong() + 1);

		mpfr::mpreal cur = mpfr::mpreal(1);
		while (cur <= st.sectionTimes.back())
		{
			st.solutionTimes.push_back(cur);
			cur += mpfr::mpreal(1);
		}

		std::cout << "[DONE]" << std::endl;
	}

	std::size_t timeOffset = 0;
	if (st.solutionTimes[0] <= mpfr::mpreal(0))
	{
		// Skip first time point
		std::cout << "Removing time point t = 0.0 since Laplace solution has a pole there" << std::endl;
		timeOffset = 1;
	}

	const mpfr::mpreal maxTime = casema::maxSimulationTime(st);

	// If summands and abscissa are given, use them
	if ((opts.abscissa > 0) && (opts.laplaceSummands > 0))
	{
		if (sys->hasValidEstimate())
		{
			opts.errorCons = casema::consistencyError(*sys, opts.abscissa, maxTime);
			opts.errorTrunc = casema::truncationError(*sys, opts.abscissa, maxTime, opts.laplaceSummands);
		}
		else
		{
			opts.errorCons = 0.0;
			opts.errorTrunc = 0.0;
		}

		opts.error = 0;
		opts.errorWeight = 0.5;
	}
	else if ((opts.errorWeight > 0) && (opts.errorWeight < 1) && (opts.error > 0))
	{
		if (!sys->hasValidEstimate())
		{
			std::cout << "WARNING: Model does not support error estimates, use abscissa and number of summands instead\n";
			std::cout << "WARNING: Consistency and truncation error estimates may be incorrect\n";
		}

		std::size_t s = 0;
		std::vector<mpfr::mpreal> nTerms(0);
		if (opts.abscissa > 0)
			nTerms = std::move(casema::abscissaSummandsFromError(*sys, opts.abscissa, maxTime, opts.errorWeight, opts.error, opts.ignoreCSTR, opts.abscissa, s));
		else
			nTerms = std::move(casema::abscissaSummandsFromError(*sys, mpfr::mpreal(0), maxTime, opts.errorWeight, opts.error, opts.ignoreCSTR, opts.abscissa, s));

		if (opts.laplaceSummands > 0)
			std::cout << "Overruling automatically determined number of summands " << s << " by " << opts.laplaceSummands << "\n";
		else
			opts.laplaceSummands = s;

		if ((opts.abscissa <= 0.0) || (opts.laplaceSummands == static_cast<std::size_t>(-1)))
		{
			std::cout << "ERROR: Incorrect abscissa (" << opts.abscissa << ") or number of summands (" << opts.laplaceSummands << "), please specify manually" << std::endl;
			::mpfr_free_cache();
			return 1;
		}

		for (int i = 0; i < nTerms.size(); ++i)
			std::cout << " -> Unit " << i << " (" << sys->unitOperation(i)->unitOperationName() << "): " << nTerms[i] << " summands\n";

		opts.errorCons = casema::consistencyError(*sys, opts.abscissa, maxTime);
		opts.errorTrunc = casema::truncationError(*sys, opts.abscissa, maxTime, opts.laplaceSummands);
	}
	else
	{
		std::cout << "ERROR: Specify either abscissa and Laplace summands or error and weight" << std::endl;
		::mpfr_free_cache();
		return 1;
	}

	std::vector<mpfr::mpreal> besselZeros(0);
	if (sys->needsBesselZeros())
	{
		if (opts.hankelSummands == 0)
		{
			std::cout << "ERROR: Please specify positive amount of Hankel summands" << std::endl;
			::mpfr_free_cache();
			return 1;
		}

		std::cout << "Computing " << opts.hankelSummands << " zeros of Bessel function... " << std::flush;

		besselZeros.resize(opts.hankelSummands);
		casema::besselZerosJ1(opts.hankelSummands, besselZeros.data());

		sys->setBesselZeros(besselZeros.data(), besselZeros.size());

		std::cout << "[DONE]" << std::endl;
	}

	writeMeta(std::cout, *sys, opts, -1.0);

	std::cout << "Starting Laplace inversion" << std::endl;

	const auto start = std::chrono::high_resolution_clock::now();

	const casema::util::SlicedVector<mpfr::mpreal> output = invert(*sys, maxTime, st.solutionTimes, timeOffset, opts);

	const double time_sim = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start).count();

	std::cout << "Writing output" << std::endl;

#ifdef ENABLE_HDF5
	// Check if input file is h5 file and if so, write to that. Otherwise, print to console
	bool inIsH5 = false;
	const std::size_t dotPos = opts.inFile.find_last_of('.');

	if (dotPos != std::string::npos && dotPos != 0)
	{
		if ((opts.inFile.substr(dotPos + 1) == "h5"))
		{
			inIsH5 = true;

			writeMetaAndResultToH5(opts, sys, st.solutionTimes, output, time_sim, timeOffset);
		}
	}
	if (!inIsH5)
#endif
	{
		if (opts.outFile.empty())
		{
			writeResult(std::cout, st.solutionTimes, output, opts.outPrecision, timeOffset, *sys);
		}
		else
		{
			std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);

			writeMeta(fs, *sys, opts, time_sim);
			writeResult(fs, st.solutionTimes, output, opts.outPrecision, timeOffset, *sys);
		}
	}

	::mpfr_free_cache();

	return 0;
}

