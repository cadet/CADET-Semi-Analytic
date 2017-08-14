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

#include "LaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "ADMoments.hpp"
#include "AnalyticalMoments.hpp"

typedef casema::LaplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;

struct ProgramOptions
{
	std::size_t precision;
	std::size_t outPrecision;
	int order;
	std::string inFile;
	std::string outFile;
	std::size_t nPoints;
	mpfr::mpreal start;
	mpfr::mpreal end;
	bool useCppAD;
	bool linearSeq;
};


void writeResult(std::ostream& fs, const std::vector<mpfr::mpreal>& pos, const std::vector<mpfr::mpreal>& mom, std::size_t precision)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "x,mom,momlog10\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < mom.size(); ++i)
	{
		fs << pos[i] << "," << mom[i] << "," << log10(abs(mom[i])) << "\n"; 
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


void calcLogLinearSpacing(casema::MomentGenerator<mpfr::mpreal>& momGen, const mpfr::mpreal& normalization, std::vector<mpfr::mpreal>& points, std::vector<mpfr::mpreal>& output, const ProgramOptions& opts)
{
	const mpfr::mpreal mpStart = log(opts.start);
	const mpfr::mpreal slope = (log(opts.end) - mpStart) / (opts.nPoints - 1);

	for (std::size_t i = 0; i < points.size(); ++i)
	{
		points[i] = exp(mpStart + slope * i);
		output[i] = momGen.moment(points[i], normalization, opts.order);
	}
}


void calcLinearSpacing(casema::MomentGenerator<mpfr::mpreal>& momGen, const mpfr::mpreal& normalization, std::vector<mpfr::mpreal>& points, std::vector<mpfr::mpreal>& output, const ProgramOptions& opts)
{
	const mpfr::mpreal slope = (opts.end - opts.start) / (opts.nPoints - 1);
	for (std::size_t i = 0; i < points.size(); ++i)
	{
		points[i] = opts.start + slope * i;
		output[i] = momGen.moment(points[i], normalization, opts.order);
	}
}


void run(casema::ModelData<mpfr::mpreal>& model, const ProgramOptions& opts)
{
	std::vector<mpfr::mpreal> points(opts.nPoints, opts.start);
	std::vector<mpfr::mpreal> output(opts.nPoints, 0);

	casema::MomentGenerator<mpfr::mpreal>* momGen = nullptr;
	Inlet_t inlet(model);

	switch (model.modelType)
	{
		case casema::GeneralRateModel:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::GeneralRateModel::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				
				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			else
			{
				typedef casema::LaplaceSolution::GeneralRateModel::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;

				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model,inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			break;
		case casema::LumpedRateModelWithPores:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithPores::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				
				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			else
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithPores::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;

				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model,inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			break;
		case casema::LumpedRateModelWithoutPores:
			if (model.kineticBinding)
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithoutPores::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;
				
				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			else
			{
				typedef casema::LaplaceSolution::LumpedRateModelWithoutPores::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;

				if (opts.useCppAD)
					momGen = new casema::CppADMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model,inlet));
				else
				{
#ifdef CASEMA_USE_FADBAD
					momGen = new casema::FadBadMomentGenerator<mpfr::mpreal, Solution_t>(Solution_t(model, inlet));
#endif
				}
			}
			break;
	}

	const mpfr::mpreal normalization = casema::detail::injectedMass(model);

	if (opts.linearSeq)
		calcLinearSpacing(*momGen, normalization, points, output, opts);
	else
		calcLogLinearSpacing(*momGen, normalization, points, output, opts);

	delete momGen;

	std::cout << "Writing output" << std::endl;
	if (opts.outFile.empty())
	{
		writeResult(std::cout, points, output, opts.outPrecision);
	}
	else
	{
		std::ofstream fs(opts.outFile, std::ofstream::out | std::ofstream::trunc);
		writeResult(fs, points, output, opts.outPrecision);
	}
}


int main(int argc, char** argv)
{
	ProgramOptions opts;

	try
	{
		TCLAP::CmdLine cmd("Uses Laplace transform to scan mu(s) of GRM models", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&opts.precision);
		cmd >> (new TCLAP::ValueArg<std::size_t>("P", "outprec", "Output precision (default: 16)", false, 16, "Int"))->storeIn(&opts.outPrecision);
		cmd >> (new TCLAP::ValueArg<int>("m", "order", "Moment order (default: 2)", false, 2, "Int"))->storeIn(&opts.order);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&opts.outFile);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("s", "start", "Start (default: 1e-20)", false, mpfr::mpreal(1e-20), "Double"))->storeIn(&opts.start);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "end", "End (default: 1e-30)", false, mpfr::mpreal(1e-30), "Double"))->storeIn(&opts.end);
		cmd >> (new TCLAP::ValueArg<std::size_t>("n", "points", "Points (default: 100)", false, 100, "Int"))->storeIn(&opts.nPoints);
#ifdef CASEMA_USE_FADBAD
		cmd >> (new TCLAP::SwitchArg("c", "cppAD", "Use CppAD (default: FADBAD++)"))->storeIn(&opts.useCppAD);
#else
		opts.useCppAD = true;
#endif
		cmd >> (new TCLAP::SwitchArg("l", "linear", "Use linear spacing (default: log-linear)"))->storeIn(&opts.linearSeq);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("model", "Model file", true, "", "File"))->storeIn(&opts.inFile);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(opts.precision));

	casema::ModelData<mpfr::mpreal> model;
	bool loaded = casema::readModel<mpfr::mpreal>(opts.inFile, model);

	if (!loaded)
		std::cout << "ERROR: Wrong input file format! Aborting..." << std::endl;
	else
	{
		if (casema::checkModelForCompatibility(model))        
			run(model, opts);
	}

	::mpfr_free_cache();

	return 0;
}
