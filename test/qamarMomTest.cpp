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

#include "GRMLaplaceSolution.hpp"
#include "LaplaceInlet.hpp"

#include "AnalyticalMoments.hpp"

#include "MPAD.hpp"

#ifdef CASEMA_USE_FADBAD
	#include <tadiff.h>
#endif
#include <cppad/cppad.hpp>

typedef casema::laplaceSolution::Inlet<mpfr::mpreal, mpfr::mpreal> Inlet_t;


casema::ModelData<mpfr::mpreal> qamarModel(const mpfr::mpreal& velocity)
{
	casema::ModelData<mpfr::mpreal> data;

	data.nComponents = 1;
	
	data.initialLiquidConcentration.push_back(mpfr::mpreal(0));
	data.initialSolidConcentration.push_back(mpfr::mpreal(0));

	data.colDispersion = mpfr::mpreal("0.002") / mpfr::mpreal(600000);
	data.colLength = mpfr::mpreal("0.017");
	data.colPorosity = mpfr::mpreal("0.4");

	data.parRadius = mpfr::mpreal("0.00004");
	data.parPorosity = mpfr::mpreal("0.333");

	data.filmDiffusion.push_back(mpfr::mpreal("0.01") / mpfr::mpreal(6000));
	data.particleDiffusion.push_back(mpfr::mpreal("1e-6") / (mpfr::mpreal(600000)*data.parPorosity));
	data.surfaceDiffusion.push_back(mpfr::mpreal(0));

	data.velocity = velocity;

	data.kineticBinding = false;

	data.bindingModel = casema::LINEAR;
	data.linearKA.push_back(mpfr::mpreal(25));
	data.linearKD.push_back(mpfr::mpreal(10));

	data.nInletSections = 2;
	data.sectionTimes.push_back(mpfr::mpreal(0));
	data.sectionTimes.push_back(mpfr::mpreal(1200));
	data.sectionTimes.push_back(mpfr::mpreal(10000));
	data.constCoeff.push_back(mpfr::mpreal(1000)); // 1 g / l in mol / m^3
	data.constCoeff.push_back(mpfr::mpreal(0));

	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));
	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));

	return data;
}

casema::ModelData<mpfr::mpreal> qamarModelCmMin(const mpfr::mpreal& velocity)
{
	casema::ModelData<mpfr::mpreal> data;

	data.nComponents = 1;
	
	data.initialLiquidConcentration.push_back(mpfr::mpreal(0));
	data.initialSolidConcentration.push_back(mpfr::mpreal(0));

	data.colDispersion = mpfr::mpreal("0.002");
	data.colLength = mpfr::mpreal("1.7");
	data.colPorosity = mpfr::mpreal("0.4");

	data.parRadius = mpfr::mpreal("0.004");
	data.parPorosity = mpfr::mpreal("0.333");

	data.filmDiffusion.push_back(mpfr::mpreal("0.01"));
	data.particleDiffusion.push_back(mpfr::mpreal("1e-6") / data.parPorosity);
	data.surfaceDiffusion.push_back(mpfr::mpreal(0));

	data.velocity = velocity;

	data.kineticBinding = false;

	data.bindingModel = casema::LINEAR;
	data.linearKA.push_back(mpfr::mpreal(25));
	data.linearKD.push_back(mpfr::mpreal(10));

	data.nInletSections = 2;
	data.sectionTimes.push_back(mpfr::mpreal(0));
	data.sectionTimes.push_back(mpfr::mpreal(20));
	data.sectionTimes.push_back(mpfr::mpreal(10000));
	data.constCoeff.push_back(mpfr::mpreal("0.001")); // 1 g / l in mol / cm^3
	data.constCoeff.push_back(mpfr::mpreal(0));

	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));
	data.linCoeff.push_back(mpfr::mpreal(0));
	data.quadCoeff.push_back(mpfr::mpreal(0));
	data.cubicCoeff.push_back(mpfr::mpreal(0));

	return data;
}


void writeResult(std::ostream& fs, const std::vector<mpfr::mpreal>& mom, std::size_t precision)
{
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "moment,outlet\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	for (std::size_t i = 0; i < mom.size(); ++i)
	{
		fs << i << "," << mom[i] << "\n"; 
	}

	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}

#ifdef CASEMA_USE_FADBAD

	template <class Model_t, class Solution_t, class Inlet_t>
	std::vector<mpfr::mpreal> momentsFADBAD(const Model_t& model, const Inlet_t& inlet, std::size_t order, std::size_t precision, const mpfr::mpreal& evalPoint)
	{
		Solution_t solution(model, inlet);
		std::cout << solution.name() << std::endl;

		// Taylor expansion
		fadbad::T<mpfr::mpreal> x;
		x = evalPoint;
		x[1] = 1;

		fadbad::T<mpfr::mpreal> val = solution(x);
		val.eval(order);

		std::cout.flags(std::ios::scientific);
		std::cout.precision(precision);

		std::vector<mpfr::mpreal> output(order+1);

		std::cout << "Taylor expansion at " << x[0] << " yields " << val[0] << ":\n";
		mpfr::mpreal fact(1);
		for (std::size_t i = 0; i < order+1; ++i)
		{
			if (i > 1)
				fact *= i;

			output[i] = fact * val[i];
			std::cout << "coeff[i] = " << val[i] << std::endl;
			std::cout << "d^" << i << " f / dx^" << i << " = " << output[i] << std::endl;;
		}
		return output;
	}

#endif

template <class Model_t, class Solution_t, class Inlet_t>
std::vector<mpfr::mpreal> momentsCppAD(const Model_t& model, const Inlet_t& inlet, std::size_t order, std::size_t precision, const mpfr::mpreal& evalPoint)
{
	Solution_t solution(model, inlet);
	std::cout << solution.name() << std::endl;

	std::vector< CppAD::AD<mpfr::mpreal> > x(1);
	x[0] = evalPoint;

	// Start recording
	CppAD::Independent(x);
	// Execute function
	std::vector< CppAD::AD<mpfr::mpreal> > y(1);
	y[0] = solution(x[0]);
	// Stop recording
	CppAD::ADFun<mpfr::mpreal> f(x, y);

	// Where to evaluate
	std::vector<mpfr::mpreal> pos(order+1);
	pos[0] = evalPoint;
	pos[1] = mpfr::mpreal(1);
	for (std::size_t i = 2; i < order+1; ++i)
		pos[i] = mpfr::mpreal(0);

	const std::vector<mpfr::mpreal> dfdx = f.Forward(order, pos);

	std::cout.flags(std::ios::scientific);
	std::cout.precision(precision);

	std::vector<mpfr::mpreal> output(order+1);

	std::cout << "Taylor expansion at " << x[0] << " yields " << dfdx[0] << ":\n";
	mpfr::mpreal fact(1);
	for (std::size_t i = 0; i < order+1; ++i)
	{
		if (i > 1)
			fact *= i;

		output[i] = fact * dfdx[i];
		std::cout << "coeff[i] = " << dfdx[i] << std::endl;
		std::cout << "d^" << i << " f / dx^" << i << " = " << output[i] << std::endl;;
	}
	return output;
}


void run(const casema::ModelData<mpfr::mpreal>& model, int order, std::size_t precision, const mpfr::mpreal& evalPoint, bool useCppAD, const std::string& outFile)
{
	std::cout << "Calculate derivatives of Laplace domain solution";

	std::vector<mpfr::mpreal> output;
	Inlet_t inlet(model);
	if (model.kineticBinding)
	{
		typedef casema::laplaceSolution::SingleComponentLinearDynamic<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;

		if (!useCppAD)
		{
#ifdef CASEMA_USE_FADBAD
			output = momentsFADBAD<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, order, precision, evalPoint);
#endif            
		}
		else
			output = momentsCppAD<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, order, precision, evalPoint);
	}
	else
	{
		typedef casema::laplaceSolution::SingleComponentLinearRapidEquilibrium<mpfr::mpreal, mpfr::mpreal, Inlet_t> Solution_t;

		if (!useCppAD)
		{
#ifdef CASEMA_USE_FADBAD
			output = momentsFADBAD<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, order, precision, evalPoint);
#endif            
		}
		else
			output = momentsCppAD<casema::ModelData<mpfr::mpreal>, Solution_t, Inlet_t>(model, inlet, order, precision, evalPoint);
	}

	std::cout << "=========== MOMENTS (analytical) ============" << std::endl; 

	const std::vector<mpfr::mpreal> moments = casema::analyticalMoments(model);

	std::cout.flags(std::ios::scientific);
	std::cout.precision(precision);
	for (std::size_t i = 0; i < moments.size(); ++i)
	{
		std::cout << "mu_" << i << " = " << moments[i] << std::endl;
	}

	std::cout << "\n";
	std::cout << "=========== MOMENTS (numerical) ============" << std::endl; 
	std::cout << "MU_0 = " << output[0] << std::endl;

	for (std::size_t i = 1; i < output.size(); ++i)
	{
		if (i % 2 == 0)
			std::cout << "MU_" << i << " = " << output[i] / output[0] << std::endl;
		else
			std::cout << "MU_" << i << " = " << -output[i] / output[0] << std::endl;
	}

	std::cout << "Writing output" << std::endl;
	if (outFile.empty())
	{
		writeResult(std::cout, output, precision);
	}
	else
	{
		std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);
		writeResult(fs, output, precision);
	}
}


int main(int argc, char** argv)
{

	std::size_t precision;
	int order;
	mpfr::mpreal velocity;
	mpfr::mpreal evalPoint;
	std::string outFile;
	bool useCppAD;
	bool useBaseSI;

	try
	{
		TCLAP::CmdLine cmd("Uses Laplace transform and algorithmic differentiation to compute moments of GRM models given in Qamar et al. (2014)", ' ', "0.0");

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision in base 10 digits (default: " + std::to_string(mpfr::bits2digits(mpfr::mpreal::get_default_prec())) + ", stay below " + std::to_string(mpfr::bits2digits(MPFR_PREC_MAX)) + ")", false, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<int>("n", "order", "Maximum moment order (default: 10)", false, 10, "Int"))->storeIn(&order);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write full precision output to file (default: disabled)", false, std::string(), "File"))->storeIn(&outFile);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("e", "eval", "Evaluation point (default: 1e-30)", false, 1.0e-30, "Double"))->storeIn(&evalPoint);
		cmd >> (new TCLAP::ValueArg<mpfr::mpreal>("v", "velocity", "Interstitial velocity in cm/min (default: 0.4 cm/min)", false, 0.4, "Double"))->storeIn(&velocity);
		cmd >> (new TCLAP::SwitchArg("s", "SI", "Use SI base units (default: cm,min)"))->storeIn(&useBaseSI);

#ifdef CASEMA_USE_FADBAD
		cmd >> (new TCLAP::SwitchArg("c", "cppAD", "Use CppAD (default: FADBAD++)"))->storeIn(&useCppAD);
#else
		useCppAD = true;
#endif

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	if (useBaseSI)
	{
		run(qamarModel(velocity / mpfr::mpreal(6000)), order, precision, evalPoint, useCppAD, outFile);
	}
	else
	{
		run(qamarModelCmMin(velocity), order, precision, evalPoint, useCppAD, outFile);
	}

	::mpfr_free_cache();

	return 0;
}

