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
#include <cmath>
#include <limits>

#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPAD.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"
#include "ExtrapolationTestFuncs.hpp"

template <typename real_t, class Seq>
void run(std::size_t maxLen, casema::ConsensusEstimator<real_t>& accMethod)
{
	accMethod.reset();

	// Build original sequence
	Seq f;
	const real_t trueVal = f.solution();

	// Accelerate
	std::size_t n = 0;
	real_t estLimit(0);
	real_t lastElement(0);
	while (n < maxLen)
	{
		lastElement = f.next();
		estLimit = accMethod.next(lastElement);
		++n;
	}

	// Output in table form
	std::cout << "========= " << f.name() << " ==========" << std::endl;
	std::cout << "Loops " << n << " calls = " << accMethod.iterations() << " iter" << std::endl;
	std::cout << "True value:   " << trueVal << std::endl;
	std::cout << "Last element: " << lastElement << " => Error " << abs(trueVal - lastElement) << std::endl;
	std::cout << "Est. Limit:   " << estLimit << " => Error " << abs(trueVal - estLimit) << std::endl;
	std::cout << "Radius:       " << accMethod.radius() << std::endl;
}


template <typename real_t>
void runAll(std::size_t maxLen, casema::ConsensusEstimator<real_t>& accMethod)
{
	std::cout << "\n######## " << accMethod.name() << " ########\n" << std::endl;

	run<real_t, EulerE<real_t>>(maxLen, accMethod);
	run<real_t, SquareRootOfTwo<real_t>>(maxLen, accMethod);
	run<real_t, RootOfCubicPolynomial<real_t>>(maxLen, accMethod);
	run<real_t, RootOfExpEquation<real_t>>(maxLen, accMethod);
	run<real_t, RootOfQuadraticPolynomial<real_t>>(maxLen, accMethod);
	run<real_t, RiemannZeta<real_t>>(maxLen, accMethod);
	run<real_t, AlternatingHarmonicSeries<real_t>>(maxLen, accMethod);
}

template <typename real_t, template <class T> class F>
void addEstimator(casema::ConsensusEstimator<real_t>& consEst)
{
	casema::ConvergenceMonitor<real_t, F>* convMon = new casema::ConvergenceMonitor<real_t, F>(F<real_t>());
	convMon->times(3);
	convMon->thresholdAbs(10 * std::numeric_limits<real_t>::epsilon());

	consEst.add(convMon);
}

int main(int argc, char** argv)
{
	const int precision = 50;
	const int maxIter = 100;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
	typedef mpfr::mpreal real_t;

	std::cout.flags(std::ios::scientific);
	std::cout.precision(10);

	casema::ConsensusEstimator<real_t> consEst;
	consEst.mustAgree(3);

	addEstimator<real_t, casema::AitkenDeltaSquaredMethod>(consEst);
	addEstimator<real_t, casema::WynnEpsilonMethod>(consEst);
	addEstimator<real_t, casema::WynnRhoMethod>(consEst);
	addEstimator<real_t, casema::IteratedAitkenDeltaSquaredMethod>(consEst);
	addEstimator<real_t, casema::LevinUMethod>(consEst);
	addEstimator<real_t, casema::LevinTMethod>(consEst);
	addEstimator<real_t, casema::IteratedBrezinskiThetaMethod>(consEst);
	addEstimator<real_t, casema::BrezinskiThetaMethod>(consEst);
	addEstimator<real_t, casema::NevilleAitkenMethod>(consEst);
	addEstimator<real_t, casema::RichardsonMethod>(consEst);

	runAll(maxIter, consEst);

	::mpfr_free_cache();

	return 0;
}
