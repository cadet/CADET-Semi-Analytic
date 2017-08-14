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
#include <cmath>
#include <limits>

#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPAD.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"
#include "ExtrapolationTestFuncs.hpp"

template <typename real_t, class Seq, template <class T> class F>
void run(std::size_t maxLen, F<real_t> accMethod)
{
	// Build original sequence
	Seq f;
	const real_t trueVal = f.solution();

	casema::ConvergenceMonitor<real_t, F> convMon(accMethod);
	convMon.times(3);
	convMon.thresholdAbs(10 * std::numeric_limits<real_t>::epsilon());

	// Accelerate
	std::size_t n = 0;
	real_t estLimit(0);
	real_t lastElement(0);
	while ((n < maxLen) && !convMon.converged())
	{
		lastElement = f.next();
		estLimit = convMon.next(lastElement);
		++n;
	}

	// Output in table form
	std::cout << "========= " << f.name() << " ==========" << std::endl;
	std::cout << "Converged: " << (convMon.converged() ? "YES" : "NO") << " -- " << n << " calls = " << convMon.iterations() << " iter" << std::endl;
	std::cout << "True value:   " << trueVal << std::endl;
	std::cout << "Last element: " << lastElement << " => Error " << abs(trueVal - lastElement) << std::endl;
	std::cout << "Est. Limit:   " << estLimit << " => Error " << abs(trueVal - estLimit) << std::endl;
}


template <typename real_t, template <class T> class F>
void runAll(std::size_t maxLen, F<real_t> accMethod)
{
	std::cout << "\n######## " << accMethod.name() << " ########\n" << std::endl;

	run<real_t, EulerE<real_t>, F>(maxLen, accMethod);
	run<real_t, SquareRootOfTwo<real_t>, F>(maxLen, accMethod);
	run<real_t, RootOfCubicPolynomial<real_t>, F>(maxLen, accMethod);
	run<real_t, RootOfExpEquation<real_t>, F>(maxLen, accMethod);
	run<real_t, RootOfQuadraticPolynomial<real_t>, F>(maxLen, accMethod);
	run<real_t, RiemannZeta<real_t>, F>(maxLen, accMethod);
	run<real_t, AlternatingHarmonicSeries<real_t>, F>(maxLen, accMethod);
}

int main(int argc, char** argv)
{
	const int precision = 50;
	const int maxIter = 100;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
	typedef mpfr::mpreal real_t;

	std::cout.flags(std::ios::scientific);
	std::cout.precision(10);

	runAll<real_t>(maxIter, casema::AitkenDeltaSquaredMethod<real_t>());
	runAll<real_t>(maxIter, casema::WynnEpsilonMethod<real_t>());
	runAll<real_t>(maxIter, casema::WynnRhoMethod<real_t>());
	runAll<real_t>(maxIter, casema::IteratedAitkenDeltaSquaredMethod<real_t>());
	runAll<real_t>(maxIter, casema::LevinUMethod<real_t>());
	runAll<real_t>(maxIter, casema::LevinTMethod<real_t>());
	runAll<real_t>(maxIter, casema::IteratedBrezinskiThetaMethod<real_t>());
	runAll<real_t>(maxIter, casema::BrezinskiThetaMethod<real_t>());
	runAll<real_t>(maxIter, casema::NevilleAitkenMethod<real_t>());
	runAll<real_t>(maxIter, casema::RichardsonMethod<real_t>());

	::mpfr_free_cache();

	return 0;
}

