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

#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPAD.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationTestFuncs.hpp"

template <typename real_t, class Seq, typename F>
void run(std::size_t len, F accMethod)
{
	// Build original sequence
	Seq f;
	const real_t trueVal = f.solution();
	std::vector<real_t> seq = f.sequence(len);

	// Accelerate
	std::vector<real_t> acc = accMethod(seq);

	// Output in table form
	std::cout << "========= " << f.name() << " ==========" << std::endl;
	std::cout << " i       seq[i]           error           acc[i]           error" << std::endl;

	for (std::size_t i = 0; i < acc.size(); ++i)
	{
		std::cout << i << " " << seq[i] << " "  << abs(trueVal - seq[i]) << " " << acc[i] << " " << abs(trueVal - acc[i]) << std::endl;
	}

	for (std::size_t i = acc.size(); i < seq.size(); ++i)
	{
		std::cout << i << " " << seq[i] << " "  << abs(trueVal - seq[i]) << std::endl;
	}
}


template <typename real_t, typename F>
void runAll(std::size_t len, F accMethod, const char* const accName)
{
	std::cout << "\n######## " << accName << " ########\n" << std::endl;

	run<real_t, EulerE<real_t> >(len, accMethod);
	run<real_t, SquareRootOfTwo<real_t> >(len, accMethod);
	run<real_t, RootOfCubicPolynomial<real_t> >(len, accMethod);
	run<real_t, RootOfExpEquation<real_t> >(len, accMethod);
	run<real_t, RootOfQuadraticPolynomial<real_t> >(len, accMethod);
	run<real_t, RiemannZeta<real_t> >(len, accMethod);
	run<real_t, AlternatingHarmonicSeries<real_t> >(len, accMethod);
}

int main(int argc, char** argv)
{
	const int precision = 50;
	const int seqLength = 10;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
	typedef mpfr::mpreal real_t;

	std::cout.flags(std::ios::scientific);
	std::cout.precision(10);

	runAll<real_t>(seqLength, casema::aitkenDeltaSquared<real_t>, "Aitken's Delta Squared");
	runAll<real_t>(seqLength, casema::wynnEpsilon<real_t>, "Wynn's Epsilon Algorithm");
	runAll<real_t>(seqLength, casema::wynnRho<real_t>, "Wynn's Rho Algorithm");
	runAll<real_t>(seqLength, casema::iteratedAitkenDeltaSquared<real_t>, "Iterated Aitken's Delta Squared");
	runAll<real_t>(seqLength, casema::levinU<real_t>, "Levin U Algorithm");
	runAll<real_t>(seqLength, casema::levinT<real_t>, "Levin T Algorithm");
	runAll<real_t>(seqLength, casema::iteratedBrezinskiTheta<real_t>, "Iterated Brezinski Theta Method");
	runAll<real_t>(seqLength, casema::brezinskiTheta<real_t>, "Brezinski Theta Method");
	runAll<real_t>(seqLength, casema::nevilleAitken<real_t>, "Neville Aitken Table");
	runAll<real_t>(seqLength, casema::richardson<real_t>, "Richardson Extrapolation");

	::mpfr_free_cache();

	return 0;
}

