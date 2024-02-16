// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-2020: Samuel Leweke¹²
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#ifdef _OPENMP
	#include <omp.h>
#endif

int main(int argc, char* argv[])
{
#ifdef _OPENMP
	omp_set_num_threads(1);
#endif

	Catch::Session session;

	// Parse the command line and check for error
	const int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0)
		return returnCode;

	// Run tests
	const int result = session.run();
	return result;
}
