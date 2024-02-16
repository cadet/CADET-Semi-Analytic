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

/**
 * @file 
 * Provides a driver for configuring and running a CASEMA simulation
 */

#ifndef CASEMA_DRIVER_HPP_
#define CASEMA_DRIVER_HPP_

#include "MPReal.hpp"

#include <vector>


namespace casema
{

	namespace io
	{
		class IParameterProvider;
	}


	struct SimulationTime
	{
		std::vector<mpfr::mpreal> sectionTimes;
		std::vector<mpfr::mpreal> solutionTimes;
	};

	SimulationTime readTimes(io::IParameterProvider& paramProvider);

	mpfr::mpreal maxSimulationTime(const SimulationTime& st);

} // namespace casema

#endif  // CASEMA_DRIVER_HPP_
