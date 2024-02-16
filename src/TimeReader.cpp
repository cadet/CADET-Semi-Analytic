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

#include "TimeReader.hpp"
#include "io/ParameterProvider.hpp"


namespace casema
{

SimulationTime readTimes(io::IParameterProvider& paramProvider)
{
	SimulationTime st;

	paramProvider.pushScope("solver");

	if (paramProvider.exists("USER_SOLUTION_TIMES"))
	{
		const std::vector<double> userTimes = paramProvider.getDoubleArray("USER_SOLUTION_TIMES");
		st.solutionTimes.reserve(userTimes.size());
		for (double t : userTimes)
			st.solutionTimes.push_back(t);
	}

	paramProvider.pushScope("sections");

	const std::vector<double> secTimes = paramProvider.getDoubleArray("SECTION_TIMES");
	st.sectionTimes.reserve(secTimes.size());
	for (double t : secTimes)
		st.sectionTimes.push_back(t);

	paramProvider.popScope(); // sections scope

	paramProvider.popScope();

	return st;
}

mpfr::mpreal maxSimulationTime(const SimulationTime& st)
{
	if (!st.solutionTimes.empty())
	{
		return max(st.solutionTimes.back(), st.sectionTimes.back());
	}
	return st.sectionTimes.back();
}

} // namespace casema
