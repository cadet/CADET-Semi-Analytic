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

#include "LaplaceError.hpp"
#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"

#include <algorithm>
#include <cstring>

namespace
{
	inline mpfr::mpreal maxTimeToT(const mpfr::mpreal& maxTime) CASEMA_NOEXCEPT { return maxTime / 2 * mpfr::mpreal("1.01"); }
}

namespace casema
{

mpfr::mpreal consistencyError(const model::ModelSystem& sys, const mpfr::mpreal& abscissa, const mpfr::mpreal& maxTime) CASEMA_NOEXCEPT
{
	const mpfr::mpreal T = maxTimeToT(maxTime);
	return sys.timeDomainUpperBound() / expm1(2 * T * abscissa);
}

mpfr::mpreal truncationError(const model::ModelSystem& sys, const mpfr::mpreal& abscissa, const mpfr::mpreal& maxTime, std::size_t nSummands) CASEMA_NOEXCEPT
{
	const mpfr::mpreal T = maxTimeToT(maxTime);
	return exp(2 * T * abscissa) / T * sys.truncationError(abscissa, T, nSummands);
}

std::vector<mpfr::mpreal> abscissaSummandsFromError(const model::ModelSystem& sys, const mpfr::mpreal& abscissaSafety, const mpfr::mpreal& maxTime, const mpfr::mpreal& errorWeight, const mpfr::mpreal& error, bool ignoreCSTRs, mpfr::mpreal& abscissa, std::size_t& laplaceSummands) CASEMA_NOEXCEPT
{
	using std::max;

	const mpfr::mpreal T = maxTimeToT(maxTime);
	abscissa = log1p(sys.timeDomainUpperBound() / (errorWeight * error)) / (2 * T) + abscissaSafety;

	const std::vector<mpfr::mpreal> nTerms = sys.inverseTruncationErrorOfUnits(abscissa, T, (1 - errorWeight) * error * T / exp(2 * abscissa * T));
	
	if (ignoreCSTRs)
	{
		laplaceSummands = 0;
		for (int i = 0; i < sys.numModels(); ++i)
		{
			if (std::strcmp(sys.unitOperation(i)->unitOperationName(), "CSTR") == 0)
				continue;

			const std::size_t m = ceil(nTerms[i]).toLong();
			laplaceSummands = max(laplaceSummands, m);
		}
	}
	else
		laplaceSummands = ceil(*std::max_element(nTerms.begin(), nTerms.end())).toLong();

	return nTerms;
}

}
