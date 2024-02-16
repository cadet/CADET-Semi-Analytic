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

#ifndef CASEMA_LAPLACEERROR_HPP_
#define CASEMA_LAPLACEERROR_HPP_

#include "casemaCompilerInfo.hpp"
#include "MPReal.hpp"

#include <vector>

namespace casema
{
	namespace model
	{
		class ModelSystem;
	}

	mpfr::mpreal consistencyError(const model::ModelSystem& sys, const mpfr::mpreal& abscissa, const mpfr::mpreal& maxTime) CASEMA_NOEXCEPT;
	mpfr::mpreal truncationError(const model::ModelSystem& sys, const mpfr::mpreal& abscissa, const mpfr::mpreal& maxTime, std::size_t nSummands) CASEMA_NOEXCEPT;

	std::vector<mpfr::mpreal> abscissaSummandsFromError(const model::ModelSystem& sys, const mpfr::mpreal& abscissaSafety, const mpfr::mpreal& maxTime, const mpfr::mpreal& errorWeight, const mpfr::mpreal& error, bool ignoreCSTRs, mpfr::mpreal& abscissa, std::size_t& laplaceSummands) CASEMA_NOEXCEPT;
}

#endif
