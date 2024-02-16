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

#ifndef CASEMA_BESSELZEROS_HPP_
#define CASEMA_BESSELZEROS_HPP_

#include "casemaCompilerInfo.hpp"

namespace mpfr
{
	class mpreal;
}

namespace casema
{
	void besselZerosJ1(int n, mpfr::mpreal* out) CASEMA_NOEXCEPT;
}

#endif
