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

#include "model/UnitOperation.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"

namespace casema
{

namespace model
{

bool UnitOperation::configure(io::IParameterProvider& paramProvider)
{
	_nComp = paramProvider.getInt("NCOMP");

	if (_nComp != 1)
		throw InvalidParameterException("Field NCOMP has to be 1");

	return true;
}

}  // namespace model

}  // namespace casema
