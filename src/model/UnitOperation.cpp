// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//
//  Copyright © 2015-present: Samuel Leweke¹² and the CADET-Semi-Analytic
//  authors, see the AUTHORS file
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

#include <iomanip>
#include <sstream>
#include <iostream>

namespace casema
{

namespace model
{

std::string getUnitName(io::IParameterProvider& paramProvider)
{
	std::string uoType = paramProvider.getString("UNIT_TYPE");

	if (uoType == "COLUMN_MODEL_1D" || uoType == "COLUMN_MODEL_2D")
	{
		const int nParType = paramProvider.getInt("NPARTYPE");

		if (nParType < 0)
			throw InvalidParameterException("Field NPARTYPE must be >= 0");

		if (nParType == 0)
		{
			if (uoType == "COLUMN_MODEL_1D")
				uoType = "LUMPED_RATE_MODEL_WITHOUT_PORES";
			else
				throw InvalidParameterException("CASEMA does not support a 2D bulk transport model without particles. Mimic this by setting film diffusion coefficients to zero");
		}
		else
		{
			for(int type = 0; type < nParType; type++)
			{
				std::ostringstream oss;
				oss << "particle_type_"  << std::setw(3) << std::setfill('0') << type;

				paramProvider.pushScope(oss.str());

				bool hasFilmDiff = paramProvider.getBool("HAS_FILM_DIFFUSION");
				bool hasPoreDiff = false;
				bool hasSurfDiff = false;
				
				if(hasFilmDiff)
				{
					hasPoreDiff = paramProvider.getBool("HAS_PORE_DIFFUSION");
					hasSurfDiff = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;
				}
				else if (nParType > 1)
				{
					throw InvalidParameterException("Multiple particle types are not supported for transport model without pores (HAS_FILM_DIFFUSION=False)");
				}

				if (uoType == "COLUMN_MODEL_1D" || uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES" || uoType == "LUMPED_RATE_MODEL_WITH_PORES" || uoType == "GENERAL_RATE_MODEL")
				{
					if(!hasFilmDiff)
					{
						if (type > 0 && uoType != "LUMPED_RATE_MODEL_WITHOUT_PORES")
							throw InvalidParameterException("Only particles of same type are supported for NPARTYPE > 1. Double check fields HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION.");

						uoType = "LUMPED_RATE_MODEL_WITHOUT_PORES";
					}
					else
					{
						if (!hasPoreDiff && !hasSurfDiff)
						{
							if (type > 0 && uoType != "LUMPED_RATE_MODEL_WITH_PORES")
								throw InvalidParameterException("Only particles of same type are supported for NPARTYPE > 1. Double check fields HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION.");

							uoType = "LUMPED_RATE_MODEL_WITH_PORES";
						}
						else
						{
						if (type > 0 && uoType != "GENERAL_RATE_MODEL")
							throw InvalidParameterException("Only particles of same type are supported for NPARTYPE > 1. Double check fields HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION.");

							uoType = "GENERAL_RATE_MODEL";
						}
					}
				}
				else if (uoType == "COLUMN_MODEL_2D" || uoType == "GENERAL_RATE_MODEL_2D")
				{
					if(!hasFilmDiff)
						throw InvalidParameterException("CASEMA does not support a 2D bulk transport model without pores (2DLRM). Mimic this by setting film diffusion coefficients to a large value");
					else if (!hasPoreDiff && !hasSurfDiff)
						throw InvalidParameterException("CASEMA does not support a 2D bulk transport model with particles without pore and surface diffusion (2DLRMP). Only 2D GRM is supported");
					else
						uoType = "GENERAL_RATE_MODEL_2D";
				}
				else
					return uoType;

				paramProvider.popScope(); // particle type group
			}
		}
	}

	return uoType;
}

bool UnitOperation::configure(io::IParameterProvider& paramProvider)
{
	_nComp = paramProvider.getInt("NCOMP");

	if (_nComp != 1)
		throw InvalidParameterException("Field NCOMP has to be 1");

	return true;
}

mpfr::mpreal UnitOperation::initialOutletValue(int port) const CASEMA_NOEXCEPT
{
	(void)port;
	return mpfr::mpreal(0.0);
}

}  // namespace model

}  // namespace casema
