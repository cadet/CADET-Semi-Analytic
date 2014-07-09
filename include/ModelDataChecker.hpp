// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_MODELDATA_CHECKER_HPP_
#define CASEMA_MODELDATA_CHECKER_HPP_

#include <iostream>

#include "ModelData.hpp"

namespace casema
{

	template <typename real_t>
	bool checkModelForCompatibility(const ModelData<real_t>& model)
	{
	    if (model.nComponents != 1)
	    {
	        std::cout << "ERROR: Only single component models allowed" << std::endl;
	        return false;
	    }

	    if ((model.initialLiquidConcentration[0] != real_t(0)) || (model.initialSolidConcentration[0] != real_t(0)))
	    {
	        std::cout << "ERROR: Initial concentrations must be 0.0" << std::endl;
	        return false;
	    }

	    if (model.bindingModel != casema::LINEAR)
	    {
	        std::cout << "ERROR: Only linear binding model allowed" << std::endl;
	        return false;
	    }
	    return true;
	}

}

#endif
