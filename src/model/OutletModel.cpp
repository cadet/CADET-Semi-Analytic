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

#include "model/OutletModel.hpp"

#include <unordered_map>
#include <string>

namespace casema
{

namespace model
{

bool OutletModel::configure(io::IParameterProvider& paramProvider)
{
	UnitOperation::configure(paramProvider);
	return true;
}

void registerOutletModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[OutletModel::identifier()] = [](int uoId) { return new OutletModel(uoId); };
}

}  // namespace model

}  // namespace casema
