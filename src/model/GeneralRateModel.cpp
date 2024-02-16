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

#include "model/GeneralRateModel.hpp"
#include "io/ParameterProvider.hpp"

#include <unordered_map>

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
GeneralRateModel::GeneralRateModel(int unitOpIdx) : ColumnWithPoreDiffusion(unitOpIdx)
{
}

GeneralRateModel::~GeneralRateModel() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool GeneralRateModel::configure(io::IParameterProvider& paramProvider)
{
	return ColumnWithPoreDiffusion::configure(paramProvider);
}

bool GeneralRateModel::hasValidEstimate() const CASEMA_NOEXCEPT
{
	if (_initC != 0.0)
		return false;

	bool v = true;
	for (int j = 0; j < _parTypeVolFrac.size(); ++j)
	{
		if ((_parSurfDiffusion[j] != 0.0) && (_kA[j] != 0.0) && _isKinetic[j])
		{
			v = false;
			break;
		}
	}

	return v;
}


void registerGeneralRateModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[GeneralRateModel::identifier()] = [](int uoId) { return new GeneralRateModel(uoId); };
	models["GRM"] = [](int uoId) { return new GeneralRateModel(uoId); };
}

}  // namespace model

}  // namespace casema
