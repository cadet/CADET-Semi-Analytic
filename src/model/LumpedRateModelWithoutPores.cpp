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

#include "model/LumpedRateModelWithoutPores.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"

#include <unordered_map>

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
LumpedRateModelWithoutPores::LumpedRateModelWithoutPores(int unitOpIdx) : ColumnLikeModel(unitOpIdx)
{
}

LumpedRateModelWithoutPores::~LumpedRateModelWithoutPores() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool LumpedRateModelWithoutPores::configure(io::IParameterProvider& paramProvider)
{
	ColumnLikeModel::configure(paramProvider);

	if (_isKinetic[0])
		_factor = (_one - _colPorosity) / _colPorosity * _kA[0];
	else
		_factor = (_one + (_one - _colPorosity) / _colPorosity * _kA[0] / _kD[0]);

	return true;
}

template <typename eval_t>
eval_t LumpedRateModelWithoutPores::phiImpl(const eval_t& s) const CASEMA_NOEXCEPT
{
	if (_isKinetic[0])
		return s * (_one + _factor / (_kD[0] + s));
	else
		return s * _factor;
}


void registerLumpedRateModelWithoutPores(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	const std::function<UnitOperation*(int)> f = [](int uoId) { return new LumpedRateModelWithoutPores(uoId); };
	models[LumpedRateModelWithoutPores::identifier()] = f;
	models["LRM"] = f;
	models["DPFR"] = f;
}

}  // namespace model

}  // namespace casema
