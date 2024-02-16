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

/**
 * @file 
 * ModelBuilder implementation
 */

#ifndef CASEMA_MODELBUILDER_IMPL_HPP_
#define CASEMA_MODELBUILDER_IMPL_HPP_

#include "casemaCompilerInfo.hpp"

#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

namespace casema
{

namespace io
{
	class IParameterProvider;
}

namespace model
{
	class UnitOperation;
	class ModelSystem;
}

/**
 * @brief Provides functionality to build a model
 */
class ModelBuilder
{
public:

	ModelBuilder();
	~ModelBuilder() CASEMA_NOEXCEPT;

	model::ModelSystem* createSystem(io::IParameterProvider& paramProvider);
	model::ModelSystem* createSystem();
	void detachSystem(model::ModelSystem const* sys);
	void destroySystem(model::ModelSystem* sys);

	model::UnitOperation* createUnitOperation(io::IParameterProvider& paramProvider, int uoId);
	model::UnitOperation* createUnitOperation(const std::string& uoType, int uoId);
	void destroyUnitOperation(model::UnitOperation* unitOp);

protected:

	/**
	 * @brief Registers an IUnitOperation
	 * @param [in] name Name of the model
	 * @tparam UnitOpModel_t Type of the model
	 */
	template <class UnitOpModel_t>
	void registerModel(const std::string& name);

	/**
	 * @brief Registers an IUnitOperation
	 * @details The name of the model is inferred from the static function IUnitOperation::identifier().
	 * @tparam UnitOpModel_t Type of the model
	 */
	template <class UnitOpModel_t>
	void registerModel();

	typedef std::unordered_map<std::string, std::function<model::UnitOperation*(int)>> ModelFactoryContainer_t;

	ModelFactoryContainer_t _modelCreators; //!< Map with factory functions for models

	std::vector<model::ModelSystem*> _models; //!< Models
};

} // namespace casema

#endif  // CASEMA_MODELBUILDER_IMPL_HPP_
