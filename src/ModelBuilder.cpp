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

#include "ModelBuilder.hpp"
#include "io/ParameterProvider.hpp"
#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"
#include "Exceptions.hpp"

#include <sstream>
#include <iomanip>

namespace casema
{
	namespace model
	{
		void registerInletModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);
		void registerOutletModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);

		void registerGeneralRateModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);
		void registerLumpedRateModelWithPores(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);
		void registerLumpedRateModelWithoutPores(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);
		void registerCSTRModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);
		void registerGeneralRateModel2D(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models);

	} // namespace model

	ModelBuilder::ModelBuilder()
	{
		// Register all available models
		model::registerInletModel(_modelCreators);
		model::registerOutletModel(_modelCreators);
		model::registerGeneralRateModel(_modelCreators);
		model::registerLumpedRateModelWithPores(_modelCreators);
		model::registerLumpedRateModelWithoutPores(_modelCreators);
		model::registerCSTRModel(_modelCreators);
		model::registerGeneralRateModel2D(_modelCreators);
	}

	ModelBuilder::~ModelBuilder() CASEMA_NOEXCEPT
	{
		for (model::ModelSystem* model : _models)
			delete model;
	}

	template <class UnitOpModel_t>
	void ModelBuilder::registerModel(const std::string& name)
	{
		_modelCreators[name] = [](int uoId) { return new UnitOpModel_t(uoId); };
	}

	template <class UnitOpModel_t>
	void ModelBuilder::registerModel()
	{
		registerModel<UnitOpModel_t>(UnitOpModel_t::identifier());
	}

	model::ModelSystem* ModelBuilder::createSystem(io::IParameterProvider& paramProvider)
	{
		model::ModelSystem* sys = new model::ModelSystem();
		
		// Create and configure all unit operations
		bool success = true;
		unsigned int i = 0;
		std::ostringstream oss;
		oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		while (paramProvider.exists(oss.str()))
		{
			// Create and configure unit operation
			paramProvider.pushScope(oss.str());

			model::UnitOperation* const unitOp = createUnitOperation(paramProvider, i);

			paramProvider.popScope();

			if (unitOp)
			{
				// Model correctly created and configured -> add to system
				sys->addModel(unitOp);
			}
			else
			{
				// Something went wrong -> abort and exit
				success = false;
				break;
			}

			++i;
			oss.str("");
			oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		}

		// Configure the whole system
		success = success && sys->configure(paramProvider);

		if (success)
		{
			_models.push_back(sys);
			return sys;
		}
		else
		{
			delete sys;
			return nullptr;
		}
	}

	model::ModelSystem* ModelBuilder::createSystem()
	{
		model::ModelSystem* sys = new model::ModelSystem();
		_models.push_back(sys);
		return sys;
	}

	void ModelBuilder::detachSystem(model::ModelSystem const* sys)
	{
		for (std::vector<model::ModelSystem*>::iterator it = _models.begin(); it != _models.end(); ++it)
		{
			if (*it == sys)
			{
				_models.erase(it);
				break;
			}
		}
	}

	void ModelBuilder::destroySystem(model::ModelSystem* sys)
	{
		delete sys;
	}

	model::UnitOperation* ModelBuilder::createUnitOperation(io::IParameterProvider& paramProvider, int uoId)
	{
		const std::string uoType = paramProvider.getString("UNIT_TYPE");
		const auto it = _modelCreators.find(uoType);
		if (it == _modelCreators.end())
		{
			// Model was not found
			return nullptr;
		}

		// Call factory function (thanks to type erasure of std::function we can store 
		// all factory functions in one container)
		model::UnitOperation* const model = it->second(uoId);

		if (!model->configure(paramProvider))
		{
			delete model;
			return nullptr;
		}

		return model;
	}

	model::UnitOperation* ModelBuilder::createUnitOperation(const std::string& uoType, int uoId)
	{
		const auto it = _modelCreators.find(uoType);
		if (it == _modelCreators.end())
		{
			// Model was not found
			return nullptr;
		}

		model::UnitOperation* const model = it->second(uoId);
		return model;
	}

	void ModelBuilder::destroyUnitOperation(model::UnitOperation* unitOp)
	{
		delete unitOp;
	}

} // namespace casema
