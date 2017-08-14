// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2017: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_MODELDATA_HPP_
#define CASEMA_MODELDATA_HPP_

#include <vector>

#include "CadetEnumeration.hpp"

namespace casema
{

	enum ChromatographyUnitType : int
	{
		GeneralRateModel,
		LumpedRateModelWithPores,
		LumpedRateModelWithoutPores
	};

	template <typename real_t>
	struct ModelData
	{
		ChromatographyUnitType modelType;
		std::size_t nComponents;
		
		std::vector<real_t> initialLiquidConcentration;
		std::vector<real_t> initialSolidConcentration;

		real_t colDispersion;
		real_t colLength;
		real_t colPorosity;

		std::vector<real_t> filmDiffusion;
		std::vector<real_t> particleDiffusion;
		std::vector<real_t> surfaceDiffusion;

		real_t parRadius;
		real_t parPorosity;

		real_t velocity;

		bool kineticBinding;

		casema::AdsorptionType bindingModel;
		std::vector<real_t> linearKA;
		std::vector<real_t> linearKD;

		std::size_t nInletSections;
		std::vector<real_t> sectionTimes;
		std::vector<real_t> constCoeff;
		std::vector<real_t> linCoeff;
		std::vector<real_t> quadCoeff;
		std::vector<real_t> cubicCoeff;

		bool writeUserTimes;
		std::vector<real_t> outletTimes;

		inline real_t totalPorosity() const { return colPorosity + (real_t(1) - colPorosity) * parPorosity; }
	};

	template <typename real_t>
	real_t maxSimulationTime(const ModelData<real_t>& model)
	{
		if (!model.outletTimes.empty())
		{
			return max(model.outletTimes[model.outletTimes.size() - 1], model.sectionTimes[model.sectionTimes.size() - 1]);
		}
		return model.sectionTimes[model.sectionTimes.size() - 1];
	}

	template <typename real_t>
	bool inletIsStep(const ModelData<real_t>& model)
	{
		// Always true for piecewise polynomials
		return true;
	}

}

#endif
