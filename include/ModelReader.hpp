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

#ifndef CASEMA_MODELREADER_HPP_
#define CASEMA_MODELREADER_HPP_

#include <iomanip>
#include <sstream>
#include <string>

#include "ModelData.hpp"

namespace casema
{

	template <class reader_t, typename real_t>
	ModelData<real_t> readModel(reader_t& reader)
	{
		ModelData<real_t> data;

		reader.setGroup(e2s(GRP_IN_MODEL));
		s2e(reader.template scalar<std::string>(e2s(ADSORPTION_TYPE)), data.bindingModel);
		data.nComponents = reader.template scalar<int>(e2s(NCOMP));

		reader.setGroup(e2s(GRP_IN_INLET));
		data.nInletSections  = reader.template scalar<int>(e2s(NSEC));

		// Reserve space
		data.initialLiquidConcentration.reserve(data.nComponents);
		data.initialSolidConcentration.reserve(data.nComponents);
		data.filmDiffusion.reserve(data.nComponents);
		data.particleDiffusion.reserve(data.nComponents);
		data.surfaceDiffusion.reserve(data.nComponents);
		data.linearKA.reserve(data.nComponents);
		data.linearKD.reserve(data.nComponents);

		data.sectionTimes.reserve(data.nInletSections);
		data.constCoeff.reserve(data.nInletSections * data.nComponents);
		data.linCoeff.reserve(data.nInletSections * data.nComponents);
		data.quadCoeff.reserve(data.nInletSections * data.nComponents);
		data.cubicCoeff.reserve(data.nInletSections * data.nComponents);

		// Read inlet specs
		reader.setGroup(e2s(GRP_IN_INLET));
		const std::vector<double> sectionTimes = reader.template vector<double>(e2s(SECTION_TIMES));
		for(std::vector<double>::const_iterator it = sectionTimes.begin(); it != sectionTimes.end(); ++it)
		{
			data.sectionTimes.push_back(real_t(*it));
		}

		for (std::size_t sec = 0; sec < data.nInletSections; ++sec)
		{
			std::ostringstream oss;
			oss.str("");
			oss << e2s(GRP_IN_INLET) << "/sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << sec;
			reader.setGroup(oss.str());
			const std::vector<double> const_coeff = reader.template vector<double>("CONST_COEFF");
			const std::vector<double> linear_coeff = reader.template vector<double>("LIN_COEFF");
			const std::vector<double> quad_coeff = reader.template vector<double>("QUAD_COEFF");
			const std::vector<double> cubic_coeff = reader.template vector<double>("CUBE_COEFF");
			for(std::size_t i = 0; i < const_coeff.size(); ++i)
			{
				data.constCoeff.push_back(real_t(const_coeff[i]));
				data.linCoeff.push_back(real_t(linear_coeff[i]));
				data.quadCoeff.push_back(real_t(quad_coeff[i]));
				data.cubicCoeff.push_back(real_t(cubic_coeff[i]));
			}
		}

		// Read binding model specs
		reader.setGroup(e2s(GRP_IN_ADSORPTION));
		data.kineticBinding = reader.template scalar<int>(e2s(IS_KINETIC));

		if (data.bindingModel == LINEAR)
		{
			const std::vector<double> linKA = reader.template vector<double>(e2s(LIN_KA));
			const std::vector<double> linKD = reader.template vector<double>(e2s(LIN_KD));

			for(std::size_t i = 0; i < data.nComponents; ++i)
			{
				data.linearKA.push_back(real_t(linKA[i]));
				data.linearKD.push_back(real_t(linKD[i]));
			}
		}

		// Chromatography model parameters
		reader.setGroup(e2s(GRP_IN_MODEL));

		data.colDispersion = real_t(reader.template scalar<double>(e2s(COL_DISPERSION)));
		data.velocity = real_t(reader.template scalar<double>(e2s(VELOCITY)));
		data.colLength = real_t(reader.template scalar<double>(e2s(COL_LENGTH)));
		data.colPorosity = real_t(reader.template scalar<double>(e2s(COL_POROSITY)));
		data.parRadius = real_t(reader.template scalar<double>(e2s(PAR_RADIUS)));
		data.parPorosity = real_t(reader.template scalar<double>(e2s(PAR_POROSITY)));

		// Vectorial parameters
		const std::vector<double> initialSolid = reader.template vector<double>(e2s(INIT_Q));
		const std::vector<double> initialLiquid = reader.template vector<double>(e2s(INIT_C));
		const std::vector<double> filmDiff = reader.template vector<double>(e2s(FILM_DIFFUSION));
		const std::vector<double> parDiff = reader.template vector<double>(e2s(PAR_DIFFUSION));
		const std::vector<double> parSurfDiff = reader.template vector<double>(e2s(PAR_SURFDIFFUSION));
		for(std::size_t i = 0; i < data.nComponents; ++i)
		{
			data.initialSolidConcentration.push_back(real_t(initialSolid[i]));
			data.initialLiquidConcentration.push_back(real_t(initialLiquid[i]));
			data.filmDiffusion.push_back(real_t(filmDiff[i]));
			data.particleDiffusion.push_back(real_t(parDiff[i]));
			data.surfaceDiffusion.push_back(real_t(parSurfDiff[i]));
		}

		// Read outlet times
		reader.setGroup(e2s(GRP_IN_SOLVER));
		data.writeUserTimes = reader.template scalar<int>(e2s(WRITE_AT_USER_TIMES));

		if (data.writeUserTimes)
		{
			const std::vector<double> out = reader.template vector<double>(e2s(USER_SOLUTION_TIMES));
			data.outletTimes.reserve(out.size());
			for(std::vector<double>::const_iterator it = out.begin(); it != out.end(); ++it)
			{
				data.outletTimes.push_back(real_t(*it));
			}

		}

		return data;			
	}
}

#endif
