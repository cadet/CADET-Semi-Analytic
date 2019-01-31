// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2019: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_ANALYTICALMOMENTS_HPP_
#define CASEMA_ANALYTICALMOMENTS_HPP_

#include <vector>

#include "ModelData.hpp"

/**
 * The analytical moments are taken from S. Qamar et al. (2014). 
 * "Analytical solutions and moment analysis of general rate model for linear liquid chromatography".
 * Chemical Engineering Science, 107, 192–205, doi: [10.1016/j.ces.2013.12.019](http://dx.doi.org/10.1016/j.ces.2013.12.019)) 
 **/

namespace casema
{

	namespace detail
	{
		template <typename T>
		T injectedMass(const ModelData<T>& model)
		{
			T mass = T(0);

			for (std::size_t i = 0; i < model.constCoeff.size(); ++i)
			{
				const T dt = model.sectionTimes[i+1] - model.sectionTimes[i];
				mass += dt * (model.constCoeff[i] + dt * (model.linCoeff[i] / T(2) + dt * (model.quadCoeff[i] / T(3) + model.cubicCoeff[i] * dt / T(4))));
			}

			return mass;
		}
	}

	template <typename T>
	std::vector<T> analyticalMoments(const ModelData<T>& model)
	{
		std::vector<T> mom(5);

		const T one = T(1);
		const T a_star = (model.parPorosity + (one - model.parPorosity) * model.linearKA[0] / model.linearKD[0]);
		const T F = (one - model.colPorosity) / model.colPorosity;
		const T LuOverD_L = model.velocity / model.colDispersion * model.colLength;
		const T Lu = model.velocity * model.colLength;
		const T onePlusAF = (one + a_star * F);
		const T D_eff = model.parPorosity * model.particleDiffusion[0] + (one - model.parPorosity) * model.surfaceDiffusion[0] * model.linearKA[0] / model.linearKD[0];
		const T k_f = model.filmDiffusion[0];

		mom[0] = detail::injectedMass(model);
		mom[1] = model.sectionTimes[1] / T(2) + model.colLength / model.velocity * onePlusAF;
		
		mom[2] = sqr(model.sectionTimes[1]) / T(3) + model.sectionTimes[1] * model.colLength / model.velocity * onePlusAF;
		mom[2] += T(2) * sqr(model.colDispersion / model.velocity) * sqr(onePlusAF / model.velocity) * (LuOverD_L * (one + LuOverD_L / T(2)) + expm1(-LuOverD_L));
		mom[2] += F * model.colLength / model.velocity * T(2) / T(3) * model.parRadius * sqr(a_star) * (one / k_f + model.parRadius / (T(5) * D_eff));

		mom[3] = sqr(onePlusAF / sqr(model.velocity)) * onePlusAF / sqr(model.velocity) * (sqr(Lu) * Lu + T(6) * model.colDispersion * sqr(Lu) + sqr(model.colDispersion) * Lu * (T(24) + T(18) * expm1(-LuOverD_L)) + T(24) * sqr(model.colDispersion) * model.colDispersion * expm1(-LuOverD_L));
		mom[3] += T(4) * model.colLength * model.colDispersion * F / model.velocity * sqr(a_star / model.velocity) * model.parRadius * onePlusAF * (one / k_f + model.parRadius / (T(5) * D_eff)) * (expm1(-LuOverD_L) / LuOverD_L + one + LuOverD_L / T(2));
		mom[3] += model.sectionTimes[1] * model.colLength / model.velocity * F * sqr(a_star) * model.parRadius * (one / k_f + model.parRadius / (T(5) * D_eff));
		mom[3] += model.colLength / model.velocity * F * sqr(a_star * model.parRadius) * a_star * T(2) / T(3) * ( T(2) / T(35) * sqr(model.parRadius / D_eff) + T(2) / T(5) * model.parRadius / (D_eff * k_f) + one / sqr(k_f) );
		mom[3] += sqr(model.sectionTimes[1]) * (model.sectionTimes[1] / T(4) + model.colLength / model.velocity * onePlusAF);
		mom[3] += T(3) * model.sectionTimes[1] * sqr((one + model.colPorosity * F) / model.velocity) * model.colLength / model.velocity * model.colDispersion * (expm1(-LuOverD_L) / LuOverD_L + one + LuOverD_L / T(2));

		mom[4] = pow(model.colDispersion * onePlusAF / sqr(model.velocity), T(4)) * (T(24) * expm1(-T(2)*LuOverD_L) + (T(312) + T(108) * sqr(LuOverD_L)) * expm1(-LuOverD_L) + T(360) * LuOverD_L * exp(-LuOverD_L) + sqr(LuOverD_L) * (T(156) + LuOverD_L * (T(12) + LuOverD_L * one )));
		mom[4] += T(4) * sqr(onePlusAF * model.colDispersion * a_star / model.velocity) * model.colDispersion / pow(model.velocity, T(4)) * F * model.parRadius * (one / k_f + model.parRadius / (T(5) * D_eff)) * ( expm1(-LuOverD_L) * (T(24) + T(18) * LuOverD_L) + LuOverD_L * (T(24) + LuOverD_L * (T(6) + LuOverD_L)));
		mom[4] += T(8) / T(3) * model.colLength * model.colDispersion * F * pow(a_star / model.velocity, T(3)) * (T(2) + T(3) * model.colPorosity * F) * sqr(model.parRadius) * (one / sqr(k_f) + T(2) / T(5) * model.parRadius / (k_f * D_eff) + sqr(model.parRadius / D_eff) * T(2) / T(35)) * (one + expm1(-LuOverD_L) / LuOverD_L + LuOverD_L);
		mom[4] += T(8) * model.colLength * F / model.velocity * sqr(a_star * model.parRadius) * model.parRadius * sqr(a_star) * (one / (T(9) * pow(k_f,T(3))) + model.parRadius / (T(15) * D_eff * sqr(k_f)) + T(3) / T(175) * sqr(model.parRadius / D_eff) / k_f + pow(model.parRadius / D_eff, T(3)) / T(175));
		mom[4] += -T(4) / T(175) * model.colLength * model.colDispersion * sqr(F / model.velocity * a_star * model.parRadius / D_eff) * sqr(a_star * model.parRadius) / model.velocity * (T(2) + LuOverD_L + T(2) * expm1(-LuOverD_L) / LuOverD_L) + sqr(model.sectionTimes[1]) * model.sectionTimes[1] * (model.colLength / model.velocity * onePlusAF + model.sectionTimes[1] / T(5));
		mom[4] += T(2) * sqr(model.sectionTimes[1]) * model.colLength * (F / model.velocity * T(2) / T(3) * model.parRadius * sqr(a_star) * (one / k_f + model.parRadius / (T(5) * D_eff)) + model.colDispersion / model.velocity * sqr(onePlusAF / model.velocity) * (T(2) + LuOverD_L + T(2) * expm1(-LuOverD_L) / LuOverD_L) );
		mom[4] += model.sectionTimes[1] * T(2) * model.colLength * sqr(model.colDispersion / model.velocity) * pow((one + model.colPorosity * F) / model.velocity, T(3)) * (T(24) * expm1(-LuOverD_L) / LuOverD_L + T(24) + T(18) * expm1(-LuOverD_L) + LuOverD_L * (T(6) + LuOverD_L));
		mom[4] += model.sectionTimes[1] * T(4) * model.colLength * model.colDispersion / model.velocity * (one + model.colPorosity * F) * F * sqr(a_star / model.velocity) * model.parRadius * (one / k_f + model.parRadius / (T(5) * D_eff)) * (T(2) + LuOverD_L + T(2) * expm1(-LuOverD_L) / LuOverD_L);
		mom[4] += model.sectionTimes[1] * T(4) / T(3) * sqr(model.parRadius*a_star) * a_star * model.colLength * F / model.velocity * (one / sqr(k_f) + T(2) / T(5) * model.parRadius / (k_f * D_eff) + T(2) / T(35) * sqr(model.parRadius / D_eff));
		mom[4] += pow(model.sectionTimes[1], T(3)) * (model.colLength / model.velocity * onePlusAF + model.sectionTimes[1] / T(5));

		return mom;
	}

	template <typename T>
	std::vector<T> analyticalCentralMoments(const ModelData<T>& model)
	{
		std::vector<T> mom(5);

		const T one = T(1);
		const T a_star = (model.parPorosity + (one - model.parPorosity) * model.linearKA[0] / model.linearKD[0]);
		const T F = (one - model.colPorosity) / model.colPorosity;
		const T D_eff = model.parPorosity * model.particleDiffusion[0] + (one - model.parPorosity) * model.surfaceDiffusion[0] * model.linearKA[0] / model.linearKD[0];
		const T onePlusAF = (one + a_star * F);
		const T LuOverD_L = model.velocity / model.colDispersion * model.colLength;
		const T Lu = model.velocity * model.colLength;
		const T k_f = model.filmDiffusion[0];

		mom[0] = detail::injectedMass(model);
		mom[1] = model.sectionTimes[1] / T(2) + model.colLength / model.velocity * onePlusAF;
		
		mom[2] = sqr(model.sectionTimes[1]) / T(12) + T(2) * model.colLength / model.velocity * model.colDispersion * sqr((one + a_star * F) / model.velocity) * (one + expm1(-LuOverD_L) / LuOverD_L );
		mom[2] += model.colLength * F / model.velocity * T(2) / T(3) * model.parRadius * sqr(a_star) * (one / k_f + model.parRadius / (T(5) * D_eff));

		mom[3] = T(12) * model.colLength * sqr(model.colDispersion * onePlusAF / sqr(model.velocity)) * onePlusAF / model.velocity * ( (one + T(2) / LuOverD_L) * expm1(-LuOverD_L) + T(2));
		mom[3] += T(4) * model.colLength / model.velocity * model.colDispersion * F * sqr(a_star / model.velocity) * onePlusAF * model.parRadius * (one / k_f + model.parRadius / (T(5) * D_eff)) * (expm1(-LuOverD_L) / LuOverD_L + one);
		mom[3] += model.colLength / model.velocity * F * sqr(a_star * model.parRadius) * a_star * T(2) / T(3) * (T(2) / T(35) * sqr(model.parRadius / D_eff) + T(2) / T(5) * model.parRadius / (D_eff * k_f) + one / sqr(k_f));

		mom[4] = T(12) * pow(onePlusAF * model.colDispersion / sqr(model.velocity), T(4)) * (T(2) * expm1(-T(2) * LuOverD_L) + expm1(-LuOverD_L) * (T(26) + LuOverD_L * (T(22) + LuOverD_L * T(4) )) + LuOverD_L * (T(30) + LuOverD_L * T(5)));
		mom[4] += T(8) * sqr(onePlusAF * model.colDispersion / sqr(model.velocity) * a_star) * model.colDispersion * F / sqr(model.velocity) * model.parRadius * (one / k_f + model.parRadius / (T(5) * D_eff)) * (expm1(-LuOverD_L) * (T(12) + T(7) * LuOverD_L) + LuOverD_L * (T(12) + LuOverD_L));
		mom[4] += T(4) / T(3) * model.colLength * model.colDispersion * F * pow(a_star / model.velocity, T(3)) * sqr(model.parRadius) * (one / sqr(k_f) + T(2) * model.parRadius / (T(5) * k_f * D_eff) + T(2) / T(35) * sqr(model.parRadius / D_eff)) * ((T(4) + T(6) * a_star * F) * (one + expm1(-LuOverD_L) / LuOverD_L) + a_star * F * LuOverD_L);
		mom[4] += T(8) * model.colLength * F * a_star * pow(model.parRadius * a_star, T(3)) / model.velocity * (one / (T(9) * pow(k_f, T(3))) + model.parRadius / (T(15) * D_eff * sqr(k_f)) + T(3) / T(175) * sqr(model.parRadius / D_eff) / k_f + pow(model.parRadius / D_eff, T(3)) / T(525));
		mom[4] += -T(4) / T(175) * model.colLength * model.colDispersion * sqr(F / D_eff * sqr(model.parRadius * a_star) / model.velocity) / model.velocity * (LuOverD_L + T(2) * (one + expm1(-LuOverD_L) / LuOverD_L)) + pow(model.sectionTimes[1], T(4)) / T(80);
		mom[4] += sqr(model.sectionTimes[1]) * (sqr(model.colDispersion * (one + model.colPorosity * F) / sqr(model.velocity)) * (LuOverD_L + expm1(-LuOverD_L)) + model.colLength * F / model.velocity * model.parRadius * sqr(a_star) / T(3) * (one / k_f + model.parRadius / (T(5) * D_eff)) );

		return mom;
	}

}

#endif
