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

#ifndef CASEMA_INLETLAPLACE_HPP_
#define CASEMA_INLETLAPLACE_HPP_

#include "ModelData.hpp"

namespace casema
{
	namespace laplaceSolution
	{
		template <typename num_t, typename data_t>
		class Inlet
		{
		public: 
			Inlet(const ModelData<data_t>& model) : _model(model), _two("2.0"), _three("3.0"), _six("6.0") { } 

			template <typename eval_t>
			eval_t operator()(const eval_t& s) const
			{
				eval_t val("0.0");
				for (std::size_t i = 0; i < _model.nInletSections; ++i)
				{
					val += evalSection(i, s);
				}
				return val;
			}

			template <typename eval_t>
			eval_t timeDomain(const eval_t& s) const
			{
				for (std::size_t i = 0; i < _model.nInletSections; ++i)
				{
					if ((_model.sectionTimes[i] <= s) && (s < _model.sectionTimes[i+1]))
						return evalTimeDomainSection(i, s);
				}

				if (s == _model.sectionTimes[_model.nInletSections])
					return evalTimeDomainSection(_model.nInletSections - 1, s);

				return 0;
			}

		protected:
			const ModelData<data_t>& _model;
			const num_t _two;
			const num_t _three;
			const num_t _six;

			template <typename eval_t>
			eval_t evalSection(std::size_t i, const eval_t& s) const
			{
				const num_t Q = _model.sectionTimes[i];
				const num_t T = _model.sectionTimes[i+1];
				const num_t dt = T - Q;
				const num_t a = _model.constCoeff[i];
				const num_t b = _model.linCoeff[i];
				const num_t c = _model.quadCoeff[i];
				const num_t d = _model.cubicCoeff[i];

				const eval_t ff = (a + (b + (_two*c + _six * d / s) / s) / s) / s;
				const eval_t tf = dt * ( (b + (_two * c + _six * d / s) / s) / s + dt * (((_three / s + dt) * d  + c) / s)) + (a + (b + (_two*c + _six * d / s) / s) / s) / s;
				return ff * exp(-s*Q) - tf * exp(-s*T);
			}

			template <typename eval_t>
			eval_t evalTimeDomainSection(std::size_t i, const eval_t& s) const
			{
				const eval_t T = s - _model.sectionTimes[i];
				
				const num_t a = _model.constCoeff[i];
				const num_t b = _model.linCoeff[i];
				const num_t c = _model.quadCoeff[i];
				const num_t d = _model.cubicCoeff[i];
				
				return a + T * (b + T * (c + T * d));
			}			
		};

	}

}

#endif
