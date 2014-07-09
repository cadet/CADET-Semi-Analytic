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

#ifndef CASEMA_GRMLAPLACESOLUTION_HPP_
#define CASEMA_GRMLAPLACESOLUTION_HPP_

#include "ModelData.hpp"

namespace casema
{
	namespace laplaceSolution
	{
		template <typename num_t, typename data_t, typename Inlet_t>
		class SingleComponentLinearRapidEquilibrium
		{
		public: 
			SingleComponentLinearRapidEquilibrium(const ModelData<data_t>& model, const Inlet_t& inlet) : _model(model), _inlet(inlet), _one("1.0"), _four("4.0")
			{ 
				precompute();
			}

			template <typename eval_t>
			eval_t operator()(const eval_t& s, const num_t& z) const
			{
				const eval_t sqrtAs = sqrt(_A_s_Factor / s);
				const eval_t f_s = _model.filmDiffusion[0] / (_model.filmDiffusion[0] - _B / _model.parRadius + _B / (sqrtAs * tanh(_model.parRadius / sqrtAs)));
				const eval_t phi_s = s + _phiFactor * (_one - f_s);
				const eval_t sqrtOnePlusFourPhiDivNuTimesU = sqrt(_one + _four * phi_s / _nuTimesU);
		        const eval_t lambda1 = _halfNu * (_one + sqrtOnePlusFourPhiDivNuTimesU);
		        const eval_t lambda2 = _halfNu * (_one - sqrtOnePlusFourPhiDivNuTimesU);

		        const eval_t commonAB = _inlet(s) / ( (lambda2 / _nu - _one) * lambda1 * exp(lambda1 * _model.colLength) + (_one - lambda1 / _nu) * lambda2 * exp(lambda2 * _model.colLength) );
		        const eval_t A = lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t B = lambda1 * exp(lambda1 * _model.colLength);
		        
		        return commonAB * (A * exp(lambda1 * z) - B * exp(lambda2 * z));
			}

			template <typename eval_t>
			eval_t operator()(const eval_t& s) const
			{
				return (*this)(s, _model.colLength);
			}

			template <typename eval_t>
			eval_t withoutInlet(const eval_t& s, const num_t& z) const
			{
				const eval_t sqrtAs = sqrt(_A_s_Factor / s);
				const eval_t f_s = _model.filmDiffusion[0] / (_model.filmDiffusion[0] - _B / _model.parRadius + _B / (sqrtAs * tanh(_model.parRadius / sqrtAs)));
				const eval_t phi_s = s + _phiFactor * (_one - f_s);
				const eval_t sqrtOnePlusFourPhiDivNuTimesU = sqrt(_one + _four * phi_s / _nuTimesU);
		        const eval_t lambda1 = _halfNu * (_one + sqrtOnePlusFourPhiDivNuTimesU);
		        const eval_t lambda2 = _halfNu * (_one - sqrtOnePlusFourPhiDivNuTimesU);

		        const eval_t commonAB = (lambda2 / _nu - _one) * lambda1 * exp(lambda1 * _model.colLength) + (_one - lambda1 / _nu) * lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t A = lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t B = lambda1 * exp(lambda1 * _model.colLength);
		        
		        return (A * exp(lambda1 * z) - B * exp(lambda2 * z)) / commonAB;
			}

			template <typename eval_t>
			eval_t withoutInlet(const eval_t& s) const
			{
				return withoutInlet(s, _model.colLength);
			}

			const char* const name() const { return "Single Component quasi-stationary Linear Binding"; }

		protected:
			const ModelData<data_t>& _model;
			const Inlet_t& _inlet;

			num_t _one;
			num_t _four;
			num_t _nu;
			num_t _nuTimesU;
			num_t _halfNu;
			num_t _phiFactor;
			num_t _B;
			num_t _A_s_Factor;

			void precompute()
			{
			    const num_t F = (_one - _model.colPorosity) / _model.colPorosity;			    
			    const num_t fPTimesKa = (_one - _model.parPorosity) / _model.parPorosity * _model.linearKA[0];
			    _nu = _model.velocity / _model.colDispersion;
			    _nuTimesU = _nu*_model.velocity;
			    _halfNu = _nu / num_t("2.0");
			    _phiFactor = _model.filmDiffusion[0] * num_t("3.0") / _model.parRadius * F;

			    _B = _model.parPorosity * _model.particleDiffusion[0] + (_one - _model.parPorosity) * _model.surfaceDiffusion[0] * _model.linearKA[0] / _model.linearKD[0];
			    _A_s_Factor = (_model.surfaceDiffusion[0] * fPTimesKa + _model.particleDiffusion[0] * _model.linearKD[0]) / (fPTimesKa + _model.linearKD[0]);
			}
		};

		template <typename num_t, typename data_t, typename Inlet_t>
		class SingleComponentLinearDynamic
		{
		public: 
			SingleComponentLinearDynamic(const ModelData<data_t>& model, const Inlet_t& inlet) : _model(model), _inlet(inlet), _one("1.0"), _four("4.0")
			{ 
				precompute();
			}

			template <typename eval_t>
			eval_t operator()(const eval_t& s, const num_t& z) const
			{
				const eval_t B_s = _bsFirst + _bsFactor / (_model.linearKD[0] + s);
				const eval_t sqrtAs = sqrt((_AsFirst + _model.particleDiffusion[0] * s) / (s * (_fPTimesKa + _model.linearKD[0] + s)));

				const eval_t f_s = _model.filmDiffusion[0] / (_model.filmDiffusion[0] - B_s / _model.parRadius + B_s / (sqrtAs * tanh(_model.parRadius / sqrtAs)));
				const eval_t phi_s = s + _phiFactor * (_one - f_s);
				const eval_t sqrtOnePlusFourPhiDivNuTimesU = sqrt(_one + _four * phi_s / _nuTimesU);
		        const eval_t lambda1 = _halfNu * (_one + sqrtOnePlusFourPhiDivNuTimesU);
		        const eval_t lambda2 = _halfNu * (_one - sqrtOnePlusFourPhiDivNuTimesU);

		        // Inlet for pulse injection: _model.constCoeff[0] * (_one - exp(-s * _model.sectionTimes[1])) / s
		        const eval_t commonAB = _inlet(s) / ( (lambda2 / _nu - _one) * lambda1 * exp(lambda1 * _model.colLength) + (_one - lambda1 / _nu) * lambda2 * exp(lambda2 * _model.colLength) );
		        const eval_t A = lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t B = lambda1 * exp(lambda1 * _model.colLength);
		        
		        return commonAB * (A * exp(lambda1 * z) - B * exp(lambda2 * z));
			}

			template <typename eval_t>
			eval_t operator()(const eval_t& s) const
			{
				return (*this)(s, _model.colLength);
			}

			template <typename eval_t>
			eval_t withoutInlet(const eval_t& s, const num_t& z) const
			{
				const eval_t B_s = _bsFirst + _bsFactor / (_model.linearKD[0] + s);
				const eval_t sqrtAs = sqrt((_AsFirst + _model.particleDiffusion[0] * s) / (s * (_fPTimesKa + _model.linearKD[0] + s)));

				const eval_t f_s = _model.filmDiffusion[0] / (_model.filmDiffusion[0] - B_s / _model.parRadius + B_s / (sqrtAs * tanh(_model.parRadius / sqrtAs)));
				const eval_t phi_s = s + _phiFactor * (_one - f_s);
				const eval_t sqrtOnePlusFourPhiDivNuTimesU = sqrt(_one + _four * phi_s / _nuTimesU);
		        const eval_t lambda1 = _halfNu * (_one + sqrtOnePlusFourPhiDivNuTimesU);
		        const eval_t lambda2 = _halfNu * (_one - sqrtOnePlusFourPhiDivNuTimesU);

		        // Inlet for pulse injection: _model.constCoeff[0] * (_one - exp(-s * _model.sectionTimes[1])) / s
		        const eval_t commonAB = (lambda2 / _nu - _one) * lambda1 * exp(lambda1 * _model.colLength) + (_one - lambda1 / _nu) * lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t A = lambda2 * exp(lambda2 * _model.colLength);
		        const eval_t B = lambda1 * exp(lambda1 * _model.colLength);
		        
		        return (A * exp(lambda1 * z) - B * exp(lambda2 * z)) / commonAB;
			}

			template <typename eval_t>
			eval_t withoutInlet(const eval_t& s) const
			{
				return withoutInlet(s, _model.colLength);
			}

			const char* const name() const { return "Single Component dynamic Linear Binding"; }

		protected:
			const ModelData<data_t>& _model;
			const Inlet_t& _inlet;

			num_t _one;
			num_t _four;
			num_t _nu;
			num_t _nuTimesU;
			num_t _halfNu;
			num_t _phiFactor;
			num_t _bsFirst;
			num_t _bsFactor;
			num_t _fPTimesKa;
			num_t _AsFirst;

			void precompute()
			{
			    const num_t F = (_one - _model.colPorosity) / _model.colPorosity;			    
			    _nu = _model.velocity / _model.colDispersion;
			    _nuTimesU = _nu*_model.velocity;
			    _halfNu = _nu / num_t("2.0");
			    _phiFactor = _model.filmDiffusion[0] * num_t("3.0") / _model.parRadius * F;

			    _bsFirst = _model.parPorosity * _model.particleDiffusion[0];
			    _bsFactor = (_one - _model.parPorosity) * _model.surfaceDiffusion[0] * _model.linearKA[0];
			    _fPTimesKa = (_one - _model.parPorosity) / _model.parPorosity * _model.linearKA[0];

			    _AsFirst = _model.surfaceDiffusion[0] * _fPTimesKa + _model.particleDiffusion[0] * _model.linearKD[0];
			}
		};

	}

}

#endif
