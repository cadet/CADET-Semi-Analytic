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
 * Defines the lumped rate model without pores (LRM).
 */

#ifndef CASEMA_LUMPEDRATEMODELWITHOUTPORES_HPP_
#define CASEMA_LUMPEDRATEMODELWITHOUTPORES_HPP_

#include "model/ColumnLikeModel.hpp"

#include <vector>

namespace casema
{

namespace model
{

/**
 * @brief Lumped rate model of liquid column chromatography without pores
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} + \frac{1 - \varepsilon_t}{\varepsilon_t} \frac{\partial q_{i}}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c, q)
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class LumpedRateModelWithoutPores : public ColumnLikeModel
{
public:

	LumpedRateModelWithoutPores(int unitOpIdx);
	virtual ~LumpedRateModelWithoutPores() CASEMA_NOEXCEPT;

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return false; }

	static const char* identifier() { return "LUMPED_RATE_MODEL_WITHOUT_PORES"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "LUMPED_RATE_MODEL_WITHOUT_PORES"; }

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:
	mpfr::mpreal _factor;

	virtual mpfr::mpreal phi(const mpfr::mpreal& s) const CASEMA_NOEXCEPT { return phiImpl(s); }
	virtual mpfr::mpcomplex phi(const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT { return phiImpl(s); }

	template <typename eval_t>
	eval_t phiImpl(const eval_t& s) const CASEMA_NOEXCEPT;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_LUMPEDRATEMODELWITHOUTPORES_HPP_
