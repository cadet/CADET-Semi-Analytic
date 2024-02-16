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
 * Defines the lumped rate model with pores (LRMP).
 */

#ifndef CASEMA_LUMPEDRATEMODELWITHPORES_HPP_
#define CASEMA_LUMPEDRATEMODELWITHPORES_HPP_

#include "model/ColumnLikeModel.hpp"

#include <array>
#include <vector>

namespace casema
{

namespace model
{

/**
 * @brief Lumped rate model of liquid column chromatography with pores
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3 k_{f,i}}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= \frac{3 k_{f,i}}{\varepsilon_p r_p} j_{f,i} \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
\end{align} @f]
@f[ \begin{align}
	j_{f,i} = c_i - c_{p,i}
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class LumpedRateModelWithPores : public ColumnWithParticles
{
public:

	LumpedRateModelWithPores(int unitOpIdx);
	virtual ~LumpedRateModelWithPores() CASEMA_NOEXCEPT;

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return false; }

	static const char* identifier() { return "LUMPED_RATE_MODEL_WITH_PORES"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "LUMPED_RATE_MODEL_WITH_PORES"; }

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:

	virtual mpfr::mpreal particleF(int j, const mpfr::mpreal& s) const CASEMA_NOEXCEPT { return particleFimpl(j, s); }
	virtual mpfr::mpcomplex particleF(int j, const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT { return particleFimpl(j, s); }

	template <typename eval_t>
	eval_t particleFimpl(int j, const eval_t& s) const CASEMA_NOEXCEPT;

	std::vector<mpfr::mpreal> _threeKfOverRp;
	std::vector<mpfr::mpreal> _parFactor;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_LUMPEDRATEMODELWITHPORES_HPP_
