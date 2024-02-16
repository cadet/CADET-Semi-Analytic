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
 * Defines the general rate model (GRM).
 */

#ifndef CASEMA_GENERALRATEMODEL_HPP_
#define CASEMA_GENERALRATEMODEL_HPP_

#include "model/ColumnLikeModel.hpp"

#include <vector>

namespace casema
{

namespace model
{

/**
 * @brief General rate model of liquid column chromatography
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= D_{p,i} \left( \frac{\partial^2 c_{p,i}}{\partial r^2} + \frac{2}{r} \frac{\partial c_{p,i}}{\partial r} \right) + D_{s,i} \frac{1 - \varepsilon_p}{\varepsilon_p} \left( \frac{\partial^2 q_{i}}{\partial r^2} + \frac{2}{r} \frac{\partial q_{i}}{\partial r} \right) \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
\end{align} @f]
@f[ \begin{align}
	j_{f,i} = k_{f,i} \left( c_i - c_{p,i} \left(\cdot, \cdot, r_p\right)\right)
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0 \\
\varepsilon_p D_{p,i} \frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, r_p) + (1-\varepsilon_p) D_{s,i} \frac{\partial q_{i}}{\partial r}(\cdot, \cdot, r_p) &= j_{f,i} \\
\frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, 0) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class GeneralRateModel : public ColumnWithPoreDiffusion
{
public:

	GeneralRateModel(int unitOpIdx);
	virtual ~GeneralRateModel() CASEMA_NOEXCEPT;

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return false; }

	static const char* identifier() { return "GENERAL_RATE_MODEL"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "GENERAL_RATE_MODEL"; }

	virtual bool configure(io::IParameterProvider& paramProvider);

	virtual bool hasValidEstimate() const CASEMA_NOEXCEPT;

protected:

};

} // namespace model
} // namespace casema

#endif  // CASEMA_GENERALRATEMODEL_HPP_
