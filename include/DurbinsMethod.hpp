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

#ifndef CASEMA_DURBINSMETHOD_HPP_
#define CASEMA_DURBINSMETHOD_HPP_

#include "casemaCompilerInfo.hpp"
#include "MPComplex.hpp"
#include "MPReal.hpp"
#include "SlicedVector.hpp"

#include <functional>
#include <vector>

namespace casema
{
	util::SlicedVector<mpfr::mpreal> invertLaplace(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime);
	util::SlicedVector<mpfr::mpreal> invertLaplaceKahan(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime);

	util::SlicedVector<mpfr::mpreal> invertLaplace(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, const std::function<void(double, int)>& callback);
	util::SlicedVector<mpfr::mpreal> invertLaplaceKahan(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, const std::function<void(double, int)>& callback);
}

#endif
