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

#ifndef CASEMA_EIGEN_HPP_
#define CASEMA_EIGEN_HPP_

#include "MPReal.hpp"
#include "MPComplex.hpp"

#include <unsupported/Eigen/MPRealSupport>
#include <Eigen/Core>

namespace Eigen
{
	template<> struct NumTraits<mpfr::mpcomplex>
		: NumTraits<mpfr::mpreal> // GenericNumTraits<mpfr::mpcomplex>
	{
		typedef mpfr::mpreal Real;
		typedef mpfr::mpcomplex NonInteger;
		typedef mpfr::mpcomplex Nested;
		enum 
		{
			IsComplex = 1,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1,
			ReadCost = HugeCost,
			AddCost = HugeCost,
			MulCost = HugeCost
		};
	};

	typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXmp;
	typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> VectorXmp;

	typedef Eigen::Matrix<mpfr::mpcomplex, Eigen::Dynamic, Eigen::Dynamic> MatrixCmp;
	typedef Eigen::Matrix<mpfr::mpcomplex, Eigen::Dynamic, 1> VectorCmp;
}

#endif  // CASEMA_EIGEN_HPP_
