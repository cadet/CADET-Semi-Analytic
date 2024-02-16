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

#include "model/InletModel.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <string>
#include <sstream>
#include <iomanip>

namespace casema
{

namespace model
{

InletModel::InletModel(int unitOpIdx) : UnitOperation(unitOpIdx)
{
}

InletModel::~InletModel() CASEMA_NOEXCEPT
{
}

bool InletModel::configure(io::IParameterProvider& paramProvider)
{
	UnitOperation::configure(paramProvider);

	_const.clear();
	_lin.clear();
	_quad.clear();
	_cub.clear();		

	unsigned int i = 0;
	std::ostringstream oss;
	oss << "sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
	while (paramProvider.exists(oss.str()))
	{
		paramProvider.pushScope(oss.str());

		if (paramProvider.exists("CONST_COEFF"))
		{
			if (paramProvider.isArray("CONST_COEFF"))
				throw InvalidParameterException("Field CONST_COEFF has to be scalar");

			_const.push_back(paramProvider.getDouble("CONST_COEFF"));
		}
		else
			_const.push_back(0);

		if (paramProvider.exists("LIN_COEFF"))
		{
			if (paramProvider.isArray("LIN_COEFF"))
				throw InvalidParameterException("Field LIN_COEFF has to be scalar");

			_lin.push_back(paramProvider.getDouble("LIN_COEFF"));
		}
		else
			_lin.push_back(0);

		if (paramProvider.exists("QUAD_COEFF"))
		{
			if (paramProvider.isArray("QUAD_COEFF"))
				throw InvalidParameterException("Field QUAD_COEFF has to be scalar");

			_quad.push_back(paramProvider.getDouble("QUAD_COEFF"));
		}
		else
			_quad.push_back(0);

		if (paramProvider.exists("CUBE_COEFF"))
		{
			if (paramProvider.isArray("CUBE_COEFF"))
				throw InvalidParameterException("Field CUBE_COEFF has to be scalar");

			_cub.push_back(paramProvider.getDouble("CUBE_COEFF"));
		}
		else
			_cub.push_back(0);

		paramProvider.popScope();

		++i;
		oss.str("");
		oss << "sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
	}

	return true;
}

void InletModel::setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT
{
	std::vector<mpfr::mpreal> v(secTimes, secTimes + nTimes);
	_sectionTimes = std::move(v);
}

template <typename eval_t>
eval_t InletModel::evalSection(int i, const eval_t& s) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal& Q = _sectionTimes[i];
	const mpfr::mpreal& T = _sectionTimes[i+1];
	const mpfr::mpreal dt = T - Q;
	const mpfr::mpreal& a = _const[i];
	const mpfr::mpreal& b = _lin[i];
	const mpfr::mpreal& c = _quad[i];
	const mpfr::mpreal& d = _cub[i];

	const eval_t ff = (a + (b + (2*c + 6 * d / s) / s) / s) / s;
	const eval_t tf = dt * ( (b + (2 * c + 6 * d / s) / s) / s + dt * (((3 / s + dt) * d  + c) / s)) + (a + (b + (2*c + 6 * d / s) / s) / s) / s;
	return ff * exp(-s*Q) - tf * exp(-s*T);
}

void InletModel::evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
//	h(0,0) = 0.0;

	g(0) = 0.0;
	for (int i = 0; i < _sectionTimes.size() - 1; ++i)
		g(0) += evalSection(i, s);
}

void InletModel::evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	evaluate(s, h, g);
}

bool InletModel::hasValidEstimate() const CASEMA_NOEXCEPT
{
	return true;
}

mpfr::mpreal InletModel::estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT
{
	mpfr::mpreal val(0);
	for (int i = 0; i < _sectionTimes.size() - 1; ++i)
	{
		const mpfr::mpreal delta = _sectionTimes[i+1] - _sectionTimes[i];
		const mpfr::mpreal upper = abs(_const[i] + delta * (_lin[i] + delta * (_quad[i] + delta * _cub[i]))) 
			+ (6 * abs(_cub[i]) / abscissa + abs(6 * _cub[i] * delta + 2 * _quad[i])) / sqr(abscissa)
			+ abs(_lin[i] + delta * (2 * _quad[i] + 3 * delta * _cub[i])) / abscissa;
		const mpfr::mpreal lower = abs(_const[i]) + (abs(_lin[i]) + (2 * abs(_quad[i]) + 6 * abs(_cub[i]) / abscissa) / abscissa) / abscissa;
		val += exp(-abscissa * _sectionTimes[i+1]) * upper + exp(-abscissa * _sectionTimes[i]) * lower;
	}

	return val;
}

mpfr::mpreal InletModel::timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT
{
	// Find maximum of nonnegative cubic spline
	mpfr::mpreal maxVal(0);
	for (int i = 0; i < _sectionTimes.size()-1; ++i)
	{
		// Assume that spline piece is nonnegative, compute zeros of first derivative
		// We have f(t) = cub * (t-t0)^3 + quad * (t-t0)^2 + lin * (t-t0) + cons
		// Therefore, f'(t) = 3 * cub * (t-t0)^2 + 2 * quad * (t-t0) + lin
		// f'(t) = 0  <=>  t-t0 = (-quad +/- sqrt( quad^2 - 3 * lin * cub) ) / (3 * cub)
		
		// Stable solution
		const mpfr::mpreal radical = sqr(_quad[i]) - 3 * _lin[i] * _cub[i];
		mpfr::mpreal maxPiece = 0;
		if (radical > 0)
		{
			// Two real solutions
			
			const mpfr::mpreal delta1 = (_quad[i] < 0) ? (-_quad[i] + sqrt(radical)) / (3 * _cub[i]) : (-_quad[i] - sqrt(radical)) / (3 * _cub[i]);
			const mpfr::mpreal delta2 = _lin[i] / (3 * _cub[i] * delta1);

			const mpfr::mpreal v1 = _const[i] + delta1 * (_lin[i] + delta1 * (_quad[i] + delta1 * _cub[i]));
			const mpfr::mpreal v2 = _const[i] + delta2 * (_lin[i] + delta2 * (_quad[i] + delta2 * _cub[i]));
			maxPiece = max(mpfr::mpreal(0), max(v1, v2));
		}
		else if (abs(radical) <= std::numeric_limits<mpfr::mpreal>::epsilon())
		{
			// One real solution
			const mpfr::mpreal delta = -_quad[i] / (3 * _cub[i]);
			maxPiece = max(mpfr::mpreal(0), _const[i] + delta * (_lin[i] + delta * (_quad[i] + delta * _cub[i])));
		}

		// Check boundaries
		const mpfr::mpreal delta = _sectionTimes[i+1] - _sectionTimes[i];
		const mpfr::mpreal right = _const[i] + delta * (_lin[i] + delta * (_quad[i] + delta * _cub[i]));
		maxPiece = max(maxPiece, max(_const[i], right));

		maxVal = max(maxVal, maxPiece);
	}
	return maxVal;
}


void registerInletModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[InletModel::identifier()] = [](int uoId) { return new InletModel(uoId); };
}

/*

void inletConcentration(double t, unsigned int sec, double* inletConc)
{
	// This function evaluates a piecewise cubic polynomial given on some intervals
	// called sections. On each section a polynomial of degree 3 is evaluated:
	// 
	//   p_i(t) = CUBE_COEFF[i] * (t - t_i)^3 + QUAD_COEFF[i] * (t - t_i)^2 + LIN_COEFF[i] * (t - t_i) + CONST_COEFF[i],
	//   
	// where p_i is the polynomial on section i given by the interval [t_i, t_{i+1}].

	casema_assert(sec < _sectionTimes.size());

	const double tShift = t - _sectionTimes[sec];
	const unsigned int wrapSec = sec % (_const.size()/_nComp);

	double const* const con = _const.data() + wrapSec * _nComp;
	double const* const lin = _lin.data() + wrapSec * _nComp;
	double const* const quad = _quad.data() + wrapSec * _nComp;
	double const* const cub = _cub.data() + wrapSec * _nComp;

	// Evaluate polynomial using Horner's scheme
	for (unsigned int comp = 0; comp < _nComp; ++comp)
		inletConc[comp] = con[comp] + tShift * (lin[comp] + tShift * (quad[comp] + tShift * cub[comp]));
}
*/

}  // namespace model

}  // namespace casema
