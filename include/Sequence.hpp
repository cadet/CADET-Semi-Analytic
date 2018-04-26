// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_SEQUENCES_HPP_
#define CASEMA_SEQUENCES_HPP_

#include <vector>

namespace casema
{
	template <typename T>
	std::vector<T> linearSequence(const T& start, const T& end, const unsigned int steps)
	{
		std::vector<T> seq(steps, start);
		for (unsigned int i = 1; i < steps; ++i)
			seq[i] = start + (end-start) / (steps-1) * i;
		return seq;
	}

	template <typename T>
	std::vector<T> logLinearSequence(const T& start, const T& end, const unsigned int steps)
	{
		const T ls = log(start);
		const T le = log(end);
		std::vector<T> seq(steps, ls);
		for (unsigned int i = 1; i < steps; ++i)
			seq[i] = exp(ls + (le-ls) / (steps-1) * i);
		return seq;
	}

	template <typename T>
	std::vector<T> geometricSequence(const T& start, const T& factor, const unsigned int steps)
	{
		std::vector<T> seq(steps, start);
		for (unsigned int i = 1; i < steps; ++i)
			seq[i] = seq[i-1] * factor;
		return seq;		
	}


	template <typename real_t>
	class Sequence
	{
	public:
		Sequence() { }
		virtual ~Sequence() { }

		virtual real_t next() = 0;
		virtual const char* name() const = 0;
	};


	template <typename real_t>
	class LinearSequence : public Sequence<real_t>
	{
	public:
		LinearSequence(const real_t& start, const real_t& slope) : _cur(start), _slope(slope) { }

		virtual real_t next()
		{
			const real_t res = _cur;
			_cur -= _slope;
			return res;
		}

		virtual const char* name() const { return "Linear"; }

	protected:
		real_t _cur;
		real_t _slope;
	};


	template <typename real_t>
	class LogLinearSequence : public Sequence<real_t>
	{
	public:
		LogLinearSequence(const real_t& start, const real_t& slope) : _cur(log(start)), _slope(slope) { }

		virtual real_t next()
		{
			const real_t res = _cur;
			_cur -= _slope;
			return exp(res);
		}

		virtual const char* name() const { return "Logarithmic Linear"; }

	protected:
		real_t _cur;
		real_t _slope;
	};


	template <typename real_t>
	class GeometricSequence : public Sequence<real_t>
	{
	public:
		GeometricSequence(const real_t& start, const real_t& factor) : _cur(start), _factor(factor) { }

		virtual real_t next()
		{
			const real_t res = _cur;
			_cur *= _factor;
			return res;
		}

		virtual const char* name() const { return "Geometric"; }

	protected:
		real_t _cur;
		real_t _factor;
	};
}

#endif
