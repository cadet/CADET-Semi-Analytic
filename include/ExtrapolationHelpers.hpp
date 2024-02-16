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

#ifndef CASEMA_EXTRAPOLATION_HELPERS_HPP_
#define CASEMA_EXTRAPOLATION_HELPERS_HPP_

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <limits>
#include "Extrapolation.hpp"

namespace casema
{
	template <typename real_t>
	class ConvergenceMonitoredExtrapolationMethod : public ExtrapolationMethod<real_t>
	{
	public:
		// Threshold for convergence test
		virtual const real_t& thresholdAbs() const = 0;
		virtual real_t thresholdAbs() = 0;
		virtual void thresholdAbs(const real_t& t) = 0;

		virtual const real_t& thresholdRel() const = 0;
		virtual real_t thresholdRel() = 0;
		virtual void thresholdRel(const real_t& t) = 0;

		// Number of times convergence test has to be successful
		virtual const std::size_t& times() const = 0;
		virtual std::size_t times() = 0;
		virtual void times(const std::size_t& t) = 0;

		virtual bool converged() const = 0;
	};

	template <typename real_t, template<class T> class ExtMethod>
	class ConvergenceMonitor : public ConvergenceMonitoredExtrapolationMethod<real_t>
	{
	public:
		ConvergenceMonitor() : _thresholdRel(0), _thresholdAbs(0), _times(2), _success(0), _lastVal(0) { _lastVal.setNan(); }
		ConvergenceMonitor(const ExtMethod<real_t>& extrapolator) : _extrapolator(extrapolator), _thresholdRel(0), _thresholdAbs(0), _times(2), _success(0), _lastVal(0) { _lastVal.setNan(); }

		// Threshold for convergence test
		virtual const real_t& thresholdAbs() const { return _thresholdAbs; }
		virtual real_t thresholdAbs() { return _thresholdAbs; }
		virtual void thresholdAbs(const real_t& t) { _thresholdAbs = t; }

		virtual const real_t& thresholdRel() const { return _thresholdRel; }
		virtual real_t thresholdRel() { return _thresholdRel; }
		virtual void thresholdRel(const real_t& t) { _thresholdRel = t; }

		// Number of times convergence test has to be successful
		virtual const std::size_t& times() const { return _times; }
		virtual std::size_t times() { return _times; }
		virtual void times(const std::size_t& t) { _times = t; }

		bool converged() const { return _success >= _times; }

		virtual real_t next(const real_t& val)
		{
			const real_t est = _extrapolator.next(val);
			checkConvergence(est);			
			return est;
		}

		virtual real_t nextPoint(const real_t& val, const real_t& x)
		{
			const real_t est = _extrapolator.nextPoint(val, x);
			checkConvergence(est);			
			return est;
		}

		static const char* const compileTimeName() { return ExtMethod<real_t>::compileTimeName(); }
		virtual const char* const name() const { return _extrapolator.name(); }
		virtual std::size_t iterations() const { return _extrapolator.iterations(); }

		virtual void reset()
		{
			_extrapolator.reset();

			_lastVal.setNan();
			_success = 0;
		}

	protected:

		void checkConvergence(const real_t& est)
		{
			if (!isnan(_lastVal))
			{
				// Not the first iteration, so perform convergence check
				const real_t diff = abs(_lastVal - est);
				if ((diff <= abs(_lastVal) * _thresholdRel) || (diff <= _thresholdAbs))
					++_success;
				else
					_success = 0;
			}
			_lastVal = est;
		}

		ExtMethod<real_t> _extrapolator;
		real_t _thresholdRel;
		real_t _thresholdAbs;
		std::size_t _times;

		std::size_t _success;
		real_t _lastVal;
	};


	namespace detail
	{
		template <typename real_t>
		void enclosingInterval(const std::vector<real_t>& ests, real_t& center, real_t& radius)
		{
			real_t minVal = std::numeric_limits<real_t>::max();
			real_t maxVal = std::numeric_limits<real_t>::lowest();

			for (typename std::vector<real_t>::const_iterator it = ests.begin(); it != ests.end(); ++it)
			{
				minVal = std::min(*it, minVal);
				maxVal = std::max(*it, maxVal);
			}

			center = (minVal + maxVal) * real_t(0.5);
			radius = (maxVal - minVal) * real_t(0.5);
		}

		template <typename real_t>
		void enclosingInterval(const std::vector<real_t>& ests, const std::vector<std::size_t>& idx, real_t& center, real_t& radius)
		{
			real_t minVal = std::numeric_limits<real_t>::max();
			real_t maxVal = std::numeric_limits<real_t>::lowest();

			for (std::size_t i = 0; i < idx.size() - 1; ++i)
			{
				const real_t& v = ests[idx[i]];
				minVal = std::min(v, minVal);
				maxVal = std::max(v, maxVal);
			}

			center = (minVal + maxVal) * real_t(0.5);
			radius = (maxVal - minVal) * real_t(0.5);
		}
	}

	template <typename real_t>
	real_t consensus(const std::vector<real_t>& ests, std::size_t mustAgree, real_t* radius = nullptr)
	{
		if (ests.empty())
		{
			if (radius) 
				*radius = 0;
			return real_t(0);
		}

		if (ests.size() == 1)
		{
			if (radius) 
				*radius = 0;
			return ests[0];
		}

		if (ests.size() <= mustAgree)
		{
			real_t curCenter;
			real_t curRadius;
			detail::enclosingInterval(ests, curCenter, curRadius);
			if (radius)
				*radius = curRadius;
			return curCenter;
		}

		std::size_t i;
		std::vector<std::size_t> bitpos(mustAgree + 1, 0);

		// Seed with the lowest word plus a sentinel
		for(i = 0; i < mustAgree; ++i)
			bitpos[i] = i;
		bitpos[i] = 0;

		real_t curCenter;
		real_t curRadius = std::numeric_limits<real_t>::max();

		do {
			real_t centerCand;
			real_t radiusCand;
			detail::enclosingInterval(ests, bitpos, centerCand, radiusCand);

			if (radiusCand < curRadius)
			{
				curRadius = radiusCand;
				curCenter = centerCand;
			}

			// Increment the least-significant series of consecutive bits
			for(i = 0; bitpos[i + 1] == bitpos[i] + 1; ++i)
				bitpos[i] = i;
		} while(++bitpos[i] != ests.size());

		if (radius)
			*radius = curRadius;
		return curCenter;
	}


	template <typename real_t>
	class ConsensusEstimator : public ExtrapolationMethod<real_t>
	{
	public:

		typedef typename std::vector<ExtrapolationMethod<real_t>*>::const_iterator const_iterator;
		typedef typename std::vector<ExtrapolationMethod<real_t>*>::iterator iterator;

		ConsensusEstimator() { }
		virtual ~ConsensusEstimator()
		{
			for (iterator it = _methods.begin(); it != _methods.end(); ++it)
				delete (*it);			
		}

		const std::size_t& mustAgree() const { return _mustAgree; }
		std::size_t mustAgree() { return _mustAgree; }
		void mustAgree(std::size_t ma) { _mustAgree = ma; }

		void add(ExtrapolationMethod<real_t>& extrapolator) { _methods.push_back(&extrapolator); }
		void add(ExtrapolationMethod<real_t>* const extrapolator) { _methods.push_back(extrapolator); }

		virtual real_t next(const real_t& val)
		{
			if (_estimates.size() != _methods.size())
			{
				// Reserve space
				_estimates.resize(_methods.size());
			}

			for (std::size_t i = 0; i < _methods.size(); ++i)
			{
				_estimates[i] = _methods[i]->next(val);
			}
			return consensus(_estimates);
		}

		virtual real_t nextPoint(const real_t& val, const real_t& x)
		{
			if (_estimates.size() != _methods.size())
			{
				// Reserve space
				_estimates.resize(_methods.size());
			}

			for (std::size_t i = 0; i < _methods.size(); ++i)
			{
				_estimates[i] = _methods[i]->nextPoint(val, x);
			}
			return consensus(_estimates);
		}

		virtual const char* const name() const
		{
			std::ostringstream oss;
			if (_methods.size() > 0)
			{
				oss << "{" << _methods[0]->name();

				for (const_iterator it = _methods.begin() + 1; it != _methods.end(); ++it)
					oss << ", " << (*it)->name();

				oss << "} ";
			}
			oss << "Consensus";
			return oss.str().c_str();
		}

		virtual std::size_t iterations() const 
		{
			std::size_t iter = 0;
			for (const_iterator it = _methods.begin(); it != _methods.end(); ++it)
			{
				iter = std::max(iter, (*it)->iterations());
			}
			return iter;
		}

		virtual void reset()
		{
			for (iterator it = _methods.begin(); it != _methods.end(); ++it)
				(*it)->reset();
		}

		const real_t& center() const { return _curCenter; }
		const real_t& radius() const { return _curRadius; }

		const size_t numEstimators() const { return _methods.size(); }

		iterator begin() { return _methods.begin(); }
		const_iterator begin() const { return _methods.begin(); }
		iterator end() { return _methods.end(); }
		const_iterator end() const { return _methods.end(); }

		const std::vector<real_t>& estimates() const { return _estimates; }
		const std::vector<ExtrapolationMethod<real_t>*>& estimators() const { return _methods; }

	protected:

		real_t consensus(const std::vector<real_t>& ests)
		{
			_curCenter = ::casema::consensus<real_t>(ests, _mustAgree, &_curRadius);
			return _curCenter;
		}

		std::vector<ExtrapolationMethod<real_t>*> _methods;
		std::vector<real_t> _estimates;
		std::size_t _mustAgree;
		real_t _curRadius;
		real_t _curCenter;
	};
}

#endif
