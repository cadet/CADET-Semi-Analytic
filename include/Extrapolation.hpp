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
 * General comments
 * ----------------
 *
 * This file contains several extrapolation algorithms all following the same
 * interface ExtrapolationMethod. There are generally two types of methods: 
 *   (a) Methods only requiring a function value sequence
 *   (b) Methods that require a function value sequence and an evaluation
 *		 point sequence
 * The evaluation point sequences (x_n) are assumed to be strict monotonically
 * decreasing and to converge to 0, i.e., lim_{n -> oo} = 0. If a method
 * requires strict monotonically increasing sequences that diverge to infinity,
 * it uses 1 / x_n instaead of x_n.
 *
 * All extrapolation methods can be used in online mode, i.e., data is added
 * iteratively via the next() and nextPoint() methods. Methods of type (a)
 * ignore the additional evaluation point given in nextPoint() and relay the
 * call to next().
 **/

#ifndef CASEMA_EXTRAPOLATION_HPP_
#define CASEMA_EXTRAPOLATION_HPP_

#include <vector>

namespace casema
{
	template <typename T>
	class ExtrapolationMethod
	{
	public:
		ExtrapolationMethod() { }
		virtual ~ExtrapolationMethod() { }

		virtual T next(const T& val) = 0;
		virtual T nextPoint(const T& val, const T& pt) { return this->next(val); }

		virtual std::vector<T> next(const std::vector<T>& vals)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = this->next(vals[i]);

			return estimates;			
		}

		virtual std::vector<T> nextPoint(const std::vector<T>& val, const std::vector<T>& pt) { return this->next(val); }

		virtual std::size_t iterations() const = 0;
		virtual const char* const name() const = 0;

		virtual void reset() = 0;
	};


	template <typename T>
	class IdentityExtrapolation : public ExtrapolationMethod<T>
	{
	public:
		IdentityExtrapolation() : _n(0) { }

		virtual T next(const T& val)
		{
			++_n;
			return val;
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { _n += vals.size(); return vals; }

		virtual std::size_t iterations() const { return _n; }
		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Identity"; }

		virtual void reset() { _n = 0; }

	protected:
		std::size_t _n;
	};


	template <typename T>
	class WynnEpsilonMethod : public ExtrapolationMethod<T>
	{
	public:
		WynnEpsilonMethod() { }
		WynnEpsilonMethod(std::size_t seqLength) { _diag.reserve(seqLength); }

		virtual T next(const T& val)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf
			
			_diag.push_back(val);

			// We need at least two elements here, so return if this is the first call
			if (_diag.size() == 1)
				return val;

			// Update counter diagonal from bottom to top
			T aux2(0);
			for (std::size_t j = _diag.size()-1; j > 0; --j)
			{
				const T aux1 = aux2;
				aux2 = _diag[j-1];

				// Guard agains under- and overflow
				const T diff = _diag[j] - aux2;
				if (abs(diff) <= std::numeric_limits<T>::min())
					_diag[j-1] = aux1 + T(1) / std::numeric_limits<T>::max();
				else
					_diag[j-1] = aux1 + T(1) / diff;
			}

			// Select estimate depending on the number of data points being odd or even
			if (_diag.size() % 2 == 1)
				return _diag[0];
			else
				return _diag[1];
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual std::size_t iterations() const { return _diag.size(); }
		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Wynn's Epsilon Algorithm"; }

		virtual void reset() { _diag.clear(); }

	protected:
		std::vector<T> _diag;
	};

	template <typename T>
	std::vector<T> wynnEpsilon(const std::vector<T>& seq)
	{
		WynnEpsilonMethod<T> wem(seq.size());
		return wem.next(seq);
	}


	template <typename T>
	class WynnRhoMethod : public ExtrapolationMethod<T>
	{
	public:
		WynnRhoMethod() { }
		WynnRhoMethod(std::size_t seqLength)
		{
			_diag.reserve(seqLength);
			_pts.reserve(seqLength);
		}

		virtual T nextPoint(const T& val, const T& pt)
		{
			// We need divergence to infinity instead of converge to 0
			const T invPt = T(1) / pt;
			return calculateNext(val, invPt);
		}

		virtual T next(const T& val)
		{
			return calculateNext(val, T(_diag.size()));
		}

		virtual std::vector<T> nextPoint(const std::vector<T>& vals, const std::vector<T>& pts)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = nextPoint(vals[i], pts[i]);

			return estimates;
		}

		virtual std::vector<T> next(const std::vector<T>& vals)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = nextPoint(vals[i], T(1) / T(i + _diag.size()));

			return estimates;
		}

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Wynn's Rho Algorithm"; }
		virtual std::size_t iterations() const { return _diag.size(); }

		virtual void reset()
		{
			_diag.clear();
			_pts.clear();
		}

	protected:

		T calculateNext(const T& val, const T& pt)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf
			
			_diag.push_back(val);
			_pts.push_back(pt);

			// We need at least two elements here, so return if this is the first call
			if (_diag.size() == 1)
				return val;

			// Update counter diagonal from bottom to top
			T aux2(0);
			for (std::size_t j = _diag.size()-1; j > 0; --j)
			{
				const T aux1 = aux2;
				aux2 = _diag[j-1];

				// Guard agains under- and overflow
				const T diff = _diag[j] - aux2;
				if (abs(diff) <= std::numeric_limits<T>::min())
					_diag[j-1] = aux1 + (pt - _pts[j-1]) / std::numeric_limits<T>::max();
				else
					_diag[j-1] = aux1 + (pt - _pts[j-1]) / diff;
			}

			// Select estimate depending on the number of data points being odd or even
			if (_diag.size() % 2 == 1)
				return _diag[0];
			else
				return _diag[1];
		}

		std::vector<T> _diag;
		std::vector<T> _pts;
	};

	template <typename T>
	std::vector<T> wynnRho(const std::vector<T>& seq)
	{
		WynnRhoMethod<T> wrm(seq.size());
		return wrm.next(seq);
	}

	template <typename T>
	std::vector<T> wynnRhoPoints(const std::vector<T>& seq, const std::vector<T>& pts)
	{
		WynnRhoMethod<T> wrm(seq.size());
		return wrm.next(seq, pts);
	}


	template <typename T>
	class IteratedAitkenDeltaSquaredMethod : public ExtrapolationMethod<T>
	{
	public:
		IteratedAitkenDeltaSquaredMethod() { }
		IteratedAitkenDeltaSquaredMethod(std::size_t seqLength) { _diag.reserve(seqLength); }

		virtual T next(const T& val)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf
			
			_diag.push_back(val);

			// We need at least three elements here, so return if this is the first or second call
			if (_diag.size() <= 2)
				return val;

			const std::size_t N = _diag.size()-1;

			for (std::size_t j = 1; j <= N / 2; ++j)
			{
				const std::size_t M = N - 2 * j;

				// Guard against under- and overflow
				const T denom = _diag[M+2] - T(2) * _diag[M+1] + _diag[M];
				if (abs(denom) <= std::numeric_limits<T>::min())
					_diag[M] = _diag[M] - sqr(_diag[M] - _diag[M+1]) / std::numeric_limits<T>::max();
				else
					_diag[M] = _diag[M] - sqr(_diag[M] - _diag[M+1]) / denom;
			}

			// Select estimate depending on the number of data points being odd or even
			if (N % 2 == 0)
				return _diag[0];
			else
				return _diag[1];
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Iterated Aitken's Delta Squared"; }
		virtual std::size_t iterations() const { return _diag.size(); }

		virtual void reset() { _diag.clear(); }

	protected:
		std::vector<T> _diag;
	};

	template <typename T>
	std::vector<T> iteratedAitkenDeltaSquared(const std::vector<T>& seq)
	{
		IteratedAitkenDeltaSquaredMethod<T> iadsm(seq.size());
		return iadsm.next(seq);
	}


	template <typename T>
	class AitkenDeltaSquaredMethod : public ExtrapolationMethod<T>
	{
	public:
		AitkenDeltaSquaredMethod() : _last2(0), _last(0), _n(0) { }

		virtual T next(const T& val)
		{
			++_n;

			if (_n <= 2)
			{
				_last2 = _last;
				_last = val;
				return val;
			}

			const T dx = _last - _last2;
			const T dx2 = _last2 - T(2) * _last + val;

			// Guard against under- and overflow
			T ret;
			if (abs(dx2) <= std::numeric_limits<T>::min())
				ret = _last2 - sqr(dx) / std::numeric_limits<T>::max();
			else
				ret = _last2 - sqr(dx) / dx2;

			_last2 = _last;
			_last = val;
			return ret;			
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Aitken's Delta Squared"; }
		virtual std::size_t iterations() const { return _n; }

		virtual void reset()
		{
			_n = 0;
			_last2 = 0;
			_last = 0;
		}

	protected:
		T _last2;
		T _last;
		std::size_t _n;
	};

	template <typename T>
	std::vector<T> aitkenDeltaSquared(const std::vector<T>& seq)
	{
		AitkenDeltaSquaredMethod<T> adsm;
		return adsm.next(seq);
	}

	template <typename T>
	std::vector<T> aitkenDeltaSquaredManual(const std::vector<T>& seq)
	{
		std::vector<T> out(seq.size()-2, T());
		for (unsigned int i = 0; i < out.size(); ++i)
		{
			const T dx = seq[i+1] - seq[i];
			const T dx2 = seq[i] - T(2) * seq[i+1] + seq[i+2];

			// Guard against under- and overflow
			if (abs(dx2) <= std::numeric_limits<T>::min())
				out[i] = seq[i] - sqr(dx) / std::numeric_limits<T>::max();
			else
				out[i] = seq[i] - sqr(dx) / dx2;
		}
		return out;
	}


	template <typename T>
	class LevinTypeMethod : public ExtrapolationMethod<T>
	{
	public:
		LevinTypeMethod() { }
		LevinTypeMethod(std::size_t seqLength) 
		{
			_num.reserve(seqLength);
			_den.reserve(seqLength);
		}

		virtual T next(const T& val, const T& beta, const T& remainder)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf

			_num.push_back(val / remainder);
			_den.push_back(T(1) / remainder);

			// We need at least two elements here, so return if this is the first call
			if (_num.size() <= 1)
				return val;

			const std::size_t N = _num.size()-1;

			if (N <= 1)
			{
				// We have two elements
				_num[N-1] = _num[N] - _num[N-1];
				_den[N-1] = _den[N] - _den[N-1];
				return _num[0] / _den[0];
			}

			// We have at least three elements
			const T bn1 = beta + T(N-1);
			const T bn2 = beta + T(N);
			const T coeff = bn1 / bn2;

			// Update numerator and denominator
			for (std::size_t j = 2; j <= N; ++j)
			{
				const T fact = (beta + T(N-j)) * pow(coeff, T(j-2)) / bn2;
				_num[N-j] = _num[N-j+1] - fact * _num[N-j];
				_den[N-j] = _den[N-j+1] - fact * _den[N-j];
			}

			// Take care of under- or overflow
			if (abs(_den[0]) <= std::numeric_limits<T>::min())
				return _num[0] / std::numeric_limits<T>::max();
			else
				return _num[0] / _den[0];

/*
			// Implementation of algorithm in Numerical Recipes 3rd ed.
			T term = T(1) / (beta + T(_num.size()));
			_den.push_back(term / remainder);
			_num.push_back(val * _den[_den.size() - 1]);
			const std::size_t n = _num.size()-1;
			if (_num.size() > 0)
			{
				const T ratio = (beta + T(n-1)) * term;
				for (int j = 1; j <= n; ++j)
				{
					const T fact = (T(n-j) + beta) * term;
					_num[n-j] = _num[n-j+1] - fact * _num[n-j];
					_den[n-j] = _den[n-j+1] - fact * _den[n-j];
					term *= ratio;
				}
			}
			return _num[0] / _den[0];
*/
		}

		virtual std::vector<T> next(const std::vector<T>& vals, const T& beta, const std::vector<T>& remainders)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = next(vals[i], beta, remainders[i]);

			return estimates;
		}

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Levin Type Method"; }
		virtual std::size_t iterations() const { return _num.size(); }

		virtual void reset()
		{
			_num.clear();
			_den.clear();
		}

	protected:
		std::vector<T> _num;
		std::vector<T> _den;
	};

	template <typename T>
	class LevinUMethod : public LevinTypeMethod<T>
	{
	public:
		LevinUMethod() : _lastVal(1), _beta(1) { }
		LevinUMethod(std::size_t seqLength) : LevinTypeMethod<T>(seqLength), _lastVal(1), _beta(1) { }

		const T& beta() const { return _beta; }
		T beta() { return _beta; }
		void beta(const T& beta) { _beta = beta; }

		virtual T next(const T& val, const T& beta)
		{
			const T retVal = LevinTypeMethod<T>::next(val, beta, val - _lastVal);
			_lastVal = val;
			return retVal;
		}

		virtual std::vector<T> next(const std::vector<T>& vals, const T& beta)
		{
			std::vector<T> estimates(vals.size(), T());
			estimates[0] = LevinTypeMethod<T>::next(vals[0], beta, T(1));
			for (std::size_t i = 1; i < estimates.size(); ++i)
				estimates[i] = LevinTypeMethod<T>::next(vals[i], beta, vals[i] - vals[i-1]);

			_lastVal = vals[vals.size()-1];
			return estimates;
		}

		virtual T next(const T& val)
		{
			return next(val, _beta);
		}

		virtual std::vector<T> next(const std::vector<T>& vals)
		{
			return next(vals, _beta);
		}

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Levin U Method"; }

		virtual void reset()
		{
			LevinTypeMethod<T>::reset();
			_lastVal = T(1);
		}

	protected:
		T _lastVal;
		T _beta;
	};

	template <typename T>
	std::vector<T> levinUbeta(const std::vector<T>& seq, const T& beta)
	{
		LevinUMethod<T> lum(seq.size());
		return lum.next(seq, beta);
	}

	template <typename T>
	std::vector<T> levinU(const std::vector<T>& seq)
	{
		return levinUbeta(seq, T(1));
	}

	template <typename T>
	class LevinTMethod : public LevinTypeMethod<T>
	{
	public:
		LevinTMethod() { _lastVal = T(1); }
		LevinTMethod(std::size_t seqLength) : LevinTypeMethod<T>(seqLength) { _lastVal = T(1); }

		virtual T next(const T& val)
		{
			const T retVal = LevinTypeMethod<T>::next(val, T(0), val - _lastVal);
			_lastVal = val;
			return retVal;
		}

		virtual std::vector<T> next(const std::vector<T>& vals)
		{
			std::vector<T> estimates(vals.size()-1, T());
			estimates[0] = LevinTypeMethod<T>::next(vals[0], T(0), T(1));
			for (std::size_t i = 1; i < estimates.size(); ++i)
				estimates[i] = LevinTypeMethod<T>::next(vals[i], T(0), vals[i] - vals[i-1]);

			_lastVal = vals[vals.size()-1];
			return estimates;
		}

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Levin T Method"; }

		virtual void reset()
		{
			LevinTypeMethod<T>::reset();
			_lastVal = T(1);
		}

	protected:
		T _lastVal;
	};

	template <typename T>
	std::vector<T> levinT(const std::vector<T>& seq)
	{
		LevinTMethod<T> lum(seq.size());
		return lum.next(seq);
	}


	template <typename T>
	class BrezinskiThetaMethod : public ExtrapolationMethod<T>
	{
	public:
		BrezinskiThetaMethod() : _N(0) { }
		BrezinskiThetaMethod(std::size_t seqLength)
		{
			_N = 0;
			_a.reserve(seqLength);
			_b.reserve(seqLength);
		}

		virtual T next(const T& val)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf

			const std::size_t jmax = (_N * 2 + 1) / 3;
			const std::size_t nMod2 = _N % 2;

			if (_N == 0)
			{
				++_N;
				_a.push_back(val);
				return val;
			}
			++_N;
			
			if (nMod2 == 0)
				return computeNext(jmax, val, _a, _b);
			else
				return computeNext(jmax, val, _b, _a);
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Brezinski's Theta Method"; }
		virtual std::size_t iterations() const { return _N; }

		virtual void reset()
		{
			_N = 0;
			_a.clear();
			_b.clear();
		}

	protected:

		T computeNext(const std::size_t jmax, const T& val, std::vector<T>& a, std::vector<T>& b)
		{
			while (a.size() <= jmax+1)
				a.push_back(T(0));
			while (b.size() <= jmax+1)
				b.push_back(T(0));

			T aux2 = T(0);
			T aux1 = a[0];
			a[0] = val;
			for (std::size_t j = 1; j <= jmax; ++j)
			{
				const T aux3 = aux2;
				aux2 = aux1;
				if (j < jmax)
					aux1 = a[j];

				if (j % 2 == 0)
				{
					const T denom = a[j-1] - T(2) * b[j-1] + aux2;

					// Take care of under- or overflow
					if (abs(denom) <= std::numeric_limits<T>::min())
						a[j] = aux3 + (b[j-2] - aux3) * (a[j-1] - b[j-1]) / std::numeric_limits<T>::max();
					else
						a[j] = aux3 + (b[j-2] - aux3) * (a[j-1] - b[j-1]) / denom;
				}
				else
				{
					const T diff = a[j-1] - b[j-1];

					// Take care of under- or overflow
					if (abs(diff) <= std::numeric_limits<T>::min())
						a[j] = aux3 + T(1) / std::numeric_limits<T>::max();
					else
						a[j] = aux3 + T(1) / diff;
				}
			}
			if (jmax % 2 == 0)
				return a[jmax];
			else
				return a[jmax-1];
		}

		std::size_t _N;
		std::vector<T> _a;
		std::vector<T> _b;
	};

	template <typename T>
	std::vector<T> brezinskiTheta(const std::vector<T>& seq)
	{
		BrezinskiThetaMethod<T> btm(seq.size());
		return btm.next(seq);
	}


	template <typename T>
	class IteratedBrezinskiThetaMethod : public ExtrapolationMethod<T>
	{
	public:
		IteratedBrezinskiThetaMethod() { }
		IteratedBrezinskiThetaMethod(std::size_t seqLength)
		{
			_arj.reserve(seqLength);
		}

		virtual T next(const T& val)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf

			_arj.push_back(val);
			const std::size_t N = _arj.size() - 1;
			if (N < 3)
				return val;

			const std::size_t lMax = N / 3;
			std::size_t M = N;
			for (size_t l = 1; l <= lMax; ++l)
			{
				M -= 3;
				
				const T diff0 = _arj[M+1] - _arj[M];
				const T diff1 = _arj[M+2] - _arj[M+1];
				const T diff2 = _arj[M+3] - _arj[M+2];
				
				const T denom = diff2 * (diff1 - diff0) - diff0 * (diff2 - diff1);

				// Take care of under- or overflow
				if (abs(denom) <= std::numeric_limits<T>::min())
					_arj[M] = _arj[M+1] - diff0 * diff1 * (diff2 - diff1) / std::numeric_limits<T>::max();
				else
					_arj[M] = _arj[M+1] - diff0 * diff1 * (diff2 - diff1) / denom;
			}
			return _arj[N % 3];
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Iterated Brezinski's Theta Method"; }
		virtual std::size_t iterations() const { return _arj.size(); }

		virtual void reset()
		{
			_arj.clear();
		}

	protected:
		std::vector<T> _arj;
	};

	template <typename T>
	std::vector<T> iteratedBrezinskiTheta(const std::vector<T>& seq)
	{
		IteratedBrezinskiThetaMethod<T> ibt(seq.size());
		return ibt.next(seq);
	}


	// TODO: Check implementation of NevilleAitkenMethod
	template <typename T>
	class NevilleAitkenMethod : public ExtrapolationMethod<T>
	{
	public:
		NevilleAitkenMethod() { }
		NevilleAitkenMethod(std::size_t seqLength) { _diag.reserve(seqLength); }

		virtual T next(const T& val)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf

			_diag.push_back(val);

			// We need at least two elements here, so return if this is the first call
			if (_diag.size() == 1)
				return val;

			// Update counter diagonal from bottom to top
			for (std::size_t j = _diag.size()-1; j > 0; --j)
			{
				const std::size_t k = _diag.size() - j;
				const std::size_t n = j - 1;
				_diag[j-1] = (T(1) + T(n) / T(k)) * _diag[j] - T(n) / T(k) * _diag[j-1];
			}

			return _diag[0];
		}

		virtual std::vector<T> next(const std::vector<T>& vals) { return ExtrapolationMethod<T>::next(vals); }

		virtual std::size_t iterations() const { return _diag.size(); }
		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Neville Aitken Table"; }

		virtual void reset() { _diag.clear(); }

	protected:
		std::vector<T> _diag;
	};

	template <typename T>
	std::vector<T> nevilleAitken(const std::vector<T>& seq)
	{
		NevilleAitkenMethod<T> nam(seq.size());
		return nam.next(seq);
	}


	// TODO: Check implementation of RichardsonMethod
	template <typename T>
	class RichardsonMethod : public ExtrapolationMethod<T>
	{
	public:
		RichardsonMethod() { }
		RichardsonMethod(std::size_t seqLength)
		{
			_diag.reserve(seqLength);
			_pts.reserve(seqLength);
		}

		virtual T nextPoint(const T& val, const T& pt)
		{
			// Implementation of algorithm in
			// Ernst Joachim Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series,
			// Computer Physics Reports, Volume 10, Issues 5–6, December 1989, Pages 189-371, ISSN 0167-7977, http://dx.doi.org/10.1016/0167-7977(89)90011-7.
			// Also available here: http://arxiv.org/pdf/math/0306302v1.pdf

			_diag.push_back(val);
			_pts.push_back(pt);

			// We need at least two elements here, so return if this is the first call
			if (_diag.size() == 1)
				return val;

			const std::size_t n = _diag.size()-1;
			// Update counter diagonal from bottom to top
			for (std::size_t j = 1; j < _diag.size(); ++j)
			{
				const T denom = _pts[n-j] - pt;
				if (abs(denom) <= std::numeric_limits<T>::min())
					_diag[n-j] = (_diag[n-j+1] * _pts[n-j] - pt * _diag[n-j]) / std::numeric_limits<T>::max();
				else
					_diag[n-j] = (_diag[n-j+1] * _pts[n-j] - pt * _diag[n-j]) / denom;
			}

			return _diag[0];
		}

		virtual T next(const T& val)
		{
			return nextPoint(val, T(_diag.size()));
		}

		virtual std::vector<T> nextPoint(const std::vector<T>& vals, const std::vector<T>& pts)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = nextPoint(vals[i], pts[i]);

			return estimates;
		}

		virtual std::vector<T> next(const std::vector<T>& vals)
		{
			std::vector<T> estimates(vals.size(), T());
			for (std::size_t i = 0; i < estimates.size(); ++i)
				estimates[i] = nextPoint(vals[i], T(i + _diag.size()));

			return estimates;
		}

		virtual const char* const name() const { return compileTimeName(); }
		static const char* const compileTimeName() { return "Richardson Extrapolation"; }
		virtual std::size_t iterations() const { return _diag.size(); }

		virtual void reset()
		{
			_diag.clear();
			_pts.clear();
		}

	protected:
		std::vector<T> _diag;
		std::vector<T> _pts;
	};

	template <typename T>
	std::vector<T> richardsonPoints(const std::vector<T>& seq, const std::vector<T>& pts)
	{
		RichardsonMethod<T> rm(seq.size());
		return rm.nextPoint(seq, pts);
	}

	template <typename T>
	std::vector<T> richardson(const std::vector<T>& seq)
	{
		RichardsonMethod<T> rm(seq.size());
		return rm.next(seq);
	}

}

#endif
