// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2017: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_CLIERRORPARSER_HPP_
#define CASEMA_CLIERRORPARSER_HPP_

#include "ErrorEstimate.hpp"

namespace casema
{
	template <typename real_t, template <class T> class Opts_t, template <class T> class Model_t>
	bool processErrorOptions(Opts_t<real_t>& opts, const Model_t<real_t>& model, const real_t& T)
	{
		// If summands and abscissa are given, use them
		if ((opts.sigma > 0) && (opts.summands > 0))
		{
			opts.errorCons = casema::consistencyError(opts.sigma, T, model);
			opts.errorTrunc = casema::truncationError(opts.sigma, opts.summands, T, model);
			opts.error = opts.errorCons + opts.errorTrunc;
			opts.errorWeight = 0.5;
		}
		else if ((opts.errorWeight > 0) && (opts.errorWeight < 1) && (opts.error > 0))
		{
			if (opts.sigma > 0)
				casema::paramsFromError(opts.sigma, T, model, opts.errorWeight, opts.error, opts.sigma, opts.summands);
			else
				casema::paramsFromError(mpfr::mpreal(0), T, model, opts.errorWeight, opts.error, opts.sigma, opts.summands);

			opts.errorCons = casema::consistencyError(opts.sigma, T, model);
			opts.errorTrunc = casema::truncationError(opts.sigma, opts.summands, T, model);
		}
		else
			return false;

		return true;
	}


	template <typename real_t, template <class T> class Opts_t, template <class T> class Model_t>
	bool processErrorOptions(Opts_t<real_t>& opts, const Model_t<real_t>& model)
	{
		// If summands and abscissa are given, use them
		if ((opts.sigma > 0) && (opts.summands > 0))
		{
			opts.errorCons = casema::consistencyError(opts.sigma, model);
			opts.errorTrunc = casema::truncationError(opts.sigma, opts.summands, model);
			opts.error = opts.errorCons + opts.errorTrunc;
			opts.errorWeight = 0.5;
		}
		else if ((opts.errorWeight > 0) && (opts.errorWeight < 1) && (opts.error > 0))
		{
			if (opts.sigma > 0)
				casema::paramsFromError(opts.sigma, model, opts.errorWeight, opts.error, opts.sigma, opts.summands);
			else
				casema::paramsFromError(mpfr::mpreal(0), model, opts.errorWeight, opts.error, opts.sigma, opts.summands);

			opts.errorCons = casema::consistencyError(opts.sigma, model);
			opts.errorTrunc = casema::truncationError(opts.sigma, opts.summands, model);
		}
		else
			return false;

		return true;
	}
}

#endif
