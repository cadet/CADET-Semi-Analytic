// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2019: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_CLIEXTRAPOLATIONPARSER_HPP_
#define CASEMA_CLIEXTRAPOLATIONPARSER_HPP_

#include <string>
#include <iostream>
#include <vector>

#include "StringUtil.hpp"

#include "Extrapolation.hpp"
#include "ExtrapolationHelpers.hpp"

namespace casema
{
	template <typename real_t, template <class> class Options>
	casema::ConvergenceMonitoredExtrapolationMethod<real_t>* createExtrapolationMethod(const std::string& name, const Options<real_t>& opts)
	{
		std::string lowName = name;
		util::toLower(lowName);

		casema::ConvergenceMonitoredExtrapolationMethod<real_t>* convMon = nullptr;

		if (lowName == "ide")
			convMon = new casema::ConvergenceMonitor<real_t, casema::IdentityExtrapolation>(casema::IdentityExtrapolation<real_t>());
		else if (lowName == "ads")
			convMon = new casema::ConvergenceMonitor<real_t, casema::AitkenDeltaSquaredMethod>(casema::AitkenDeltaSquaredMethod<real_t>());
		else if (lowName == "wem")
			convMon = new casema::ConvergenceMonitor<real_t, casema::WynnEpsilonMethod>(casema::WynnEpsilonMethod<real_t>());
		else if (lowName == "wrm")
			convMon = new casema::ConvergenceMonitor<real_t, casema::WynnRhoMethod>(casema::WynnRhoMethod<real_t>());
		else if (lowName == "iad")
			convMon = new casema::ConvergenceMonitor<real_t, casema::IteratedAitkenDeltaSquaredMethod>(casema::IteratedAitkenDeltaSquaredMethod<real_t>());
		else if (lowName == "lum")
			convMon = new casema::ConvergenceMonitor<real_t, casema::LevinUMethod>(casema::LevinUMethod<real_t>());
		else if (lowName == "ltm")
			convMon = new casema::ConvergenceMonitor<real_t, casema::LevinTMethod>(casema::LevinTMethod<real_t>());
		else if (lowName == "ibt")
			convMon = new casema::ConvergenceMonitor<real_t, casema::IteratedBrezinskiThetaMethod>(casema::IteratedBrezinskiThetaMethod<real_t>());
		else if (lowName == "btm")
			convMon = new casema::ConvergenceMonitor<real_t, casema::BrezinskiThetaMethod>(casema::BrezinskiThetaMethod<real_t>());
		else if (lowName == "nam")
			convMon = new casema::ConvergenceMonitor<real_t, casema::NevilleAitkenMethod>(casema::NevilleAitkenMethod<real_t>());
		else if (lowName == "rem")
			convMon = new casema::ConvergenceMonitor<real_t, casema::RichardsonMethod>(casema::RichardsonMethod<real_t>());

		if (convMon)
		{
			// Configure
			convMon->times(opts.convTimes);
			convMon->thresholdAbs(opts.convThresholdAbs);
			convMon->thresholdRel(opts.convThresholdRel);
		}

		return convMon;
	}


	template <typename real_t, template <class> class Options>
	bool addExtrapolationMethods(std::ostream* os, const std::string& emList, casema::ConsensusEstimator<real_t>& est, const Options<real_t>& opts)
	{
		std::vector<std::string> methods = casema::util::split(emList, ',');
		for (std::vector<std::string>::iterator it = methods.begin(); it != methods.end(); ++it)
		{
			casema::ConvergenceMonitoredExtrapolationMethod<real_t>* ptrMethod = createExtrapolationMethod<real_t>(*it, opts);
			if (ptrMethod)
				est.add(ptrMethod);
			else
			{
				if (*os)
					(*os) << "WARNING: Ignoring unknown extrapolation method \"" << *it << "\"" << std::endl;
			}
		}

		bool retVal = est.numEstimators() > 0;
		if (!retVal)
		{
			casema::ConvergenceMonitoredExtrapolationMethod<real_t>* convMon = new casema::ConvergenceMonitor<real_t, casema::IdentityExtrapolation>(casema::IdentityExtrapolation<real_t>());
			convMon->times(opts.convTimes);
			convMon->thresholdAbs(opts.convThresholdAbs);
			convMon->thresholdRel(opts.convThresholdRel);

			est.add(convMon);
		}

		// Configure ConsensusEstimator
		est.mustAgree(opts.consAgree);
		return retVal;
	}


	template <typename real_t, template <class> class Options>
	bool addExtrapolationMethods(const std::string& emList, casema::ConsensusEstimator<real_t>& est, const Options<real_t>& opts)
	{
		return addExtrapolationMethods(&std::cout, emList, est, opts);
	}


	template <typename real_t, template <class> class Options>
	bool addExtrapolationMethods(const std::string& emList, casema::ConsensusEstimator<real_t>& est, const Options<real_t>& opts, bool silent)
	{
		if (silent)
			return addExtrapolationMethods(nullptr, emList, est, opts);
		else
			return addExtrapolationMethods(&std::cout, emList, est, opts);
	}
}

#endif
